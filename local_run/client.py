import argparse
import socket
import numpy as np
import os
import psutil
import time
from scipy.sparse import load_npz, diags, csr_matrix
import GWAS_lib
import communication
import random
from concurrent.futures import ThreadPoolExecutor
import utilities
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)
np.seterr(divide='ignore', invalid='ignore', over='ignore', under='ignore')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--number_of_parties", type=int, required=True)
    parser.add_argument("--party_id", type=int, required=True)
    parser.add_argument("--base_port", type=int, required=True)
    parser.add_argument("--number_of_samples", type=int, required=True)
    parser.add_argument("--number_of_snps", type=int, required=True)
    parser.add_argument("--number_of_covariates", type=int, required=True)
    parser.add_argument("--number_of_blocks", type=int, required=True)
    parser.add_argument("--number_of_folds", type=int, required=True)
    parser.add_argument("--number_of_blocks_per_run", type=int, required=True)
    args = parser.parse_args()
    B, K, N, M, P, p, C = args.number_of_blocks, args.number_of_folds, args.number_of_samples, args.number_of_snps, args.number_of_parties, args.party_id, args.number_of_covariates
    bulk=args.number_of_blocks_per_run
    loader_times = 0
    while not os.path.exists(f'Data/server_ready_{p}.txt'):
        time.sleep(0.1)
    print('SERVER IS READY')
    with open('Data/ip_address_file.txt', 'r') as file:
        server_hostname = file.read().strip()
    party_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    party_socket.setsockopt(socket.SOL_SOCKET, socket.SO_KEEPALIVE, 1)
    print(f'party {p} connecting via {args.base_port + p}')
    party_socket.connect((server_hostname, args.base_port + p))
    party_start_time = time.time()
    reduction_count=0
    total_comm_size=0
    np.random.seed(0)
    tt = time.time()
    y = np.load(
        'Data/N{}_M{}_C{}_P{}_B{}/Party_{}/y.npy'.format(N, M, C, P, B, p))
    loader_times += (time.time() - tt)
    v_6, sum_6 = GWAS_lib.standardization_normalizers(1, P, p, 6)
    v_7, sum_7 = GWAS_lib.standardization_normalizers(1, P, p, 7)
    masked_y_sum = np.sum(y, axis=0) + v_6
    red_time=time.time()
    total_comm_size+=utilities.get_size_in_gb(masked_y_sum)
    reduction_count+=time.time()-red_time
    communication.send_data_to_server(party_socket, masked_y_sum, p)
    aggregated_y = communication.receive_data_from_server(party_socket)
    red_time = time.time()
    total_comm_size += utilities.get_size_in_gb(aggregated_y)
    reduction_count += time.time() - red_time
    y_mean = (aggregated_y - sum_6) / N
    std_y = (np.sum(np.square(y - y_mean)) / (N - 1))
    std_y = std_y + v_7
    red_time = time.time()
    total_comm_size += utilities.get_size_in_gb(std_y)
    reduction_count += time.time() - red_time
    communication.send_data_to_server(party_socket, std_y, p)
    agg_std_y = communication.receive_data_from_server(party_socket)
    red_time = time.time()
    total_comm_size += utilities.get_size_in_gb(agg_std_y)
    reduction_count += time.time() - red_time
    agg_std_y = agg_std_y.reshape(1, -1)
    agg_std_y = np.sqrt(agg_std_y - sum_7)
    S_y = 1 / agg_std_y
    GWAS_lib.generate_Z_mask(M, B, p, P, N, C)

    start_index = (p - 1) * int(N / P)
    end_index = start_index + int(N / P)
    end_index = start_index + (N // P if p < P  else N - (N // P) * (P - 1))
    tt = time.time()
    Z = np.load(
        'Data/N{}_M{}_C{}_P{}_B{}/Party_{}/Z.npy'.format(N, M, C, P, B, p))
    loader_times += (time.time() - tt)
    Z = np.hstack((np.ones((Z.shape[0], 1)), Z))
    file_loaded = False
    max_attempts = 1000
    attempts = 0
    delay = 0.025
    Z_mask = None
    while not file_loaded and attempts < max_attempts:
        try:
            Z_mask = load_npz('Data/N{}_M{}_C{}_P{}_B{}/Masks/Z_mask.npz'.format(N, M, C, P,B))
            file_loaded = True
        except Exception as e:
            attempts += 1
            time.sleep(delay)

    if not file_loaded:
        print("Failed to load the file after multiple attempts. Increase time for Z_mask.")

    Z_mask_party = Z_mask[:, start_index:end_index]
    Z_masked = Z_mask_party @ Z
    red_time = time.time()
    total_comm_size += utilities.get_size_in_gb(Z_masked)
    reduction_count += time.time() - red_time
    communication.send_data_to_server(party_socket, Z_masked, p)

    k_y = GWAS_lib.generate_a_number(0)
    masked_y = Z_mask_party @ y
    masked_y = k_y * masked_y
    red_time = time.time()
    total_comm_size += utilities.get_size_in_gb(masked_y)
    reduction_count += time.time() - red_time
    communication.send_data_to_server(party_socket, masked_y, p)
    Y_tilde = communication.receive_data_from_server(party_socket)
    red_time = time.time()
    total_comm_size += utilities.get_size_in_gb(Y_tilde)
    reduction_count += time.time() - red_time
    Y_tilde = (Z_mask.transpose() @ Y_tilde @ S_y)[start_index:end_index]
    Y_tilde = ((1 / k_y) * Y_tilde).astype(np.float32)
    del Z_masked, masked_y, S_y, y

    GWAS_lib.generate_O(M, B, K, p, P, N, C)

    file_loaded = False
    max_attempts = 1000
    attempts = 0
    delay = 0.025
    O = 0
    while not file_loaded and attempts < max_attempts:
        try:
            O = load_npz(
                'Data/N{}_M{}_C{}_P{}_B{}/Masks/O.npz'.format(N, M, C, P,B))
            file_loaded = True
        except Exception as e:
            attempts += 1
            time.sleep(delay)

    if not file_loaded:
        print("Failed to load the file after multiple attempts. Increase time for O.")

    O2 = O[:, start_index:end_index]
    k_1 = GWAS_lib.generate_a_number(1)
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

    indices = np.arange(Y_tilde.shape[0])
    np.random.seed(int(time.time()))

    np.random.shuffle(indices)

    chunks = np.array_split(indices, K)

    for k in range(K):
        Y_tilde_copy = Y_tilde.copy()

        Y_tilde_copy[chunks[k], :] = 0
        O2 = O2.astype(np.float32)
        Y_tilde_copy = Y_tilde_copy.astype(np.float32)
        masked_y_tilde = (O2@Y_tilde_copy)
        masked_y_tilde *= k_1
        red_time = time.time()
        total_comm_size += utilities.get_size_in_gb(masked_y_tilde)
        reduction_count += time.time() - red_time
        communication.send_data_to_server(party_socket, masked_y_tilde, p)

        Y_tilde_copy = Y_tilde.copy()

        mask = np.ones(Y_tilde.shape[0], bool)
        mask[chunks[k]] = 0
        Y_tilde_copy[mask, :] = 0
        O2 = O2.astype(np.float32)
        Y_tilde_copy = Y_tilde_copy.astype(np.float32)
        masked_y_tilde = (O2@Y_tilde_copy)
        masked_y_tilde *= k_1
        red_time = time.time()
        total_comm_size += utilities.get_size_in_gb(masked_y_tilde)
        reduction_count += time.time() - red_time
        communication.send_data_to_server(party_socket, masked_y_tilde, p)

    del Y_tilde, Y_tilde_copy, masked_y_tilde
    extra_elements = [np.random.randint(1, 100) for _ in range(B)]
    extra_elements[B-1]+= M%(int(M/B))
    B_blocked = GWAS_lib.split_number(B, bulk)
    num_loops = len(B_blocked)
    initial_block_size = int(M / B)
    for n in range(num_loops):

        #####################################################################################################
        #------------------------------------------QUALITY CONTROL------------------------------------------#
        #####################################################################################################

        if n>0:
            while not os.path.exists(f'Data/server_ready_loop_{n}.txt'):
                time.sleep(0.1)
        tt = time.time()
        X_list = [None] * len(B_blocked[n])
        with ThreadPoolExecutor(max_workers=16) as executor:
            future_to_idx = {executor.submit(GWAS_lib.load_X, N, M, C, P, B, p, b, idx, X_list): idx for idx, b in
                             enumerate(B_blocked[n])}
        tt = time.time()
        utilities.print_memory_usage()

        ones_list = []
        twos_list = []
        sums_list = []

        for X in X_list:
            ones = np.array(
                (X == 1).sum(axis=0)).ravel()
            twos = np.array(
                (X == 2).sum(axis=0)).ravel()
            sums = np.array(X.sum(axis=0)).ravel()

            ones_list.append(ones)
            twos_list.append(twos)
            sums_list.append(sums)

        ones = np.hstack(ones_list)
        twos = np.hstack(twos_list)
        sums = np.hstack(sums_list)
        del ones_list, twos_list, sums_list
        num_new_elements = random.randint(1, 100)
        rand_indices = np.random.randint(0, len(ones), size=num_new_elements)

        ones_random = ones[rand_indices] + np.random.choice([-1, 1], size=num_new_elements)
        twos_random = twos[rand_indices] + np.random.choice([-1, 1], size=num_new_elements)
        sums_random = sums[rand_indices] + np.random.choice([-1, 1], size=num_new_elements)


        ones = np.append(ones, ones_random)
        twos = np.append(twos, twos_random)
        sums = np.append(sums, sums_random)

        del ones_random, twos_random, sums_random
        v_2, sum_2 = GWAS_lib.standardization_normalizers(len(ones), P, p, 2)
        v_3, sum_3 = GWAS_lib.standardization_normalizers(len(ones), P, p, 3)
        v_4, sum_4 = GWAS_lib.standardization_normalizers(len(ones), P, p, 4)
        ones = ones + v_2
        twos = twos + v_3
        sums = sums + v_4
        red_time = time.time()
        total_comm_size += utilities.get_size_in_gb(ones)
        reduction_count += time.time() - red_time
        communication.send_data_to_server(party_socket, ones, p)
        del ones, v_2, v_3, v_4
        aggregated_vector = communication.receive_data_from_server(party_socket)
        red_time = time.time()
        total_comm_size += utilities.get_size_in_gb(aggregated_vector)
        reduction_count += time.time() - red_time
        ones_count = (aggregated_vector - sum_2)[:-num_new_elements]
        red_time = time.time()
        total_comm_size += utilities.get_size_in_gb(twos)
        reduction_count += time.time() - red_time
        communication.send_data_to_server(party_socket, twos, p)
        del twos
        aggregated_vector = communication.receive_data_from_server(party_socket)
        red_time = time.time()
        total_comm_size += utilities.get_size_in_gb(aggregated_vector)
        reduction_count += time.time() - red_time
        twos_count = (aggregated_vector - sum_3)[:-num_new_elements]
        red_time = time.time()
        total_comm_size += utilities.get_size_in_gb(sums)
        reduction_count += time.time() - red_time
        communication.send_data_to_server(party_socket, sums, p)
        del sums
        aggregated_vector = communication.receive_data_from_server(party_socket)
        red_time = time.time()
        total_comm_size += utilities.get_size_in_gb(aggregated_vector)
        reduction_count += time.time() - red_time
        means = ((aggregated_vector - sum_4) / N)[:-num_new_elements]
        removable_column_indices = []
        maf = ((twos_count * 2) + ones_count) / (2 * N)
        obs_counts_array = np.vstack([(N - ones_count - twos_count), ones_count, twos_count]).T
        chi2_values = np.apply_along_axis(GWAS_lib.hwe_chi_square, 1, obs_counts_array)
        removable_column_indices.extend(np.where((maf < 0.001) & (chi2_values > 23928))[0].tolist())
        tt = time.time()
        start_col = 0
        stds = np.zeros(sum(X.shape[1] for X in X_list))
        for X in X_list:
            end_col = start_col + X.shape[1]
            local_means = means[start_col:end_col]
            sq_diffs = csr_matrix(X - local_means)
            sum_sq = np.array(sq_diffs.power(2).sum(axis=0)).ravel()
            stds[start_col:end_col] += sum_sq
            start_col = end_col
        stds = (stds / (N - 1))
        del start_col, end_col, local_means, sq_diffs, sum_sq
        stds_random = stds[rand_indices] + np.random.choice([-0.1, 0.1], size=num_new_elements)
        stds = np.append(stds, stds_random)
        v_5, sum_5 = GWAS_lib.standardization_normalizers(len(stds), P, p, 5)
        stds = stds + v_5
        red_time = time.time()
        total_comm_size += utilities.get_size_in_gb(stds)
        reduction_count += time.time() - red_time
        communication.send_data_to_server(party_socket, stds, p)
        del stds
        aggregated_vector = communication.receive_data_from_server(party_socket)
        red_time = time.time()
        total_comm_size += utilities.get_size_in_gb(aggregated_vector)
        reduction_count += time.time() - red_time

        stds = np.sqrt(np.maximum(aggregated_vector - sum_5, 0))[:-num_new_elements]

        del aggregated_vector
        removed_so_far = 0

        global_col_counter = 0

        new_stds = []
        kept_indices = [j for j in range(len(stds)) if j not in removable_column_indices]

        #####################################################################################################
        #-----------------------DISTRIBUTED PROJECTION OF COVARIATES AND STANDARDIZING----------------------#
        #####################################################################################################

        stds = stds[kept_indices].astype(np.float32)

        del kept_indices
        for i, X in enumerate(X_list):
            num_cols = X.shape[1]
            valid_removable_column_indices = {global_idx - global_col_counter for global_idx in removable_column_indices
                                              if global_col_counter <= global_idx < global_col_counter + num_cols}
            keep_indices = [j for j in range(num_cols) if j not in valid_removable_column_indices]
            X_list[i] = X[:, keep_indices]
            if (i+1)%P==0:
                GWAS_lib.generate_O_X(M, B, p, P, N, C, X_list[i-p+1].shape[1],extra_elements[(i-p+1) + (n * bulk)] + initial_block_size -X_list[i-p+1].shape[1],(i-p+1) + (n * bulk))
            if i ==len(X_list)-1:
                if (i+1)%P!=0:
                    temp =(i+1)%P
                    if p <= temp:
                        GWAS_lib.generate_O_X(M, B, p, P, N, C, X_list[i-p+1].shape[1], extra_elements[(i - p + 1) + (n * bulk)] + initial_block_size - X_list[i-p+1].shape[1], (i - p + 1) + (n * bulk))
            global_col_counter += num_cols
            X_list[i] = X_list[i].astype(np.float32)
            if GWAS_lib.compute_sparsity(X_list[i]) < 0.7:
                X_list[i] = X_list[i].toarray()
        start_col=0
        S=[]
        for X in X_list:
            end_col=start_col+X.shape[1]
            s=stds[start_col:end_col]
            S.append(diags(1/s).tocsr())
            start_col=end_col
        O_X = [0 for _ in range(len(X_list))]

        for _ in range(len(X_list)):
            file_loaded = False
            max_attempts = 1000
            attempts = 0
            delay = 0.025
            while not file_loaded and attempts < max_attempts:
                try:
                    O_X[_] = load_npz('Data/N{}_M{}_C{}_P{}_B{}/Masks/O_X_block_{}.npz'.format(N, M, C, P, B, (_ + (n * bulk)) + 1))
                    file_loaded = True
                except Exception as e:
                    attempts += 1
                    time.sleep(delay)
            t2 = time.time()
            if not file_loaded:
                print(f"Failed to load the file after multiple attempts. Increase time for O_X.")
        current_count = len(X_list)
        right = [[] for _ in range(current_count)]
        masked_X = [[] for _ in range(current_count)]
        for i, X in enumerate(X_list):
            X = X.astype(np.float32)
            O_X[i] = O_X[i].astype(np.float32)

            right[i] = (X@O_X[i].T)
            Z_mask_party=Z_mask_party.astype(np.float32)

            masked_X[i] = (Z_mask_party@right[i])
        for _ in range(current_count):
            red_time = time.time()
            total_comm_size += utilities.get_size_in_gb(masked_X[_])
            reduction_count += time.time() - red_time
            communication.send_data_to_server(party_socket, masked_X[_], p)
            a=communication.receive_data_from_server(party_socket)
        del masked_X
        X_tilde = [[] for _ in range(current_count)]
        for _ in range(current_count):
            X = communication.receive_data_from_server(party_socket)
            red_time = time.time()
            total_comm_size += utilities.get_size_in_gb(X)
            reduction_count += time.time() - red_time
            S[_] = S[_].astype(np.float32)
            right = (O_X[_]@S[_])
            X = X.astype(np.float32)
            right = right.astype(np.float32)
            intermediate = (X@right)
            del right
            Z_mask=Z_mask.astype(np.float32)
            X_tilde[_] = (Z_mask.T@intermediate)[start_index:end_index, :]
            del intermediate

            #####################################################################################################
            #--------------------------------------LEVEL 1 RIDGE REGRESSION-------------------------------------#
            #####################################################################################################

            if (_ +1) % P == 0:
                GWAS_lib.generate_O_b(M, B, p, P, N, C, X_tilde[_-p+1].shape[1], int((_-p+1) + (n * bulk)))
            if _==current_count-1:
                if (_+1)%P != 0:
                    temp = (_+1)%P
                    if p<= temp:
                        GWAS_lib.generate_O_b(M,B,p,P,N,C,X_tilde[_-p+1].shape[1],(_-p+1)+(n*bulk))
        O_b=[[] for _ in range(current_count)]
        sending=[[] for _ in range(current_count)]
        for _ in range(current_count):
            file_loaded = False
            max_attempts = 1000
            attempts = 0
            delay = 0.025
            while not file_loaded and attempts < max_attempts:
                try:
                    O_b[_] = load_npz('Data/N{}_M{}_C{}_P{}_B{}/Masks/O_b_block_{}.npz'.format(N, M, C, P, B, (_ + (n * bulk)) + 1))
                    file_loaded = True
                except Exception as e:
                    attempts += 1
                    time.sleep(delay)
        for _ in range(current_count):
            X_tilde[_]=X_tilde[_].astype(np.float32)
            O_b[_]=O_b[_].astype(np.float32)
            intermediate = (X_tilde[_]@ O_b[_])
            sending[_] = (O2@ intermediate)
        for _ in range(current_count):
            red_time = time.time()
            total_comm_size += utilities.get_size_in_gb(sending[_])
            reduction_count += time.time() - red_time
            communication.send_data_to_server(party_socket, sending[_], p)

        #####################################################################################################
        #--------------------------------------LEVEL 2 RIDGE REGRESSION-------------------------------------#
        #####################################################################################################

        for _ in range(current_count):
            for k in range(K):
                X_tilde_k = X_tilde[_].copy()
                idx = chunks[k]
                X_tilde_k[idx, :] = 0
                X_tilde_k=X_tilde_k.astype(np.float32)
                O_b[_]=O_b[_].astype(np.float32)
                intermediate = (X_tilde_k@ O_b[_])
                X_tilde_k = (O2@ intermediate)
                red_time = time.time()
                total_comm_size += utilities.get_size_in_gb(X_tilde_k)
                reduction_count += time.time() - red_time
                communication.send_data_to_server(party_socket, X_tilde_k, p)
                del X_tilde_k, intermediate
        D=[[] for _ in range(current_count)]
        result_matrix=[[] for _ in range(current_count)]

        #####################################################################################################
        #-----------------------------DISTRIBUTED SINGLE SNP ASSOCIATION TESTING----------------------------#
        #####################################################################################################


        for _ in range(current_count):
            shape1 = X_tilde[_].shape[1]
            diagonal_entries = [GWAS_lib.generate_a_number(i + int(shape1 * ((n*bulk)+_))) for i in range(shape1)]
            diagonal_entries = np.array(diagonal_entries)
            diagonal_entries = diagonal_entries.astype(np.float32)
            D = diags(diagonal_entries).tocsr()
            intermediate = (O2@ X_tilde[_])
            D = D.astype(np.float32)
            result_matrix= (intermediate@D)
            red_time = time.time()
            total_comm_size += utilities.get_size_in_gb(result_matrix)
            reduction_count += time.time() - red_time
            communication.send_data_to_server(party_socket, result_matrix, p)
        del result_matrix

    print(f'TIME TAKEN FOR {B} blocks is {time.time()-party_start_time-reduction_count}')
    #print(f'reduce from server {reduction_count}')
    print(f'SIZE OF DATA IN TOTAL FOR {B} blocks is {total_comm_size}')
    party_socket.close()

if __name__ == "__main__":
    main()
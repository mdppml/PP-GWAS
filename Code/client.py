import shutil
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
from sparse_dot_mkl import dot_product_mkl
import utilities
os.environ["MKL_NUM_THREADS"] = "16"
os.environ["MKL_DYNAMIC"] = "TRUE"

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
    while not os.path.exists(f'../test_site/Data/server_ready_{p}.txt'):
        time.sleep(0.1)
    print('SERVER IS READY')
    with open('../test_site/Data/ip_address_file.txt', 'r') as file:
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
        '../test_site/Data/N{}_M{}_C{}_P{}_B{}/Party_{}/y.npy'.format(N, M, C, P, B, p))
    loader_times += (time.time() - tt)
    r_y, sum_r_y = GWAS_lib.randomized_encoding_addtion(1, P, p, 6)
    r_y_std, sum_r_y_std = GWAS_lib.randomized_encoding_addtion(1, P, p, 7)
    masked_y_sum = np.sum(y, axis=0) + r_y
    red_time=time.time()
    total_comm_size+=utilities.get_size_in_gb(masked_y_sum)
    reduction_count+=time.time()-red_time
    communication.send_data_to_server(party_socket, masked_y_sum, p)
    aggregated_y = communication.receive_data_from_server(party_socket)
    red_time = time.time()
    total_comm_size += utilities.get_size_in_gb(aggregated_y)
    reduction_count += time.time() - red_time
    y_mean = (aggregated_y - sum_r_y) / N
    std_y = (np.sum(np.square(y - y_mean)) / (N - 1))
    std_y = std_y + r_y_std
    red_time = time.time()
    total_comm_size += utilities.get_size_in_gb(std_y)
    reduction_count += time.time() - red_time
    communication.send_data_to_server(party_socket, std_y, p)
    agg_std_y = communication.receive_data_from_server(party_socket)
    red_time = time.time()
    total_comm_size += utilities.get_size_in_gb(agg_std_y)
    reduction_count += time.time() - red_time

    agg_std_y = agg_std_y.reshape(1, -1)
    agg_std_y = np.sqrt(agg_std_y - sum_r_y_std)
    S_y = 1 / agg_std_y
    S_y_scalar = float(S_y)
    GWAS_lib.generate_O_Z(M, B, p, P, N, C)
    GWAS_lib.generate_O_Z_prime(M, B, p, P, N, C)

    start_index = (p - 1) * int(N / P)
    end_index = start_index + int(N / P)
    end_index = start_index + (N // P if p < P  else N - (N // P) * (P - 1))
    tt = time.time()
    Z = np.load(
        '../test_site/Data/N{}_M{}_C{}_P{}_B{}/Party_{}/Z.npy'.format(N, M, C, P, B, p))
    loader_times += (time.time() - tt)
    Z = np.hstack((np.ones((Z.shape[0], 1)), Z))
    file_loaded = False
    max_attempts = 1000
    attempts = 0
    delay = 0.025
    O_Z = None
    while not file_loaded and attempts < max_attempts:
        try:
            O_Z = load_npz(
                '../test_site/Data/N{}_M{}_C{}_P{}_B{}/Party_{}/Masks/O_Z.npz'.format(N, M, C, P, B, p))
            file_loaded = True
        except Exception as e:
            attempts += 1
            time.sleep(delay)

    if not file_loaded:
        print("Failed to load the file after multiple attempts. Increase time for O_Z.")

    file_loaded = False
    attempts = 0

    O_Z_prime = None
    while not file_loaded and attempts < max_attempts:
        try:
            O_Z_prime = load_npz(
                '../test_site/Data/N{}_M{}_C{}_P{}_B{}/Party_{}/Masks/O_Z_prime.npz'.format(N, M, C, P, B, p))
            file_loaded = True
        except Exception as e:
            attempts += 1
            time.sleep(delay)

    if not file_loaded:
        print("Failed to load the file after multiple attempts. Increase time for O_Z_prime.")

    Q = O_Z.toarray()  # or use .A

    QtQ = Q.T @ Q
    diag_diff = np.abs(np.diag(QtQ) - 1)
    off_diag = np.abs(QtQ - np.diag(np.diag(QtQ)))
    off_diag_max = np.max(off_diag)

    print("max |diag(QᵀQ) − 1| =", diag_diff.max())
    print("max |off-diag(QᵀQ)| =", off_diag_max)


    O_Z_p = O_Z[:, start_index:end_index]
    Z_masked = O_Z_p @ Z
    Z_masked = dot_product_mkl(Z_masked, O_Z_prime.transpose())
    red_time = time.time()
    total_comm_size += utilities.get_size_in_gb(Z_masked)
    reduction_count += time.time() - red_time
    communication.send_data_to_server(party_socket, Z_masked, p)

    num_dummy_cols = 10
    U = 1 + num_dummy_cols

    np.random.seed(12345)
    M_y = np.random.randn(N, num_dummy_cols).astype(np.float32)
    M_y_p = M_y[start_index:end_index, :]

    Y_concat = np.hstack((y, M_y_p))

    np.random.seed(54321)
    rho = np.random.permutation(U)
    y_col_index = int(np.where(rho == 0)[0][0])
    Y_permuted = Y_concat[:, rho]

    O_y = GWAS_lib.get_O(U, 1, 999)
    O_y = O_y.astype(np.float32)

    k_y = GWAS_lib.generate_a_number(0)

    O_Z_p = O_Z_p.astype(np.float32)
    Y_permuted = Y_permuted.astype(np.float32)

    masked_Y = dot_product_mkl(O_Z_p, Y_permuted)
    masked_Y = dot_product_mkl(masked_Y, O_y.transpose())
    masked_Y = k_y * masked_Y

    red_time = time.time()
    total_comm_size += utilities.get_size_in_gb(masked_Y)
    reduction_count += time.time() - red_time
    communication.send_data_to_server(party_socket, masked_Y, p)

    Y_tilde = communication.receive_data_from_server(party_socket)
    red_time = time.time()
    total_comm_size += utilities.get_size_in_gb(Y_tilde)
    reduction_count += time.time() - red_time

    O_Z = O_Z.astype(np.float32)
    Y_tilde = Y_tilde.astype(np.float32)

    left = dot_product_mkl(O_Z.transpose(), Y_tilde)
    O_y = O_y.astype(np.float32)
    Y_full = dot_product_mkl(left, O_y)

    Y_full = (1.0 / k_y) * Y_full
    Y_full = Y_full * S_y_scalar

    Y_tilde = Y_full[start_index:end_index, y_col_index:y_col_index + 1].astype(np.float32)

    del Z_masked, M_y, M_y_p, Y_concat, Y_permuted, Y_full, S_y, y

    GWAS_lib.generate_O_y_tilde(M, B, K, p, P, N, C)

    file_loaded = False
    max_attempts = 1000
    attempts = 0
    delay = 0.025
    O_y_tilde = 0
    while not file_loaded and attempts < max_attempts:
        try:
            O_y_tilde = load_npz(
                '../test_site/Data/N{}_M{}_C{}_P{}_B{}/Party_{}/Masks/O_y_tilde.npz'.format(N, M, C, P, B, p))
            file_loaded = True
        except Exception as e:
            attempts += 1
            time.sleep(delay)

    if not file_loaded:
        print("Failed to load the file after multiple attempts. Increase time for O.")

    O_y_tilde_p = O_y_tilde[:, start_index:end_index]
    k_1 = GWAS_lib.generate_a_number(1)
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

    indices = np.arange(Y_tilde.shape[0])
    np.random.seed(int(time.time()))

    np.random.shuffle(indices)

    chunks = np.array_split(indices, K)

    for k in range(K):
        Y_tilde_copy = Y_tilde.copy()

        Y_tilde_copy[chunks[k], :] = 0
        O_y_tilde_p = O_y_tilde_p.astype(np.float32)
        Y_tilde_copy = Y_tilde_copy.astype(np.float32)
        masked_y_tilde = dot_product_mkl(O_y_tilde_p, Y_tilde_copy)
        masked_y_tilde *= k_1
        red_time = time.time()
        total_comm_size += utilities.get_size_in_gb(masked_y_tilde)
        reduction_count += time.time() - red_time
        communication.send_data_to_server(party_socket, masked_y_tilde, p)

        Y_tilde_copy = Y_tilde.copy()

        mask = np.ones(Y_tilde.shape[0], bool)
        mask[chunks[k]] = 0
        Y_tilde_copy[mask, :] = 0
        O_y_tilde_p = O_y_tilde_p.astype(np.float32)
        Y_tilde_copy = Y_tilde_copy.astype(np.float32)
        masked_y_tilde = dot_product_mkl(O_y_tilde_p, Y_tilde_copy)
        masked_y_tilde *= k_1
        red_time = time.time()
        total_comm_size += utilities.get_size_in_gb(masked_y_tilde)
        reduction_count += time.time() - red_time
        communication.send_data_to_server(party_socket, masked_y_tilde, p)

    del Y_tilde, Y_tilde_copy, masked_y_tilde


    np.random.seed(1234)
    extra_elements = [np.random.randint(1, 100) for _ in range(B)]
    extra_elements[B - 1] += M % int(M / B)
    B_blocked = GWAS_lib.split_number(B, bulk)
    num_loops = len(B_blocked)
    initial_block_size = int(M / B)

    for n in range(num_loops):

        #####################################################################################################
        #------------------------------------------QUALITY CONTROL------------------------------------------#
        #####################################################################################################

        if n>0:
            while not os.path.exists(f'../test_site/Data/server_ready_loop_{n}.txt'):
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
        r_one_freq, sum_r_one_freq = GWAS_lib.randomized_encoding_addtion(len(ones), P, p, 2)
        r_two_freq, sum_r_two_freq = GWAS_lib.randomized_encoding_addtion(len(ones), P, p, 3)
        r_total_freq, sum_r_total_freq = GWAS_lib.randomized_encoding_addtion(len(ones), P, p, 4)
        ones = ones + r_one_freq
        twos = twos + r_two_freq
        sums = sums + r_total_freq
        red_time = time.time()
        total_comm_size += utilities.get_size_in_gb(ones)
        reduction_count += time.time() - red_time
        communication.send_data_to_server(party_socket, ones, p)
        del ones, r_one_freq, r_two_freq, r_total_freq
        aggregated_vector = communication.receive_data_from_server(party_socket)
        red_time = time.time()
        total_comm_size += utilities.get_size_in_gb(aggregated_vector)
        reduction_count += time.time() - red_time
        ones_count = (aggregated_vector - sum_r_one_freq)[:-num_new_elements]
        red_time = time.time()
        total_comm_size += utilities.get_size_in_gb(twos)
        reduction_count += time.time() - red_time
        communication.send_data_to_server(party_socket, twos, p)
        del twos
        aggregated_vector = communication.receive_data_from_server(party_socket)
        red_time = time.time()
        total_comm_size += utilities.get_size_in_gb(aggregated_vector)
        reduction_count += time.time() - red_time
        twos_count = (aggregated_vector - sum_r_two_freq)[:-num_new_elements]
        red_time = time.time()
        total_comm_size += utilities.get_size_in_gb(sums)
        reduction_count += time.time() - red_time
        communication.send_data_to_server(party_socket, sums, p)
        del sums
        aggregated_vector = communication.receive_data_from_server(party_socket)
        red_time = time.time()
        total_comm_size += utilities.get_size_in_gb(aggregated_vector)
        reduction_count += time.time() - red_time
        means = ((aggregated_vector - sum_r_total_freq) / N)[:-num_new_elements]
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
            sum_r_sq = np.array(sq_diffs.power(2).sum(axis=0)).ravel()
            stds[start_col:end_col] += sum_r_sq
            start_col = end_col
        stds = (stds / (N - 1))
        del start_col, end_col, local_means, sq_diffs, sum_r_sq
        stds_random = stds[rand_indices] + np.random.choice([-0.1, 0.1], size=num_new_elements)
        stds = np.append(stds, stds_random)
        r_stds, sum_r_stds = GWAS_lib.randomized_encoding_addtion(len(stds), P, p, 5)
        stds = stds + r_stds
        red_time = time.time()
        total_comm_size += utilities.get_size_in_gb(stds)
        reduction_count += time.time() - red_time
        communication.send_data_to_server(party_socket, stds, p)
        del stds
        aggregated_vector = communication.receive_data_from_server(party_socket)
        red_time = time.time()
        total_comm_size += utilities.get_size_in_gb(aggregated_vector)
        reduction_count += time.time() - red_time

        stds = np.sqrt(np.maximum(aggregated_vector - sum_r_stds, 0))[:-num_new_elements]

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
            block_number = B_blocked[n][i]
            block_index = block_number - 1
            block_size = X_list[i].shape[1]
            additional_cols = extra_elements[block_index] + initial_block_size - block_size
            GWAS_lib.generate_mask_X_block(M, B, p, P, N, C, block_size, additional_cols, block_index)
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
            block_number = B_blocked[n][_]
            while not file_loaded and attempts < max_attempts:
                try:
                    O_X[_] = load_npz(
                        '../test_site/Data/N{}_M{}_C{}_P{}_B{}/Party_{}/Masks/O_X_block_{}.npz'.format(
                            N, M, C, P, B, p, block_number
                        )
                    )
                    file_loaded = True
                except Exception as e:
                    attempts += 1
                    time.sleep(delay)
            if not file_loaded:
                print("Failed to load the file after multiple attempts. Increase time for O_X.")

        current_count = len(X_list)
        right = [[] for _ in range(current_count)]
        masked_X = [[] for _ in range(current_count)]
        for i, X in enumerate(X_list):
            X = X.astype(np.float32)
            O_X[i] = O_X[i].astype(np.float32)

            right[i] = dot_product_mkl(X, O_X[i].T)
            O_Z_p=O_Z_p.astype(np.float32)

            masked_X[i] = dot_product_mkl(O_Z_p, right[i])
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
            right = dot_product_mkl(O_X[_], S[_])
            X = X.astype(np.float32)
            right = right.astype(np.float32)
            intermediate = dot_product_mkl(X, right)
            del right
            O_Z = O_Z.astype(np.float32)
            X_tilde[_] = dot_product_mkl(O_Z.T, intermediate)[start_index:end_index, :]
            del intermediate

            #####################################################################################################
            # --------------------------------------LEVEL 1 RIDGE REGRESSION-------------------------------------#
            #####################################################################################################

            block_number = B_blocked[n][_]  # 1-based block id for this X_tilde[_]
            block_index = block_number - 1  # 0-based index used for filenames / seeds

            GWAS_lib.generate_mask_beta_block(
                M, B, p, P, N, C,
                X_tilde[_].shape[1],
                block_index
            )
        O_b=[[] for _ in range(current_count)]
        sending=[[] for _ in range(current_count)]
        for _ in range(current_count):
            file_loaded = False
            max_attempts = 1000
            attempts = 0
            delay = 0.025
            while not file_loaded and attempts < max_attempts:
                try:
                    O_b[_] = load_npz(
                        '../test_site/Data/N{}_M{}_C{}_P{}_B{}/Party_{}/Masks/O_b_block_{}.npz'.format(N, M, C, P, B, p,
                                                                                                       (_ + (
                                                                                                                   n * bulk)) + 1))
                    file_loaded = True
                except Exception as e:
                    attempts += 1
                    time.sleep(delay)
        for _ in range(current_count):
            X_tilde[_]=X_tilde[_].astype(np.float32)
            O_b[_]=O_b[_].astype(np.float32)
            intermediate = dot_product_mkl(X_tilde[_], O_b[_])
            sending[_] = dot_product_mkl(O_y_tilde_p, intermediate)
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
                intermediate = dot_product_mkl(X_tilde_k, O_b[_])
                X_tilde_k = dot_product_mkl(O_y_tilde_p, intermediate)
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
            intermediate = dot_product_mkl(O_y_tilde_p, X_tilde[_])
            D = D.astype(np.float32)
            result_matrix= dot_product_mkl(intermediate,D)
            red_time = time.time()
            total_comm_size += utilities.get_size_in_gb(result_matrix)
            reduction_count += time.time() - red_time
            communication.send_data_to_server(party_socket, result_matrix, p)
        del result_matrix

    print(f'TIME TAKEN FOR {B} blocks is {time.time()-party_start_time-reduction_count}')
    #print(f'reduce from server {reduction_count}')
    print(f'SIZE OF DATA IN TOTAL FOR {B} blocks is {total_comm_size}')
    directory = '../test_site/Data/N{}_M{}_C{}_P{}_B{}/Party_{}/Masks'.format(N, M, C, P, B, p)
    if os.path.exists(directory):
        shutil.rmtree(directory)
    party_socket.close()

if __name__ == "__main__":
    main()

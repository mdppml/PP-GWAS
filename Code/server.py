import mkl
import argparse
import socket
import numpy as np
import gc
import time
from scipy.sparse.linalg import cg
import os
import GWAS_lib
import shutil
import utilities
import communication
from scipy.stats import chi2
import warnings

os.environ["MKL_NUM_THREADS"] = "16"
os.environ["MKL_DYNAMIC"] = "TRUE"

with warnings.catch_warnings():
    warnings.simplefilter("ignore")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--number_of_clients", type=int, required=True)
    parser.add_argument("--base_port", type=int, required=True)
    parser.add_argument("--number_of_samples", type=int, required=True)
    parser.add_argument("--number_of_snps", type=int, required=True)
    parser.add_argument("--number_of_covariates", type=int, required=True)
    parser.add_argument("--number_of_blocks", type=int, required=True)
    parser.add_argument("--number_of_folds", type=int, required=True)
    parser.add_argument("--number_of_blocks_per_run", type=int, required=True)
    args = parser.parse_args()
    B, K, N, M, P, C, R = args.number_of_blocks, args.number_of_folds, args.number_of_samples, args.number_of_snps, args.number_of_clients, args.number_of_covariates, 5
    bulk=args.number_of_blocks_per_run

    server_starting = time.time()
    pre_gwas_time=0
    level_0_time=0
    level_1_time=0
    server_sockets = []
    client_sockets = []
    total_save=0
    total_removed=0
    for _ in range(P):
        print(f'loop number {_}')
        server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        hostname = socket.gethostbyname(socket.gethostname())
        if _==0:
            with open('/home/swaminathan/ppREGENIE/Data/ip_address_file.txt', 'w') as file:
                file.write(hostname)
        server.setsockopt(socket.SOL_SOCKET, socket.SO_KEEPALIVE, 1)
        server.bind((hostname, args.base_port + _ + 1))
        with open(f'/home/swaminathan/ppREGENIE/Data/server_ready_{_+1}.txt', 'w') as f:
            f.write('Server is ready')
        server.listen(1)
        print(f"Server is listening on port {args.base_port + _ + 1}")
        client_socket, addr = server.accept()
        print(f"Accepted connection from {addr}")
        client_sockets.append(client_socket)
        server_sockets.append(server)
    print(f'SERVER: time taken to establish connections is {time.time() - server_starting}')
    utilities.print_memory_usage()

    directory = '/home/swaminathan/ppREGENIE/Data/N{}_M{}_C{}_P{}_B{}/Masks'.format(N, M, C, P, B)
    if os.path.exists(directory):
        shutil.rmtree(directory)
    os.makedirs(directory, exist_ok=True)
    directory = '/home/swaminathan/ppREGENIE/Data/N{}_M{}_C{}_P{}_B{}/Server'.format(N, M, C, P, B)
    if os.path.exists(directory):
        shutil.rmtree(directory)
    os.makedirs(directory, exist_ok=True)
    
    collection_time = time.time()
    first_two = time.time()
    aggregated_vector = communication.receive_and_sum(client_sockets)
    for client_socket in client_sockets:
        communication.send_data_to_client(client_socket, aggregated_vector)
    aggregated_vector = communication.receive_and_sum(client_sockets)
    for client_socket in client_sockets:
        communication.send_data_to_client(client_socket, aggregated_vector)

    #####################################################################################################
    #-----------------------DISTRIBUTED PROJECTION OF COVARIATES AND STANDARDIZING----------------------#
    #####################################################################################################

    z_receive_time = time.time()
    total_Z_matrix, masked_Z = communication.receive_sum_and_append(client_sockets)
    ready_Z = []
    z_cal_time = time.time()
    for p in range(P):
        ZTZinv = np.linalg.inv(total_Z_matrix.T @ total_Z_matrix)
        ready_Z.append(masked_Z[p] @ ZTZinv @ total_Z_matrix.T)
        del ZTZinv
    del masked_Z, total_Z_matrix
    gc.collect()

    y_receive_time = time.time()
    total_Y_matrix, masked_Y = communication.receive_sum_and_append(client_sockets)
    Y_tilde = []
    y_cal_time = time.time()
    for p in range(P):
        Y_tilde.append(masked_Y[p] - ready_Z[p] @ total_Y_matrix)
    del masked_Y, total_Y_matrix
    y_send_time = time.time()
    communication.send_data_to_clients(client_sockets, Y_tilde)
    del Y_tilde

    gamma_Y_tilde = []
    gamma_Y_tilde2 = np.zeros((N + 1, 1))
    k_fold_time = time.time()
    for k in range(K):
        Y_tilde_masked = communication.receive_and_sum(client_sockets)
        gamma_Y_tilde.append(Y_tilde_masked)
        Y_tilde_masked = communication.receive_and_sum(client_sockets)
        gamma_Y_tilde2 += Y_tilde_masked
    saving=0

    B_blocked = GWAS_lib.split_number(B, bulk)
    num_loops = len(B_blocked)

    X_cal=[]
    W_1 = [np.zeros((N + 1, B * R)) for _ in range(K)]
    W_2 = [np.zeros((N + 1, B * R)) for _ in range(K)]

    l = np.linspace(0.01, 0.99, R)
    l = M * (1 - l ** 2) / l ** 2

    sum_all = [np.zeros((N + 1, R)) for _ in range(K)]
    sum2_all = [np.zeros((N + 1, R)) for _ in range(K)]

    pre_gwas_time+=(time.time()-server_starting)
    for n in range(num_loops):
        if n>0:
            with open(f'/home/swaminathan/ppREGENIE/Data/server_ready_loop_{n}.txt', 'w') as f:
                f.write('Server is ready')
        timing = time.time()
            
        #####################################################################################################
        #------------------------------------------QUALITY CONTROL------------------------------------------#
        #####################################################################################################

        aggregated_vector = communication.receive_and_sum(client_sockets)
        for client_socket in client_sockets:
            communication.send_data_to_client(client_socket, aggregated_vector)
        aggregated_vector = communication.receive_and_sum(client_sockets)
        for client_socket in client_sockets:
            communication.send_data_to_client(client_socket, aggregated_vector)
        aggregated_vector = communication.receive_and_sum(client_sockets)
        for client_socket in client_sockets:
            communication.send_data_to_client(client_socket, aggregated_vector)
        aggregated_vector = communication.receive_and_sum(client_sockets)
        for client_socket in client_sockets:
            communication.send_data_to_client(client_socket, aggregated_vector)

        #####################################################################################################
        #-----------------------DISTRIBUTED PROJECTION OF COVARIATES AND STANDARDIZING----------------------#
        #####################################################################################################
        
        current_count=len(B_blocked[n])

        total_X_matrix=[[] for _ in range(current_count)]
        masked_X=[[] for _ in range(current_count)]
        for _ in range(current_count):
            total_X_matrix[_], masked_X[_]=communication.receive_sum_and_append(client_sockets)
            communication.send_data_to_clients(client_sockets,aggregated_vector)
        for _ in range(current_count):
            utilities.print_sys_info()
            X_tilde = []

            for idx in range(P):
                result = masked_X[_][idx] - (ready_Z[idx] @ total_X_matrix[_])
                X_tilde.append(result)

            communication.send_data_to_clients(client_sockets, X_tilde)
        del total_X_matrix, masked_X

        #####################################################################################################
        #-----------------------------LEVEL 1 RIDGE REGRESSION: DISTRIBUTED-ADMM----------------------------#
        #####################################################################################################

        gamma_X_beta=[[] for _ in range(current_count)]
        gamma_X_beta2=[[] for _ in range(current_count)]

        total = [[] for _ in range(current_count)]
        for _ in  range(current_count):
            total[_]=communication.receive_and_sum(client_sockets)
        X_tilde_masked_k=[[] for _ in range(current_count)]
        for _ in range(current_count):
            for k in range(K):
                X_tilde_masked_k[_] = communication.receive_and_sum(client_sockets)
                gamma_X_beta[_].append(X_tilde_masked_k[_])
                gamma_X_beta2[_].append(total[_] - X_tilde_masked_k[_])
        del total

        for _ in range(current_count):
            gamma=communication.receive_and_sum(client_sockets)
            X_cal.append(gamma)
        del gamma

        utilities.print_memory_usage()
        pre_gwas_time+=time.time()-timing

        #####################################################################################################
        #-----------------------------LEVEL 2 RIDGE REGRESSION: DISTRIBUTED-CGD-----------------------------#
        #####################################################################################################
        
        calc = time.time()
        for b in range(current_count):
            precomputed_right = [gb.T @ gy for gb, gy in zip(gamma_X_beta[b], gamma_Y_tilde)]
            precomputed_inside = [gb.T @ gb for gb in gamma_X_beta[b]]

            I = np.eye(precomputed_inside[0].shape[0])
            for k in range(K):
                right = precomputed_right[k]
                inside = precomputed_inside[k]
                W2_all = []
                for r in l:
                    left = inside + r * I
                    beta, info = cg(left, right)
                    ins = gamma_X_beta2[b][k] @ beta.reshape(-1, 1)
                    W2_all.append(ins)
                    del left, beta, ins, info

                W2_all = np.concatenate(W2_all, axis=1)
                for k2, (s, s2) in enumerate(zip(sum_all, sum2_all)):
                    if k2 != k:
                        s += W2_all
                    else:
                        s2 += W2_all
                del W2_all
            del precomputed_inside
        level_0_time+=time.time()-calc
    for server in server_sockets:
        server.close()

    #####################################################################################################
    #-----------------------------DISTRIBUTED SINGLE SNP ASSOCIATION TESTING----------------------------#
    #####################################################################################################

    utilities.print_memory_usage()
    level1=time.time()
    y_hat = []
    l = np.linspace(0.01, 0.99, R)
    for _ in range(len(l)):
        l[_] = (B * R) * (1 - (l[_] ** 2)) / (l[_] ** 2)
    for r in range(R):
        y_hat.append(np.zeros((N + 1, 1)))
    for k in range(K):
        inside = (W_1[k].T @ W_1[k])
        rightside = W_1[k].T @ gamma_Y_tilde[k]
        for r in range(R):
            y_hat[r] += W_2[k] @ (np.linalg.inv(inside + (l[r] * np.eye(B * R))) @ (rightside))
    del W_1, W_2
    score = 1000000
    index = 0
    y_hatty = np.empty((N + 1, 1))
    for r in range(R):
        temp = np.linalg.norm(gamma_Y_tilde2 - y_hat[r])
        if temp < score:
            score = temp
            index = r
            y_hatty = y_hat[r]
    y = gamma_Y_tilde2
    sigma = (1 / (N - C)) * (np.linalg.norm(y - y_hatty)) ** 2
    sigma_inv = 1 / sigma
    t = []
    total_columns= __builtins__.sum(matrix.shape[1] for matrix in X_cal)
    gc.collect()
    for b in range(B):
        for i in range(X_cal[b].shape[1]):
            xi = X_cal[b][:, i]
            t.append(sigma_inv * (xi.T @ (y - y_hatty)) ** 2 / (xi.T @ xi))
    t=np.array(t)
    df = 1
    p_values = 1 - chi2.cdf(t, df)
    positions = np.arange(len(t))
    neg_log_p = -np.log10(p_values)
    np.save('/home/swaminathan/ppREGENIE/Data/N{}_M{}_C{}_P{}_B{}/neg_log_transfer.npy'.format(N, M, C, P, B),neg_log_p)
    level_1_time=time.time()-level1
    print(f'TIME RESULTS')
    print(f'Time taken for pre_gwas_communication and computations - {pre_gwas_time}')
    print(f'Time taken for level 0 - {level_0_time}')
    print(f'Time taken for level 1 - {level_1_time}')
    print(f'TOTAL TIME FOR PROGRAM is {pre_gwas_time+level_0_time+level_1_time}')
    print(f'Formatted - {utilities.convert_seconds_to_dhms(pre_gwas_time+level_0_time+level_1_time)}')
    directory = '/home/swaminathan/ppREGENIE/Data/N{}_M{}_C{}_P{}_B{}/Masks'.format(N, M, C, P, B)
    if os.path.exists(directory):
        shutil.rmtree(directory)
    os.makedirs(directory, exist_ok=True)


if __name__ == "__main__":
    main()

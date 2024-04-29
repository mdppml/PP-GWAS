import mkl
import argparse
import socket
import numpy as np
import gc
import time
from scipy.sparse.linalg import cg
import os
import GWAS_lib
from scipy.linalg import lu_factor, lu_solve, solve_triangular
import psutil
from sparse_dot_mkl import dot_product_mkl
import shutil
import cProfile
from scipy.stats import chi2
import dask.array as da
import dask
import matplotlib.pyplot as plt
from scipy.sparse import load_npz, diags, issparse, save_npz, csr_matrix, csc_matrix, isspmatrix, hstack
import warnings
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import ProcessPoolExecutor
os.environ["MKL_NUM_THREADS"] = "16"
os.environ["MKL_DYNAMIC"] = "TRUE"
from scipy.sparse.linalg import cg, gmres, bicgstab

pid = os.getpid()
py = psutil.Process(pid)

def load_file(filename):
    return np.load(filename, allow_pickle=True)

def print_matrix_info(matrix):
    is_sparse = isspmatrix(matrix)
    dtype = matrix.dtype
    size = matrix.size
    shape = matrix.shape
    if is_sparse:
        non_zero_count = matrix.nnz  
    else:
        non_zero_count = np.count_nonzero(matrix)

    sparsity = 1 - (non_zero_count / size)
    print(f"Matrix = {type(matrix)}")
    print(f"Matrix Type: {dtype}")
    print(f"Matrix Size: {size}")
    print(f"Matrix Shape: {shape}")
    print(f"Matrix Sparsity: {sparsity * 100:.2f}%")
    if dtype == np.float32:
        print("The matrix is of type float32.")
    elif dtype == np.float64:
        print("The matrix is of type float64.")
    else:
        print("The matrix is neither float32 nor float64.")


def print_sys_info():
    cpu_percent = psutil.cpu_percent(interval=1)
    memory_info = psutil.virtual_memory()

    p = psutil.Process()
    num_threads = p.num_threads()

    num_cpu_cores = psutil.cpu_count()

    mkl_threads = os.environ.get('MKL_NUM_THREADS', 'Not Set')

    num_threads_not_mkl = num_threads if mkl_threads == 'Not Set' else num_threads - int(mkl_threads)

    print(f"--- System Info ---")
    print(f"CPU usage: {cpu_percent}%")
    print(f"Total memory: {memory_info.total / (1024 * 1024)} MB")
    print(f"Available memory: {memory_info.available / (1024 * 1024)} MB")
    print(f"Used memory: {memory_info.used / (1024 * 1024)} MB")
    print(f"Memory percent used: {memory_info.percent}%")
    print(f"Total threads used by process: {num_threads}")
    print(f"Total available threads (cores): {num_cpu_cores}")
    print(f"MKL threads: {mkl_threads}")
    print(f"Threads not used by MKL: {num_threads_not_mkl}")
    max_mkl_threads = mkl.get_max_threads()
    print(f"Max threads MKL is allowed to use: {max_mkl_threads}")
    print("--------------------")


with warnings.catch_warnings():
    warnings.simplefilter("ignore")

def compute_sparsity(matrix):
    if isinstance(matrix, np.ndarray):
        total_elements = matrix.size
        zero_elements = np.count_nonzero(matrix == 0)
    elif hasattr(matrix, "nnz") and hasattr(matrix, "shape"):
        total_elements = matrix.shape[0] * matrix.shape[1]
        zero_elements = total_elements - matrix.nnz

    else:
        raise ValueError("Unknown matrix type. Expected numpy ndarray or scipy sparse matrix.")

    sparsity = zero_elements / total_elements

    return sparsity

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
    print(f'Starting program for P={P}, N={N}, M={M}, C={C}')


    server_starting = time.time()
    pre_gwas_time=0
    level_0_time=0
    level_1_time=0
    server_sockets = []
    client_sockets = []
    total_save=0
    total_removed=0
    for _ in range(P):
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
    GWAS_lib.print_memory_usage()

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
    aggregated_vector = GWAS_lib.receive_and_sum(client_sockets)
    for client_socket in client_sockets:
        GWAS_lib.send_data_to_client(client_socket, aggregated_vector)
    aggregated_vector = GWAS_lib.receive_and_sum(client_sockets)
    for client_socket in client_sockets:
        GWAS_lib.send_data_to_client(client_socket, aggregated_vector)

    z_receive_time = time.time()
    total_Z_matrix, masked_Z = GWAS_lib.receive_sum_and_append(client_sockets)
    ready_Z = []
    z_cal_time = time.time()
    for p in range(P):
        ZTZinv = np.linalg.inv(total_Z_matrix.T @ total_Z_matrix)
        ready_Z.append(masked_Z[p] @ ZTZinv @ total_Z_matrix.T)
        del ZTZinv
    del masked_Z, total_Z_matrix
    gc.collect()

    y_receive_time = time.time()
    total_Y_matrix, masked_Y = GWAS_lib.receive_sum_and_append(client_sockets)
    Y_tilde = []
    y_cal_time = time.time()
    for p in range(P):
        Y_tilde.append(masked_Y[p] - ready_Z[p] @ total_Y_matrix)
    del masked_Y, total_Y_matrix
    y_send_time = time.time()
    GWAS_lib.send_data_to_clients(client_sockets, Y_tilde)
    del Y_tilde

    gamma_Y_tilde = []
    gamma_Y_tilde2 = np.zeros((N + 1, 1))
    k_fold_time = time.time()
    for k in range(K):
        Y_tilde_masked = GWAS_lib.receive_and_sum(client_sockets)
        gamma_Y_tilde.append(Y_tilde_masked)
        Y_tilde_masked = GWAS_lib.receive_and_sum(client_sockets)
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

        timing = time.time()

        aggregated_vector = GWAS_lib.receive_and_sum(client_sockets)
        for client_socket in client_sockets:
            GWAS_lib.send_data_to_client(client_socket, aggregated_vector)
        aggregated_vector = GWAS_lib.receive_and_sum(client_sockets)
        for client_socket in client_sockets:
            GWAS_lib.send_data_to_client(client_socket, aggregated_vector)
        aggregated_vector = GWAS_lib.receive_and_sum(client_sockets)
        for client_socket in client_sockets:
            GWAS_lib.send_data_to_client(client_socket, aggregated_vector)
        aggregated_vector = GWAS_lib.receive_and_sum(client_sockets)
        for client_socket in client_sockets:
            GWAS_lib.send_data_to_client(client_socket, aggregated_vector)

        current_count=len(B_blocked[n])


        total_X_matrix=[[] for _ in range(current_count)]
        masked_X=[[] for _ in range(current_count)]
        for _ in range(current_count):
            total_X_matrix[_], masked_X[_]=GWAS_lib.receive_sum_and_append(client_sockets)
        for _ in range(current_count):
            X_tilde = list(map(lambda p: masked_X[_][p] - (ready_Z[p] @ total_X_matrix[_]), range(P)))
            GWAS_lib.send_data_to_clients(client_sockets, X_tilde)
        del total_X_matrix, masked_X

        gamma_X_beta=[[] for _ in range(current_count)]
        gamma_X_beta2=[[] for _ in range(current_count)]

        total = [[] for _ in range(current_count)]
        for _ in  range(current_count):
            total[_]=GWAS_lib.receive_and_sum(client_sockets)

        X_tilde_masked_k=[[] for _ in range(current_count)]
        for _ in range(current_count):
            for k in range(K):
                X_tilde_masked_k[_] = GWAS_lib.receive_and_sum(client_sockets)
                gamma_X_beta[_].append(X_tilde_masked_k[_])
                gamma_X_beta2[_].append(total[_] - X_tilde_masked_k[_])
        del total


        for _ in range(current_count):
            gamma=GWAS_lib.receive_and_sum(client_sockets)
            X_cal.append(gamma)
        del gamma

        GWAS_lib.print_memory_usage()
        print(f'DONE! Took for all {bulk} blocks, {time.time()-timing}')
        pre_gwas_time+=time.time()-timing


        calc = time.time()
        for b in range(current_count):
            sum = []
            sum2 = []
            for k in range(K):
                sum.append(np.zeros((N + 1, R)))
                sum2.append(np.zeros((N + 1, R)))
            for k in range(K):
                right = gamma_X_beta[b][k].transpose() @ gamma_Y_tilde[k]
                inside = (gamma_X_beta[b][k].transpose() @ gamma_X_beta[b][k])
                W2 = np.empty((N + 1, 0))
                I = np.eye(inside.shape[0])
                for r in range(R):
                    left = (inside + (l[r] * I))
                    beta, info = cg(left, right)
                    beta = np.reshape(beta, (beta.shape[0], 1))
                    del left
                    ins = gamma_X_beta2[b][k] @ beta
                    del beta
                    W2 = np.concatenate((W2, ins), axis=1)
                    del ins
                del right, inside
                for k2 in range(K):
                    if k2 != k:
                        sum[k2] += W2
                    else:
                        sum2[k2] += W2
                del W2
            for k in range(K):
                W_1[k] = np.concatenate((W_1[k], sum[k]), axis=1)
                W_2[k] = np.concatenate((W_2[k], sum2[k]), axis=1)
        level_0_time+=time.time()-calc
        print(f'END OF 60 blocks -')
    for server in server_sockets:
        server.close()

    GWAS_lib.print_memory_usage()
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
    np.save('/home/swaminathan/ppREGENIE/Data/N{}_M{}_C{}_P{}_B{}/neg_log_transfer.npy'.format(N, M, C, P, B),
            neg_log_p)
    level_1_time=time.time()-level1
    print(f'TIME RESULTS')
    print(f'Time taken for pre_gwas_communication and computations - {pre_gwas_time}')
    print(f'Time taken for level 0 - {level_0_time}')
    print(f'Time taken for level 1 - {level_1_time}')
    print(f'TOTAL TIME FOR PROGRAM is {pre_gwas_time+level_0_time+level_1_time}')
    print(f'Formatted - {GWAS_lib.convert_seconds_to_dhms(pre_gwas_time+level_0_time+level_1_time)}')
    directory = '/home/swaminathan/ppREGENIE/Data/N{}_M{}_C{}_P{}_B{}/Masks'.format(N, M, C, P, B)
    if os.path.exists(directory):
        shutil.rmtree(directory)
    os.makedirs(directory, exist_ok=True)


if __name__ == "__main__":
    cProfile.run('main()', '/home/swaminathan/ppREGENIE/Data/server2.pstats')

import numpy as np
from scipy.sparse import csr_matrix, vstack, hstack, diags, save_npz, csc_matrix, issparse, bmat
import gc
import random
import os
import struct
import psutil
from concurrent.futures import ThreadPoolExecutor
import concurrent.futures
import time
import pickle
from multiprocessing import Pool
from threading import Thread, Lock
import mkl
CHUNK_SIZE = 60000
from scipy.stats import chisquare

def serialize_matrix_client_pickle(client_id,matrix):
    # Function to serialize matrix on the client side using pickle
    data_to_send=[client_id,matrix]
    data_to_send=pickle.dumps(data_to_send)
    return data_to_send

def serialize_matrix_client(client_id, matrix):
    # Custom function to serialize data on the client side using various strategies based on data type
    if isinstance(matrix, np.ndarray):
        if len(matrix.shape) == 1:
            matrix_type = 0  # 1D array
            rows = len(matrix)
            cols = 0
            matrix_data = matrix.astype(np.float32).tobytes()
        else:
            rows, cols = matrix.shape
            matrix_data = matrix.astype(np.float32).tobytes()
            if rows > 1 and cols > 1:
                matrix_type = 3  # General 2D array
            elif rows == 1:
                matrix_type = 2  # Row vector
            else:
                matrix_type = 1  # Column vector
    elif isinstance(matrix, (csr_matrix, csc_matrix)):
        rows, cols = matrix.shape
        matrix_type = 4 if isinstance(matrix, csr_matrix) else 5  # CSR or CSC
        matrix_data = pickle.dumps(matrix)
    elif isinstance(matrix, list):
        if all(isinstance(x, np.ndarray) and len(x.shape) == 1 for x in matrix):
            matrix_type = 6  # list of 1D arrays
            serialized_matrices = []
            for mat in matrix:
                serialized_matrices.append(serialize_matrix_client(client_id, mat))
            matrix_data = b''.join(serialized_matrices)
            rows = len(matrix)
            cols = 0  # not applicable for list of 1D arrays
        else:
            raise ValueError("Unsupported list elements.")
    else:
        raise ValueError("Unsupported matrix type.")

    header = struct.pack('!4I1Q', client_id, rows, cols, matrix_type, len(matrix_data))
    return header + matrix_data


def serialize_matrix_server(matrix):
    # Custom function to serialize data on the server side using pickle
    return pickle.dumps(matrix)


def send_data_to_server(client_socket, matrix, client_id):
    # Function to send both the matrix data and client_id to server
    tt=time.time()
    data = serialize_matrix_client(client_id, matrix)
    data_size = struct.pack('!Q', len(data))
    tt=time.time()
    client_socket.sendall(data_size)
    client_socket.send(data)

def send_data_to_client(client_socket, matrix):
    # Function to send the aggregated results to one client via their respective socket
    data = serialize_matrix_server(matrix)
    data_size = struct.pack('!Q', len(data))
    client_socket.sendall(data_size)
    chunk_count = 0
    for i in range(0, len(data), CHUNK_SIZE):
        chunk = data[i:i + CHUNK_SIZE]
        client_socket.sendall(chunk)
        chunk_count += len(chunk)

def send_data_to_clients(client_sockets,matrices):
    # Function to send the aggregated results to all clients via their respective sockets
    mkl.set_num_threads(4)
    with ThreadPoolExecutor(max_workers=len(client_sockets)) as executor:
        futures = []
        for p in range(len(client_sockets)):
            futures.append(executor.submit(send_data_to_client, client_sockets[p], matrices[p]))

        for future in futures:
            future.result()
    mkl.set_num_threads(16)

def deserialize_matrix_server_pickle(data):
    # Function to deserialize data on the server side using pickle
    data=pickle.loads(data)
    client_id=data[0]
    matrix=data[1]
    return client_id, matrix

def deserialize_matrix_server(serialized_data):
    # Function to deserialize data on the server side using custom deserialization based on data type
    header_length = 4 * 4 + 1 * 8  
    client_id, rows, cols, matrix_type, data_length = struct.unpack('!4I1Q', serialized_data[:header_length])
    matrix_data = serialized_data[header_length:]
    del serialized_data
    if matrix_type == 0:  # 1D array
        matrix = np.frombuffer(matrix_data, dtype=np.float32)
    elif matrix_type in [1, 2, 3]:  # 2D array, Row vector, Column vector
        matrix = np.frombuffer(matrix_data, dtype=np.float32).reshape(rows, cols if cols else 1)
    elif matrix_type == 4:  # CSR
        matrix = pickle.loads(matrix_data)

    elif matrix_type == 5:  # CSC
        matrix = pickle.loads(matrix_data)

    elif matrix_type == 6:  # list of 1D arrays
        matrix = []
        offset = 0
        for _ in range(rows):  
            sub_client_id, sub_rows, _, sub_matrix_type, sub_data_length = struct.unpack(
                '!4I1Q', matrix_data[offset:offset + header_length]
            )
            offset += header_length

            sub_matrix_data = matrix_data[offset:offset + sub_data_length]
            offset += sub_data_length
            matrix.append(np.frombuffer(sub_matrix_data, dtype=np.float32))
    else:
        raise ValueError("Unsupported matrix type.")

    return client_id, matrix


def deserialize_matrix_client(data):
    # Function to deserialize data on the client side using pickle
    return pickle.loads(data)


def receive_and_sum(client_sockets):
    # Function to receive data from all clients and aggregate them
    raw_data = [None] * len(client_sockets)
    threads = []
    mkl.set_num_threads(1)

    def threaded_receive_data(client_socket, results, index):
        results[index] = receive_data(client_socket)

    for i, client_socket in enumerate(client_sockets):
        t = Thread(target=threaded_receive_data, args=(client_socket, raw_data, i))
        t.start()
        threads.append(t)

    for t in threads:
        t.join()

    mkl.set_num_threads(16)

    total_X = None  # for matrices
    total_X_arrays = None  # for lists of 1D arrays
    masked_X = []  # to hold matrices for summing

    for data in raw_data:
        client_id, received_item = deserialize_matrix_server(data)

        if isinstance(received_item, list):  
            if total_X_arrays is None:
                total_X_arrays = [np.zeros_like(arr) for arr in received_item]

            for i, arr in enumerate(received_item):
                total_X_arrays[i] += arr

        else:  
            masked_X.append(received_item)

    if total_X_arrays is not None and masked_X:
        raise ValueError("Inconsistent data types received: both matrices and lists of arrays.")

    if masked_X:
        total_X = np.sum(masked_X, axis=0)

    return total_X_arrays if total_X_arrays is not None else total_X


def receive_data(client_socket):
    # Function to receive data from one client
    tt=time.time()
    data_size_bytes = client_socket.recv(8)
    if len(data_size_bytes) != 8:
        raise ConnectionError("Client disconnected prematurely while sending data size.")
    tt=time.time()
    data_size = struct.unpack('!Q', data_size_bytes)[0]
    tt=time.time()
    chunks = []
    received_size = 0
    while received_size < data_size:
        chunk = client_socket.recv(min(8192, data_size - received_size))
        if not chunk:
            break
        chunks.append(chunk)
        received_size += len(chunk)
    if received_size < data_size:
        raise ConnectionError(f"Expected {data_size} bytes but received only {received_size} bytes.")

    return b''.join(chunks)

def receive_sum_and_append(client_sockets):
    # Function to receive data from all clients, aggregate them and also store them
    mkl.set_num_threads(1)
    raw_data = [None] * len(client_sockets)
    threads = []

    def threaded_receive_data(client_socket, results, index):
        results[index] = receive_data(client_socket)

    for i, client_socket in enumerate(client_sockets):
        t = Thread(target=threaded_receive_data, args=(client_socket, raw_data, i))
        t.start()
        threads.append(t)

    for t in threads:
        t.join()
    mkl.set_num_threads(16)

    masked_X = [None] * len(raw_data)

    for data in raw_data:
        client_id, matrix = deserialize_matrix_server(data)
        masked_X[client_id - 1] = matrix
        del data

    del raw_data

    total_X_matrix = np.sum(masked_X, axis=0)

    return total_X_matrix, masked_X


def receive_data_from_server(server_socket):
    # Function to receive data from the server
    data_size_bytes = server_socket.recv(8)  
    data_size = struct.unpack('!Q', data_size_bytes)[0]

    received_data = bytearray()
    while len(received_data) < data_size:
        chunk = server_socket.recv(min(CHUNK_SIZE, data_size - len(received_data)))
        if not chunk:
            raise ConnectionError("Connection closed before all data was received.")
        received_data.extend(chunk)

    matrix = deserialize_matrix_client(received_data)
    return matrix


#### DATA GENERATION TOOLS

np.random.seed(420)

def get_delta(C, seed):
    # QR Decomposition
    A = np.random.rand(C + 1, C)
    Q, R = np.linalg.qr(A)
    return Q

def get_gamma(N, K, seed):
    # Randomized Mask Generation 
    np.random.seed(seed)
    min_block_size = 95
    max_block_size = 100

    num_blocks = N // max_block_size
    while num_blocks <= K:
        min_block_size -= 5
        max_block_size -= 5
        num_blocks = N // max_block_size
    remaining_columns = N % max_block_size

    blocks = []
    for i in range(K):
        block_size = np.random.randint(min_block_size, max_block_size + 1)

        block = generate_semi_orthogonal_matrix(block_size + 1, block_size)
        blocks.append(block)

    for i in range(num_blocks - 1 - K):
        block_size = np.random.randint(min_block_size, max_block_size + 1)

        block = generate_random_orthogonal_matrix(block_size)

        blocks.append(block)

    last_block_size = N - sum(block.shape[1] for block in blocks)
    last_block = generate_random_orthogonal_matrix(last_block_size)
    blocks.append(last_block)
    gamma_blocks = []
    for i, block in enumerate(blocks):
        num_zero_rows = sum(block.shape[0] for block in blocks[:i])
        num_zero_cols = sum(block.shape[1] for block in blocks[:i])
        upper_padding = csr_matrix((num_zero_rows, N))
        lower_padding = csr_matrix((N + K - num_zero_rows - block.shape[0], N))
        left_padding = csr_matrix((block.shape[0], num_zero_cols))
        right_padding = csr_matrix((block.shape[0], N - num_zero_cols - block.shape[1]))
        padded_block = hstack((left_padding, block, right_padding)).tocsc()
        padded_block = vstack((upper_padding, padded_block, lower_padding)).tocsc()

        gamma_blocks.append(padded_block)

    gamma = random_permute_columns(sum(gamma_blocks), 0).tocsr()
    gamma = gamma.astype(np.float32)
    return gamma


def get_square_gamma(N, seed):
    # Randomized Square Mask Generation
    np.random.seed(seed)
    min_block_size = 95
    max_block_size = 100

    num_blocks = N // max_block_size
    remaining_columns = N % max_block_size

    blocks = []
    for i in range(num_blocks):
        if i != num_blocks - 1 or remaining_columns == 0:
            block_size = max_block_size
        else:
            block_size = remaining_columns

        block = generate_random_orthogonal_matrix(block_size)

        blocks.append(block)

    block_matrix = [[None] * len(blocks) for _ in range(len(blocks))]
    for i, block in enumerate(blocks):
        block_matrix[i][i] = block

    gamma = bmat(block_matrix)

    gamma = random_permute_square_columns(gamma, seed)
    gamma = gamma.tocsr()
    gamma = gamma.astype(np.float32)
    return gamma


def random_permute_square_columns(matrix, seed):
    # Function to randomly permute columns of a square matrix
    num_columns = matrix.shape[1]
    np.random.seed(seed)
    permuted_indices = np.random.permutation(num_columns)

    if not isinstance(matrix, np.ndarray):
        matrix = matrix.tocsr()

    permuted_matrix = matrix[:, permuted_indices]

    return permuted_matrix


def random_permute_columns(matrix, seed):
    # Function to randomly permute columns of a matrix
    num_columns = matrix.shape[1]
    np.random.seed(seed)
    permuted_indices = np.random.permutation(num_columns)

    permuted_matrix = matrix[:, permuted_indices]

    return permuted_matrix


def generate_random_orthogonal_matrix(size):
    # Function to generate a random orthogonal matrix
    random_matrix = np.random.randn(size, size)
    tol = 1e-7
    Q, R = np.linalg.qr(random_matrix)
    Q[abs(Q) < tol] = 0
    return Q


def generate_semi_orthogonal_matrix(rows, cols):
    # Funciton to generate a random semi-orthogonal matrix
    A = np.random.randn(rows, cols)
    for i in range(cols):
        for j in range(i):
            A[:, i] -= np.dot(A[:, i], A[:, j]) / np.dot(A[:, j], A[:, j]) * A[:, j]

        A[:, i] /= np.linalg.norm(A[:, i])

    for i in range(rows):
        for j in range(i + 1, cols):
            dot_product = abs(np.dot(A[i, :], A[j, :]))
            assert dot_product < 1, f"Dot product between row {i} and row {j} is {dot_product}, not less than 1"

    return A


def standardization_normalizers(M, P, p, seed):
    np.random.seed(seed)
    matrix = np.random.randint(-10, 10, size=(P, M), dtype=np.int8)

    column_sum_vector = np.sum(matrix, axis=0, dtype=np.int32)
    return matrix[p - 1, :], column_sum_vector


def generate_beta(M, B, p, P, N, C, block_size,b):
    # Function to generate betas
    if not os.path.exists(
            '/Data/N{}_M{}_C{}_P{}_B{}/Masks'.format(N, M, C, P, B)):
        os.makedirs('/Data/N{}_M{}_C{}_P{}_B{}/Masks'.format(N, M, C, P, B),
                    exist_ok=True)
    gamma = get_gamma(block_size, 0, b)
    filename = os.path.join(
        '/Data/N{}_M{}_C{}_P{}_B{}/Masks'.format(N, M, C, P, B),
        "beta_block_{}.npz".format(b + 1))
    save_npz(filename, gamma)
    del gamma
    gc.collect()

def save_matrices(X,b,N,M,C,P,B):
    if not os.path.exists(
            '/N{}_M{}_C{}_P{}_B{}/Server'.format(N, M, C, P, B)):
        os.makedirs('/Data/N{}_M{}_C{}_P{}_B{}/Server'.format(N, M, C, P, B),
                    exist_ok=True)
    np.save('/Data/N{}_M{}_C{}_P{}_B{}/Server/X_block_{}.npy'.format(N, M, C, P, B,b+1),X)

def generate_delta(M, B, p, P, N, C, block_size, additional_cols,b):
    directory = '/Data/N{}_M{}_C{}_P{}_B{}/Masks'.format(N, M, C, P, B)
    if not os.path.exists(directory):
        os.makedirs(directory, exist_ok=True)
    if additional_cols<0:
        print('additional_col = 0 found')
        additional_cols=0
    gamma = get_gamma(block_size, additional_cols, b).tocsc()
    filename = os.path.join(directory, "delta_block_{}.npz".format(b + 1))
    save_npz(filename, gamma)
    del gamma

    gc.collect()


def generate_Z_mask(M, B, p, P, N, C):
    directory = '/Data/N{}_M{}_C{}_P{}_B{}/Masks'.format(N, M, C, P, B)
    if not os.path.exists(directory):
        os.makedirs(directory, exist_ok=True)
    if p == 1:
        gamma = get_gamma(N, 1, 0)
        save_npz(
            '/Data/N{}_M{}_C{}_P{}_B{}/Masks/Z_mask.npz'.format(N, M, C, P, B),
            gamma)
        del gamma
        gc.collect()


def generate_gamma(M, B, K, p, P, N, C):
    # Function to generate Gamma
    if p == 1:
        gamma = get_gamma(N, 1, 1)
        save_npz(
            '/Data/N{}_M{}_C{}_P{}_B{}/Masks/gamma.npz'.format(N, M, C, P, B),
            gamma)
        del gamma
        gc.collect()


def generate_a_number(seed):
    random.seed(seed)
    value = random.random() + 1
    return value


#### DATA MANIPULATION TOOLS

def convert_to_numpy(matrix):
    if issparse(matrix):
        return matrix.toarray()
    else:
        return matrix


def calculate_param(k, N, K, p, P):
    k1 = k * (N / K)
    k2 = (k + 1) * (N / K)
    p1 = (p - 1) * (N / P)
    p2 = p * (N / P)
    end = 0
    if k1 <= p1:
        start = 0
        if k2 < p1:
            end = 0
        elif k2 < p2:
            end = k2 % (N / P)
        elif k2 > p2:
            end = N / P
    elif k1 > p2:
        start = 0
        end = 0
    else:
        start = k1 % (N / P)
        if k2 < p2:
            end = k2 % (N / P)
        else:
            end = N / P
    return int(start), int(end)


def hwe_chi_square(obs_counts):
    # Function to compute chi-square values
    total = sum(obs_counts)
    p = (2 * obs_counts[0] + obs_counts[1]) / (2 * total)
    q = 1 - p

    expected_counts = np.array([total * p ** 2, 2 * total * p * q, total * q ** 2])
    chi2, _ = chisquare(f_obs=obs_counts, f_exp=expected_counts)


    chi2, _ = chisquare(f_obs=obs_counts, f_exp=expected_counts, axis=0)

    return chi2


def split_by_block_sizes(std, block_sizes):
    std_listed = []
    start_idx = 0
    for size in block_sizes:
        end_idx = start_idx + size
        std_listed.append(std[start_idx:end_idx])
        start_idx = end_idx
    return std_listed


#### MISCELLANEOUS

def print_memory_usage():
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    print(f"RAM usage: {mem_info.rss / 1024 ** 2:.2f} MB")

def compute_element(p, masked_X, ready_Z, total_X_matrix):
    return masked_X[p] - (ready_Z[p] @ total_X_matrix)


def split_number(B, L):
    B_blocks = []
    start = 1
    end = L

    while start <= B:
        block = list(range(start, min(end + 1, B + 1)))
        B_blocks.append(block)
        start += L
        end += L

    return B_blocks


def get_size_in_gb(data):
    if isinstance(data, list):
        return sum([get_size_in_gb(item) for item in data])
    elif isinstance(data, np.ndarray):
        return data.nbytes / (1024 ** 3)
    elif issparse(data):
        return (data.data.nbytes + data.indices.nbytes + data.indptr.nbytes) / (1024 ** 3)
    else:
        raise TypeError("Unsupported data type")


def convert_seconds_to_dhms(seconds):
    days = seconds // 86400  
    remaining_seconds = seconds % 86400
    hours = remaining_seconds // 3600  
    remaining_seconds %= 3600
    minutes = remaining_seconds // 60  
    remaining_seconds %= 60

    return f"{days}d {hours}h {minutes}m {remaining_seconds}s"

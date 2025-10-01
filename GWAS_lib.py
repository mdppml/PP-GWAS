import numpy as np
from scipy.sparse import csr_matrix, vstack, hstack, save_npz, issparse, load_npz
import gc
import random
from scipy.sparse import bmat
import os
from scipy.stats import chisquare

np.random.seed(420)
def get_O_X(C, seed):
    A = np.random.rand(C + 1, C)
    Q, R = np.linalg.qr(A)
    return Q

def load_X(N, M, C, P, B, p, b, idx, X_list):
    try:
        loaded_data = load_npz('Data/N{}_M{}_C{}_P{}_B{}/Party_{}/X_block_{}.npz'.format(N, M, C, P, B, p, b))
        if loaded_data is not None and loaded_data.size > 0:
            X_list[idx] = loaded_data
        else:
            print(f"Warning: Loaded data for block {idx} is empty or None.")
            X_list[idx] = None
    except Exception as e:
        print(f"Error loading data for block {idx}: {e}")
        X_list[idx] = None

def get_O(N, K, seed):
    def generate_semi_orthogonal_matrix(rows, cols):
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
    np.random.seed(seed)
    min_block_size = 95
    max_block_size = 100
    num_blocks = N // max_block_size
    while num_blocks<=K:
        min_block_size -= 5
        max_block_size -= 5
        num_blocks = N // max_block_size
    remaining_columns = N % max_block_size
    blocks = []
    for i in range(K):
        block_size = np.random.randint(min_block_size, max_block_size + 1)
        block = generate_semi_orthogonal_matrix(block_size+1,block_size)
        blocks.append(block)
    for i in range(num_blocks - 1 - K):
        block_size = np.random.randint(min_block_size, max_block_size + 1)
        block = generate_random_orthogonal_matrix(block_size)
        blocks.append(block)
    last_block_size =max(N - sum(block.shape[1] for block in blocks),0)
    last_block = generate_random_orthogonal_matrix(last_block_size)
    blocks.append(last_block)
    O_blocks = []
    for i, block in enumerate(blocks):
        num_zero_rows = sum(block.shape[0] for block in blocks[:i])
        num_zero_cols = sum(block.shape[1] for block in blocks[:i])
        upper_padding = csr_matrix((num_zero_rows, N))
        lower_padding = csr_matrix((N + K - num_zero_rows - block.shape[0], N))
        left_padding = csr_matrix((block.shape[0], num_zero_cols))
        right_padding = csr_matrix((block.shape[0], N - num_zero_cols - block.shape[1]))
        padded_block = hstack((left_padding, block, right_padding)).tocsc()
        padded_block = vstack((upper_padding, padded_block, lower_padding)).tocsc()
        O_blocks.append(padded_block)
    O = random_permute_columns(sum(O_blocks),0)
    return O

def get_square_O(N, seed):
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

    O = bmat(block_matrix)
    O = random_permute_square_columns(O, seed)
    O = O.tocsr()
    O = O.astype(np.float32)
    return O


def random_permute_square_columns(matrix, seed):
    num_columns = matrix.shape[1]
    np.random.seed(seed)
    permuted_indices = np.random.permutation(num_columns)
    if not isinstance(matrix, np.ndarray):
        matrix = matrix.tocsr()
    permuted_matrix = matrix[:, permuted_indices]
    return permuted_matrix


def random_permute_columns(matrix, seed):
    num_columns = matrix.shape[1]
    np.random.seed(seed)
    permuted_indices = np.random.permutation(num_columns)
    permuted_matrix = matrix[:, permuted_indices]
    return permuted_matrix


def generate_random_orthogonal_matrix(size):
    random_matrix = np.random.randn(size, size)
    tol = 1e-7
    Q, R = np.linalg.qr(random_matrix)
    Q[abs(Q) < tol] = 0
    return Q

def standardization_normalizers(M, P, p, seed):
    np.random.seed(seed)
    matrix = np.random.randint(-10, 10, size=(P, M), dtype=np.int8)
    column_sum_vector = np.sum(matrix, axis=0, dtype=np.int32)
    return matrix[p - 1, :], column_sum_vector

def generate_O_b(M, B, p, P, N, C, block_size,b):
    if not os.path.exists('Data/N{}_M{}_C{}_P{}_B{}/Masks'.format(N, M, C, P, B)):os.makedirs('Data/N{}_M{}_C{}_P{}_B{}/Masks'.format(N, M, C, P, B),exist_ok=True)
    O = get_O(block_size, 0, b)
    filename = os.path.join('Data/N{}_M{}_C{}_P{}_B{}/Masks'.format(N, M, C, P, B),"O_b_block_{}.npz".format(b + 1))
    save_npz(filename, O)
    del O
    gc.collect()

def save_matrices(X,b,N,M,C,P,B):
    if not os.path.exists('Data/N{}_M{}_C{}_P{}_B{}/Server'.format(N, M, C, P, B)): os.makedirs('Data/N{}_M{}_C{}_P{}_B{}/Server'.format(N, M, C, P, B), exist_ok=True)
    np.save('Data/N{}_M{}_C{}_P{}_B{}/Server/X_block_{}.npy'.format(N, M, C, P, B,b+1),X)

def generate_O_X(M, B, p, P, N, C, block_size, additional_cols,b):
    directory = 'Data/N{}_M{}_C{}_P{}_B{}/Masks'.format(N, M, C, P, B)
    if not os.path.exists(directory):
        os.makedirs(directory, exist_ok=True)
    if additional_cols<0:
        print('additional_col = 0 found')
        additional_cols=0
    O = get_O(block_size, additional_cols, b).tocsc()
    filename = os.path.join(directory, "O_X_block_{}.npz".format(b + 1))
    save_npz(filename, O)
    del O
    gc.collect()

def generate_Z_mask(M, B, p, P, N, C):
    directory = 'Data/N{}_M{}_C{}_P{}_B{}/Masks'.format(N, M, C, P, B)
    if not os.path.exists(directory):
        os.makedirs(directory, exist_ok=True)
    if p == 1:
        O = get_O(N, 1, 0)
        save_npz(
            'Data/N{}_M{}_C{}_P{}_B{}/Masks/Z_mask.npz'.format(N, M, C, P, B),O)
        del O
        gc.collect()

def generate_O(M, B, K, p, P, N, C):
    if p == 1:
        O = get_O(N, 1, 1)
        save_npz(
            'Data/N{}_M{}_C{}_P{}_B{}/Masks/O.npz'.format(N, M, C, P, B),O)
        del O
        gc.collect()

def generate_a_number(seed):
    random.seed(seed)
    value = random.random() + 1
    return value

def hwe_chi_square(obs_counts):
    total = sum(obs_counts)
    p = (2 * obs_counts[0] + obs_counts[1]) / (2 * total)
    q = 1 - p
    expected_counts = np.array([total * p ** 2, 2 * total * p * q, total * q ** 2])
    chi2, _ = chisquare(f_obs=obs_counts, f_exp=expected_counts, axis=0)
    return chi2

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
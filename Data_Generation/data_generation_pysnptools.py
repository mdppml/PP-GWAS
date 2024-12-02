import numpy as np
from scipy.sparse import lil_matrix, save_npz
import multiprocessing as mp
from pysnptools.util import snp_gen
from scipy.sparse import csr_matrix
import gc
from scipy.stats import chisquare

def make_matrix(N, columns, seed):
    chr_count = 1
    snpdata = snp_gen(fst=.1, dfr=.25, iid_count=N, sid_count=columns, maf_low=.05, seed=seed, chr_count=chr_count)
    if (int(snpdata.iid_count))!= N:
        make_matrix(N,columns,seed+1500)
    return snpdata.valde


def hwe_chi_square(obs_counts):
    total = sum(obs_counts)
    p = (2 * obs_counts[0] + obs_counts[1]) / (2 * total)
    q = 1 - p
    expected_counts = np.array([total * p ** 2, 2 * total * p * q, total * q ** 2])
    chi2, _ = chisquare(f_obs=obs_counts, f_exp=expected_counts)
    chi2, _ = chisquare(f_obs=obs_counts, f_exp=expected_counts, axis=0)
    return chi2


def make_matrix2(N, columns, seed):
    chr_count = 1
    chunk_size = int(columns/3)
    total_columns = 0
    valid_columns = []
    X = None
    while len(valid_columns) < columns:
        x = 5000
        X_chunk = None
        while X_chunk is None or X_chunk.shape[0] < N:
            snpdata = snp_gen(fst=.1, dfr=.25, iid_count=N+x, sid_count=chunk_size, maf_low=.05, seed=seed,
                              chr_count=chr_count)
            X_chunk = snpdata.val
            x = x + 5000
            del snpdata
        X_chunk = X_chunk[:N, :]
        print('X_chunk ready')
        zero_counts = np.sum(X_chunk == 0, axis=0)
        one_counts = np.sum(X_chunk == 1, axis=0)
        two_counts = np.sum(X_chunk == 2, axis=0)
        maf = ((two_counts * 2) + one_counts) / (2 * N)
        obs_counts_array = np.vstack([zero_counts, one_counts, two_counts]).T
        chi2_values = np.apply_along_axis(hwe_chi_square, 1, obs_counts_array)

        valid_col_indices = np.where((chi2_values < 25) & (maf > 0.095))[0]
        if X is None:
            X = X_chunk[:, valid_col_indices]
        else:
            X = np.hstack((X, X_chunk[:, valid_col_indices]))

        valid_columns.extend(total_columns + valid_col_indices)

        if len(valid_columns) >= columns:
            break

        total_columns += chunk_size
        print(f'columns ready so far - {len(valid_columns)}')
        seed = seed + 1

    if len(valid_columns) > columns:
        X = X[:, :columns]

    return X.astype(np.float32)


def worker_block(args):
    j, N, M, C, P, B, vertical_split_size, horizontal_split_size, seed = args

    columns = horizontal_split_size if j != B - 1 else M - j * horizontal_split_size
    block = make_matrix2(N, columns, seed)

    for i in range(P):
        start_row = i * vertical_split_size
        end_row = start_row + vertical_split_size if i != P - 1 else N
        party_block = block[start_row:end_row, :]

        if isinstance(party_block, np.ndarray):
            party_block = csr_matrix(party_block)

        save_npz(f'../Code/Data/N{N}_M{M}_C{C}_P{P}_B{B}/Party_{i + 1}/X_block_{j + 1}.npz', party_block)
        print(f'Data generated for Party {i + 1} for block {j + 1}')

    del block

def generate_G(N, M, C, P, B):
    vertical_split_size = N // P
    horizontal_split_size = M // B
    base_seed = 0

    args_list = []
    for j in range(B):
        seed = base_seed + j
        args_list.append((j, N, M, C, P, B, vertical_split_size, horizontal_split_size, seed))

    with mp.Pool(processes=min(mp.cpu_count(), len(args_list))) as pool:
        pool.map(worker_block, args_list)

    print('Done')

import psutil

def print_memory_usage():
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    print(f"RAM usage: {mem_info.rss / 1024 ** 2:.2f} MB")


def generate_Z(N, C, P, filename_prefix):
    dtype = np.float32
    random_values = np.column_stack([np.random.normal(mean, std_dev, N) for mean, std_dev in zip([1.5]*C, [2]*C)])
    for i in range(P):
        start = i * (N // P)
        end = start + (N // P if i < P - 1 else N - (N // P) * (P - 1)) 
        np.save(f"{filename_prefix}/Party_{i+1}/Z.npy", random_values[start:end, :].astype(dtype))
    print('Generated Z')


def generate_y(N, P, directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

    y = np.random.normal(0.0125, 1, (N, 1))

    split_size = N // P

    for i in range(P):
        start = i * split_size
        end = (i + 1) * split_size if i != P - 1 else N   
        part = y[start:end, :]
        directory_path = os.path.join(directory, f'Party_{i + 1}')
        if not os.path.exists(directory_path):
            os.makedirs(directory_path)

        np.save(os.path.join(directory_path, 'y.npy'), part)
    print('Generated y')

import os
import time
import argparse
def main():
    print("Starting Data Generation")
    parser = argparse.ArgumentParser()
    parser.add_argument("--number_of_samples", type=int, required=True)
    parser.add_argument("--number_of_snps", type=int, required=True)
    parser.add_argument("--number_of_covariates", type=int, required=True)
    parser.add_argument("--number_of_parties", type=int, required=True)
    parser.add_argument("--number_of_blocks", type=int, required=True)
    args = parser.parse_args()
    P = args.number_of_parties
    N = args.number_of_samples
    M = args.number_of_snps
    C = args.number_of_covariates
    B = args.number_of_blocks
    tt=time.time()
    directory = {}
    for p in range(P):
        directory[p] = "/home/swaminathan/ppREGENIE/Data/N{}_M{}_C{}_P{}_B{}/Party_{}".format(N, M, C,
                                                                                                          P,
                                                                                                          B, p + 1)
        if not os.path.exists(directory[p]):
            os.makedirs(directory[p])
    generate_Z(N, C, P, '../Code/Data/N{}_M{}_C{}_P{}_B{}'.format(N, M, C, P, B))
    generate_y(N, P, '../Code/Data/N{}_M{}_C{}_P{}_B{}'.format(N, M, C, P, B))
    generate_G(N, M, C, P, B)
    print(f'time taken to run the program = {time.time() - tt}')

if __name__ == "__main__":
    main()


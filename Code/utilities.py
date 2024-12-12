import os
import psutil
import numpy as np
from scipy.sparse import issparse
import mkl

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

def print_memory_usage():
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    print(f"RAM usage: {mem_info.rss / 1024 ** 2:.2f} MB")

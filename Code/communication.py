
import numpy as np
from scipy.sparse import csr_matrix, csc_matrix
import pickle
import struct
import time
from threading import Thread
from concurrent.futures import ThreadPoolExecutor
import mkl

CHUNK_SIZE = 60000

#####################################################################################################
#--------------------------------------------SENDING DATA-------------------------------------------#
#####################################################################################################

def send_data_to_server(party_socket, matrix, party_id):
    data = serialize_matrix_party(party_id, matrix)
    data_size = struct.pack('!Q', len(data))
    party_socket.sendall(data_size)
    max_chunk_size = 8192
    total_sent = 0
    while total_sent < len(data):
        end = min(total_sent + max_chunk_size, len(data))
        chunk = data[total_sent:end]
        party_socket.sendall(chunk)
        total_sent += len(chunk)

def send_data_to_party(party_socket, matrix):
    data = serialize_matrix_server(matrix)
    data_size = struct.pack('!Q', len(data))
    party_socket.sendall(data_size)
    chunk_count = 0
    for i in range(0, len(data), CHUNK_SIZE):
        chunk = data[i:i + CHUNK_SIZE]
        party_socket.sendall(chunk)
        chunk_count += len(chunk)

def send_data_to_parties(party_sockets,matrices):
    mkl.set_num_threads(4)
    with ThreadPoolExecutor(max_workers=len(party_sockets)) as executor:
        futures = []
        for p in range(len(party_sockets)):
            futures.append(executor.submit(send_data_to_party, party_sockets[p], matrices[p]))
        for future in futures:
            future.result()
    mkl.set_num_threads(16)

#####################################################################################################
#-------------------------------------------SERIALIZATION-------------------------------------------#
#####################################################################################################

def serialize_matrix_party(party_id, matrix):
    def serialize_1d_array(array):
        return array.astype(np.float32).tobytes()
    if isinstance(matrix, np.ndarray):
        rows, cols = matrix.shape if matrix.ndim > 1 else (len(matrix), 0)
        matrix_type = (
            0 if cols == 0 else
            2 if rows == 1 else
            1 if cols == 1 else
            3
        )
        matrix_data = serialize_1d_array(matrix)
    elif isinstance(matrix, (csr_matrix, csc_matrix)):
        rows, cols = matrix.shape
        matrix_type = 4 if isinstance(matrix, csr_matrix) else 5
        matrix_data = pickle.dumps(matrix)
    elif isinstance(matrix, list) and all(isinstance(x, np.ndarray) and x.ndim == 1 for x in matrix):
        rows, cols = len(matrix), 0
        matrix_type = 6
        serialized_matrices = [serialize_matrix_party(party_id, mat) for mat in matrix]
        matrix_data = b''.join(serialized_matrices)
    else:
        raise ValueError("Unsupported matrix type or invalid list elements.")
    header = struct.pack('!4I1Q', party_id, rows, cols, matrix_type, len(matrix_data))
    return header + matrix_data

def serialize_matrix_server(matrix):
    return pickle.dumps(matrix)

#####################################################################################################
#-------------------------------------------RECEIVING DATA------------------------------------------#
#####################################################################################################


def receive_data(party_socket):
    tt=time.time()
    data_size_bytes = party_socket.recv(8)
    if len(data_size_bytes) != 8:
        raise ConnectionError("party disconnected prematurely while sending data size.")
    tt=time.time()
    data_size = struct.unpack('!Q', data_size_bytes)[0]
    tt=time.time()
    chunks = []
    received_size = 0
    while received_size < data_size:
        chunk = party_socket.recv(min(16384, data_size - received_size))
        if not chunk:
            break
        chunks.append(chunk)
        received_size += len(chunk)
    if received_size < data_size:
        raise ConnectionError(f"Expected {data_size} bytes but received only {received_size} bytes.")
    data = b''.join(chunks)
    return data

def receive_and_sum(party_sockets):
    def threaded_receive_data(party_socket, results, index):
        results[index] = receive_data(party_socket)

    mkl.set_num_threads(1)
    raw_data = [None] * len(party_sockets)
    threads = []

    for i, party_socket in enumerate(party_sockets):
        thread = Thread(target=threaded_receive_data, args=(party_socket, raw_data, i))
        thread.start()
        threads.append(thread)

    for thread in threads:
        thread.join()

    mkl.set_num_threads(16)
    masked_X = [None] * len(raw_data)

    for data in raw_data:
        party_id, matrix = deserialize_matrix_server(data)
        masked_X[party_id - 1] = matrix
    del raw_data
    total_X_matrix = np.sum(masked_X, axis=0)
    del masked_X
    return total_X_matrix

def receive_sum_and_append(party_sockets):
    def threaded_receive_data(party_socket, results, index):
        results[index] = receive_data(party_socket)

    mkl.set_num_threads(1)
    raw_data = [None] * len(party_sockets)
    threads = []

    for i, party_socket in enumerate(party_sockets):
        thread = Thread(target=threaded_receive_data, args=(party_socket, raw_data, i))
        thread.start()
        threads.append(thread)

    for thread in threads:
        thread.join()

    mkl.set_num_threads(16)
    masked_X = [None] * len(raw_data)

    for data in raw_data:
        party_id, matrix = deserialize_matrix_server(data)
        masked_X[party_id - 1] = matrix
    del raw_data
    total_X_matrix = np.sum(masked_X, axis=0)
    return total_X_matrix, masked_X

def receive_data_from_server(server_socket):
    data_size_bytes = server_socket.recv(8)
    data_size = struct.unpack('!Q', data_size_bytes)[0]
    received_data = bytearray()
    while len(received_data) < data_size:
        chunk = server_socket.recv(min(CHUNK_SIZE, data_size - len(received_data)))
        if not chunk:
            raise ConnectionError("Connection closed before all data was received.")
        received_data.extend(chunk)
    matrix = deserialize_matrix_party(received_data)
    return matrix

#####################################################################################################
#-------------------------------------------DESERIALIZATION-------------------------------------------#
#####################################################################################################

def deserialize_matrix_server(serialized_data):
    header_length = 4 * 4 + 1 * 8
    party_id, rows, cols, matrix_type, data_length = struct.unpack('!4I1Q', serialized_data[:header_length])
    matrix_data = serialized_data[header_length:]
    del serialized_data
    if matrix_type == 0:
        matrix = np.frombuffer(matrix_data, dtype=np.float32)
    elif matrix_type in [1, 2, 3]:
        matrix = np.frombuffer(matrix_data, dtype=np.float32).reshape(rows, cols if cols else 1)
    elif matrix_type == 4:
        matrix = pickle.loads(matrix_data)
    elif matrix_type == 5:
        matrix = pickle.loads(matrix_data)
    elif matrix_type == 6:
        matrix = []
        offset = 0
        for _ in range(rows):
            sub_party_id, sub_rows, _, sub_matrix_type, sub_data_length = struct.unpack(
                '!4I1Q', matrix_data[offset:offset + header_length]
            )
            offset += header_length
            sub_matrix_data = matrix_data[offset:offset + sub_data_length]
            offset += sub_data_length
            matrix.append(np.frombuffer(sub_matrix_data, dtype=np.float32))
    else:
        raise ValueError("Unsupported matrix type.")
    return party_id, matrix

def deserialize_matrix_party(data):
    return pickle.loads(data)

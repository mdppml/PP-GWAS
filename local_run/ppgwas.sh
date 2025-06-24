#!/usr/bin/env bash

set -euo pipefail

BASE_PORT=$1
N=$2
M=$3
C=$4
P=$5
B=$6
FOLDS=$7
BPR=$8

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_ROOT="logs/ppgwas_${TIMESTAMP}"
SERVER_LOG="${LOG_ROOT}/server"
CLIENTS_LOG="${LOG_ROOT}/clients"

mkdir -p "${SERVER_LOG}"
mkdir -p "${CLIENTS_LOG}"

rm -f Data/server_ready_*.txt Data/ip_address_file.txt

eval "$(conda shell.bash hook)"
conda activate ppREGENIE

echo "Launching server on port ${BASE_PORT} (logs → ${SERVER_LOG})"
python -u server.py \
  --number_of_parties "${P}" \
  --base_port       "${BASE_PORT}" \
  --number_of_samples   "${N}" \
  --number_of_snps      "${M}" \
  --number_of_covariates "${C}" \
  --number_of_blocks    "${B}" \
  --number_of_folds     "${FOLDS}" \
  --number_of_blocks_per_run "${BPR}" \
  > "${SERVER_LOG}/output.txt" 2> "${SERVER_LOG}/error.txt" &

sleep 4

echo "Spawning ${P} clients (logs → ${CLIENTS_LOG})"
for (( i=1; i<=P; i++ )); do
  THIS_LOG="${CLIENTS_LOG}/party${i}"
  mkdir -p "${THIS_LOG}"

  python -u client.py \
    --number_of_parties "${P}" \
    --party_id           "${i}" \
    --base_port          "${BASE_PORT}" \
    --number_of_samples   "${N}" \
    --number_of_snps      "${M}" \
    --number_of_covariates "${C}" \
    --number_of_blocks    "${B}" \
    --number_of_folds     "${FOLDS}" \
    --number_of_blocks_per_run "${BPR}" \
    > "${THIS_LOG}/output.txt" 2> "${THIS_LOG}/error.txt" &

  sleep 0.1
done

wait

echo "All processes finished. Logs in ${LOG_ROOT}"

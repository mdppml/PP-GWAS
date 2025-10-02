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
QT_XCB_GL_INTEGRATION=none

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_ROOT="logs/ppgwas_${TIMESTAMP}"
SERVER_LOG="${LOG_ROOT}/server"
CLIENTS_LOG="${LOG_ROOT}/clients"
RUN_LOG="${LOG_ROOT}/run.log"
RUN_ERR="${LOG_ROOT}/run.err"

mkdir -p "${SERVER_LOG}" "${CLIENTS_LOG}"

rm -f Data/server_ready_*.txt Data/ip_address_file.txt

exec > >(tee -a "${RUN_LOG}") 2> >(tee -a "${RUN_ERR}" >&2)

if command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)"
  conda activate ppgwas_test
fi

echo "Please refer to the following folders to read outputs and errors from both the server and the clients."
echo "Server output: ${SERVER_LOG}/output.txt"
echo "Server error: ${SERVER_LOG}/error.txt"
echo "Clients output: ${CLIENTS_LOG}/party*/output.txt"
echo "Clients error: ${CLIENTS_LOG}/party*/error.txt"
echo "Run output: ${RUN_LOG}"
echo "Run error: ${RUN_ERR}"

echo "Launching server on port ${BASE_PORT}"
python -u server.py \
  --number_of_parties "${P}" \
  --base_port "${BASE_PORT}" \
  --number_of_samples "${N}" \
  --number_of_snps "${M}" \
  --number_of_covariates "${C}" \
  --number_of_blocks "${B}" \
  --number_of_folds "${FOLDS}" \
  --number_of_blocks_per_run "${BPR}" \
  > "${SERVER_LOG}/output.txt" 2> "${SERVER_LOG}/error.txt" &

SERVER_PID=$!

sleep 4

echo "Spawning ${P} clients"
PIDS=("${SERVER_PID}")
NAMES=("Server")
for (( i=1; i<=P; i++ )); do
  THIS_LOG="${CLIENTS_LOG}/party${i}"
  mkdir -p "${THIS_LOG}"
  python -u client.py \
    --number_of_parties "${P}" \
    --party_id "${i}" \
    --base_port "${BASE_PORT}" \
    --number_of_samples "${N}" \
    --number_of_snps "${M}" \
    --number_of_covariates "${C}" \
    --number_of_blocks "${B}" \
    --number_of_folds "${FOLDS}" \
    --number_of_blocks_per_run "${BPR}" \
    > "${THIS_LOG}/output.txt" 2> "${THIS_LOG}/error.txt" &
  PIDS+=("$!")
  NAMES+=("Client_${i}")
  sleep 0.1
done

FAIL=0
for idx in "${!PIDS[@]}"; do
  pid="${PIDS[$idx]}"
  name="${NAMES[$idx]}"
  set +e
  wait "$pid"
  rc=$?
  set -e
  echo "${name} exited with ${rc}"
  if [[ $rc -ne 0 ]]; then
    FAIL=1
  fi
done

if [[ $FAIL -ne 0 ]]; then
  echo "At least one process failed."
  echo "Refer to ${SERVER_LOG}/error.txt and ${CLIENTS_LOG}/party*/error.txt for errors."
  exit 1
fi

echo "All processes finished. Logs in ${LOG_ROOT}."

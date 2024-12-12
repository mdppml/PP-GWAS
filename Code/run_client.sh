#!/bin/bash

client_id=$8
DIR=${10}

mkdir -p $DIR

client_dir="${DIR}/client${client_id}"
mkdir -p $client_dir

echo "Starting job script for client $client_id" > ${client_dir}/output.txt 2> ${client_dir}/error.txt

eval "$(conda shell.bash hook)" >> ${client_dir}/output.txt 2>> ${client_dir}/error.txt
echo "Activating conda environment for client $client_id" >> ${client_dir}/output.txt 2>> ${client_dir}/error.txt
conda activate "INSERT_CONDA_ENVIRONMENT_HERE" >> ${client_dir}/output.txt 2>> ${client_dir}/error.txt

python -u client.py --number_of_clients $7 --party_id $8 --base_port $1 --number_of_samples $2 --number_of_snps $3 --number_of_covariates $4 --number_of_blocks $5 --number_of_folds $6 --number_of_blocks_per_run $9 >> ${client_dir}/output.txt 2>> ${client_dir}/error.txt

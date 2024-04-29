#!/bin/bash

echo "Starting job script for server"

eval "$(conda shell.bash hook)"
echo "Activating conda environment for server"
conda activate /home/swaminathan/miniconda3/envs/ppregenie_env

python -u server.py --number_of_clients $7 --base_port $1 --number_of_samples $2 --number_of_snps $3 --number_of_covariates $4 --number_of_blocks $5 --number_of_folds $6 --number_of_blocks_per_run $8

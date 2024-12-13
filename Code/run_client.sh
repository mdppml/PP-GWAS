#!/bin/bash

party_id=$9
DIR=${10}

mkdir -p $DIR

party_dir="${DIR}/party${party_id}"
mkdir -p $party_dir

echo "Starting job script for party $party_id" > ${party_dir}/output.txt 2> ${party_dir}/error.txt

eval "$(conda shell.bash hook)" >> ${party_dir}/output.txt 2>> ${party_dir}/error.txt
echo "Activating conda environment for party $party_id" >> ${party_dir}/output.txt 2>> ${party_dir}/error.txt
conda activate "INSERT_CONDA_ENVIRONMENT_HERE" >> ${party_dir}/output.txt 2>> ${party_dir}/error.txt

python -u party.py --number_of_partys $7 --party_id $9 --base_port $1 --number_of_samples $2 --number_of_snps $3 --number_of_covariates $4 --number_of_blocks $5 --number_of_folds $6 --number_of_blocks_per_run $8 >> ${party_dir}/output.txt 2>> ${party_dir}/error.txt

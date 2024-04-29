#!/bin/bash
#SBATCH --job-name=data_gen_pysnp
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100000
#SBATCH --time=06:00:00
#SBATCH --partition=day
#SBATCH --output=job_output_%j.txt


echo "Starting job script"

eval "$(conda shell.bash hook)"

echo "Activating conda environment"
conda activate 'enter_conda_directory_here'

echo "Running data_generation.sh script"
python -u data_generation_pysnptools.py --number_of_samples $1 --number_of_snps $2 --number_of_covariates $3 --number_of_parties $4 --number_of_blocks $5
wait

echo "Job script finished"









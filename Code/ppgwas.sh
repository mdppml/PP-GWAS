#!/bin/bash
#SBATCH --job-name=ppREGENIE
#SBATCH --nodes=4  
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=24
#SBATCH --time=03:59:00
#SBATCH --partition=day
#SBATCH --mem=240000
#SBATCH --output=test_%j.txt
# Define the log folder with the job ID
client_DIR="logs/job_$SLURM_JOB_ID/"
server_DIR="logs/job_$SLURM_JOB_ID/server"
mkdir -p $client_DIR
mkdir -p $server_DIR

rm -f ../test_site/Data/server_ready_*.txt
rm ../test_site/Data/ip_address_file.txt
chmod +x run_server.sh
chmod +x run_client.sh

export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun --exclusive -u --nodes=1 --ntasks=1 run_server.sh $1 $2 $3 $4 $5 $6 $7 $8 > $server_DIR/output.txt 2> $server_DIR/error.txt &

sleep 4

for i in $(seq 1 $7)
do
  srun --exclusive -u --nodes=1 --ntasks=1 run_client.sh $1 $2 $3 $4 $5 $6 $7 $8 $i $client_DIR & 
  sleep 0.1
done

wait

echo "Job script finished"

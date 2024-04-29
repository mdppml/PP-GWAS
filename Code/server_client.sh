#!/bin/bash
#SBATCH --job-name=PPF-GWAS
#SBATCH --nodes=7
#SBATCH --ntasks=7
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G 
#SBATCH --time=03:00:00
#SBATCH --partition=day
#SBATCH --output=job_output_%j.txt
rm -f /home/swaminathan/ppREGENIE/Data/server_ready_*.txt
rm /home/swaminathan/ppREGENIE/Data/ip_address_file.txt
chmod +x run_server.sh
chmod +x run_client.sh

export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun --exclusive -u --nodes=1 --ntasks=1 run_server.sh $1 $2 $3 $4 $5 $6 $7 $8 &

sleep 4

for i in $(seq 1 $7)
do
  srun --exclusive -u --nodes=1 --ntasks=1 run_client.sh $1 $2 $3 $4 $5 $6 $7 $i $8 &
  sleep 0.15
done

wait

echo "Job script finished"


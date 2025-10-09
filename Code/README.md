# 🧬🔒 PP-GWAS 

## 📂 Files and Structure

#### SLURM Script: `ppgwas.sh`
This SLURM batch script orchestrates the execution of a distributed privacy-preserving GWAS workflow by allocating resources from the cluster, starting the server, and launching client processes on separate nodes while managing the necessary environment and logging. It launches both `run_server.sh` and `run_client.sh`. 

#### Bash Script: `run_server.sh`
This script sets up the server environment, activates the appropriate Conda environment, and executes the server.py script with the necessary command-line arguments for the workflow.

#### Python Script: `server.py`
This Python file implements the server-side logic for the PP-GWAS framework. It establishes socket-based connections with multiple clients, coordinates data exchange, and executes the distributed Alternating Direction Method of Multipliers (ADMM) and distributed Conjugate Gradient Descent (CGD) algorithms. It then performs SNP association testing and outputs the negative logarithm of p-values for each SNP to the file `../test_site/Data/N{N}_M{M}_C{C}_P{P}_B{B}/neg_log_transfer.npy`.

#### Bash Script: `run_client.sh`
This script sets up the client environments, activates the appropriate Conda environment, and executes multiple instances of client.py script with the necessary command-line arguments for the workflow.

#### Python Script: `client.py`
This Python script implements the client-side logic of PP-GWAS. Each client connects to the server, processes local genomic data, and coordinates performing the distributed computations. 

#### Python Script: `GWAS_lib.py`
This script contains essential functions for genomic data preprocessing and manipulation, including generating and loading randomized encodings, handling sparse matrices, performing standardizations, and calculating metrics such as sparsity and Hardy-Weinberg equilibrium chi-square values. It also includes utility functions for splitting data into blocks and randomizing matrices.

#### Python Script: `communication.py`
This script provides functions for serializing, deserializing, sending, and receiving data between the server and clients. It supports various data types, handles multi-threaded data transfer, and ensures efficient communication using optimizations like chunking and parallel processing.

#### Python Script: `utilities.py`
This script provides utility functions for system monitoring and data management.

#### Conda Environments File: `ppgwas_environment.yml`
This file specifies the Conda environment required for running the Python script. It includes dependencies optimized for systems with Intel chips, particularly leveraging MKL (Math Kernel Library) for efficient sparse matrix computations. It should be used to create a Conda environment for running the scripts.

#### Pip Requirements File: `ppgwas_pip.txt`
This file lists additional Python packages that should be installed using pip within the Conda environment.

---

## ⚙️ Requirements

### 📦 Conda Environment
The SLURM script requires a Conda environment. Replace `"INSERT_CONDA_ENVIRONMENT_HERE"` in the scripts `ppgwas.sh`, `run_server.sh`, and `run_client.sh` with the path to your Conda environment.

Intel MKL Requirement:
This code requires Intel MKL (Math Kernel Library) for sparse matrix computations, which is typically pre-installed on Linux-based clusters with Intel processors. If your system uses AMD chips or does not have MKL, consider modifying dependencies to use OpenBLAS as an alternative.

Ensure the `run_server.sh` and `run_client.sh` scripts are executable:
   ```bash
   chmod +x run_server.sh
   chmod +x run_client.sh
```

---

## ▶️ Usage

#### Command-Line Arguments
The slurm script requires the following arguments:
- `--base-port` (int): Base port for communication..
- `--number_of_samples` (int): Number of samples (N).
- `--number_of_snps` (int): Number of SNPs (M).
- `--number_of_covariates` (int): Number of covariates (C).
- `--number_of_blocks` (int): Number of blocks to split the data into (B).
- `--number_of_folds` (int): Number of folds to use for k-fold CV (k).
- `--number_of_parties` (int): Number of computational nodes (P).
- `--number_of_blocks_per_run` (int): Number of blocks to process at once. Depends on memory allocation. 


#### Example Command
## 🖥️ Using SLURM

Please edit the `ppgwas.sh` file with appropriate SLURM parameters based on your HPC configuration. 

```bash
sbatch ppgwas.sh 8110 1000 10000 5 2 5 3 2 
```
### SLURM Argument Template
The format for SLURM submission is as follows:
```bash
sbatch ppgwas.sh base_port number_of_samples number_of_snps number_of_covariates number_of_blocks number_of_folds number_of_parties number_of_blocks_per_run 
```

## 🛠️ Running without SLURM
If you're not using SLURM, you need to manually start the server and client scripts. First, run the server:
```bash
bash run_server.sh 8110 1000 10000 5 2 5 3 2
```
### Server Argument Template
The format for starting the server manually is:
```bash
bash run_server.sh base_port number_of_samples number_of_snps number_of_covariates number_of_blocks number_of_folds number_of_parties number_of_blocks_per_run
```
Next, activate the computational nodes by calling `run_client.sh`. For each computational node (party), start a separate instance of `run_client.sh`:
```bash
bash run_client.sh 8110 1000 10000 5 2 5 3 2 1 "../test_site/"
bash run_client.sh 8110 1000 10000 5 2 5 3 2 2 "../test_site/"
bash run_client.sh 8110 1000 10000 5 2 5 3 2 3 "../test_site/"

```
### Client Argument Template
The format for starting each client manually is:
```bash
bash run_client.sh base_port number_of_samples number_of_snps number_of_covariates number_of_blocks number_of_folds number_of_blocks_per_run number_of_parties party_id output_folder
```




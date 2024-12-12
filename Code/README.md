# PP-GWAS
# UNDER CONSTRUCTION

## Files and Structure

### Python Scripts:

### Bash Scripts: 
---

## Requirements

### Dependencies
Ensure the following Python libraries are installed:
- `numpy`
- `scipy`
- `mkl`
- `psutil`


### Conda Environment
The SLURM script requires a Conda environment. Replace `"INSERT_CONDA_ENVIRONMENT_HERE"` in the scripts `ppgwas.sh`, `run_server.sh`, and `run_client.sh` with the path to your Conda environment.

Ensure the `run_server.sh` and `run_client.sh` scripts are executable:
   ```bash
   chmod +x run_server.sh
   chmod +x run_client.sh
```

---

## Usage

#### Command-Line Arguments
The slurm script requires the following arguments:
- - `--base-port` (int): Base port for communication..
- `--number_of_samples` (int): Number of samples (N).
- `--number_of_snps` (int): Number of SNPs (M).
- `--number_of_covariates` (int): Number of covariates (C).
- `--number_of_blocks` (int): Number of blocks to split the data into (B).
- `--number_of_folds` (int): Number of folds to use for k-fold CV (k).
- `--number_of_parties` (int): Number of computational nodes (P).
- `--number_of_blocks_per_run` (int): Number of blocks to process at once. Depends on memory allocation. 


#### Example Command
```bash
sbatch ppgwas.sh 8110 1000 10000 5 2 5 3 2 
```
```bash
sbatch ppgwas.sh base_port number_of_samples number_of_snps number_of_covariates number_of_blocks number_of_folds number_of_parties number_of_blocks_per_run 
```

If you're not working with slurm, you must activate multiple instances of `run_client.sh` depending on the number of computational nodes apart from the server. 

```bash
python run_server.sh 8110 1000 10000 5 2 5 3 2
```
```bash
python run_server.sh base_port number_of_samples number_of_snps number_of_covariaets number_of_blocks number_of_folds number_of_parties number_of_blocks_per_run
```

```bash
python run_client.sh 8110 1000 10000 5 2 5 3 p 2
```
```bash
python run_client.sh base_port number_of_samples number_of_snps number_of_covariaets number_of_blocks number_of_folds number_of_parties number_of_blocks_per_run party_id 
```




---

## Output Structure
The generated data will be saved in the following directory structure:
```
../Code/Data/N{N}_M{M}_C{C}_P{P}_B{B}/
    ├── Party_1/
    │   ├── X_block_1.npz
    │   ├── ...
    │   ├── X_block_B.npz
    │   ├── y.npy
    │   └── Z.npy
    ├── Party_2/
    │   ├── X_block_1.npz
    │   ├── ...
    │   ├── X_block_B.npz
    │   ├── y.npy
    │   └── Z.npy
    ├── ...
    └── Party_P/
        ├── X_block_B.npz
        └── y.npy
```
### Files Generated:
- **`X_block_{j}.npz`**: Genomic data each block
- **`Z.npy`**: Covariate data.
- **`y.npy`**: Phenotype data.

---

## Notes
- Ensure sufficient memory and disk space to handle large datasets.
- The SLURM script assumes a cluster environment with Conda installed.
- Modify paths in both scripts as per your directory structure.

---


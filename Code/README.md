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
The SLURM script requires a Conda environment. Replace `"INSERT_CONDA_ENVIRONMENT_HERE"` in the script with the path to your Conda environment.

---

## Usage

#### Command-Line Arguments
The Python script requires the following arguments:
- `--number_of_samples` (int): Number of samples (N).
- `--number_of_snps` (int): Number of SNPs (M).
- `--number_of_covariates` (int): Number of covariates (C).
- `--number_of_parties` (int): Number of data-generating parties (P).
- `--number_of_blocks` (int): Number of SNP blocks (B).

#### Example Command
```bash
python data_generation_pysnptools.py --number_of_samples 100000 --number_of_snps 10000 --number_of_covariates 5 --number_of_parties 3 --number_of_blocks 2
```

The SLURM script accepts the same arguments as the Python script. Example:
```bash
sbatch data_generation_pysnptools_slurm.sh 100000 10000 5 3 2
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


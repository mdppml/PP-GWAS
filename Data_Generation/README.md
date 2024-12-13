# Synthetic Genomic Data Generation

## Files and Structure

### Python Script: `data_generation_pysnptools.py`
This Python script generates synthetic genomic data using the `pysnptools` library. The script supports parallel processing and allows for flexible configurations such as the number of samples, SNPs, covariates, computational nodes, and data blocks.

### Bash Script: `data_generation_pysnptools_slurm.sh`
A SLURM batch script to automate the execution of `data_generation_pysnptools.py` in a high-performance computing environment. The script specifies resource requirements, activates a Conda environment, and runs the Python script with user-defined arguments.

### Conda Environments File: `data_generation_environment.yml`
This file specifies the Conda environment required for running the Python script. It includes dependencies optimized for systems with Intel chips, particularly leveraging MKL (Math Kernel Library) for efficient sparse matrix computations. It should be used to create a Conda environment for running the scripts.

### Pip Requirements File: `data_generation_environment_pip.txt`
This file lists additional Python packages that should be installed using pip within the Conda environment. These packages extend functionality for genomic data processing.

---

## Requirements

### Dependencies
Ensure the following Python libraries are installed:
- `numpy`
- `scipy`
- `pysnptools`
- `psutil`
- `argparse`

### Conda Environment
The SLURM script requires a Conda environment. Replace `'enter_conda_directory_here'` in the script with the path to your Conda environment.

---

## Usage

#### Command-Line Arguments
The Python script requires the following arguments:
- `--number_of_samples` (int): Number of samples (N).
- `--number_of_snps` (int): Number of SNPs (M).
- `--number_of_covariates` (int): Number of covariates (C).
- `--number_of_parties` (int): Number of computational nodes (P).
- `--number_of_blocks` (int): Number of blocks to split the data into (B).

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


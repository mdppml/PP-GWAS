# Synthetic Genomic Data Generation

## Files and Structure

### Python Script: `data_generation_pysnptools.py`
This Python script generates synthetic genomic data using the `pysnptools` library and other utilities. The script supports parallel processing and allows for flexible configurations such as the number of samples, SNPs, covariates, parties, and data blocks.

### Bash Script: `data_generation_pysnptools_slurm.sh`
A SLURM batch script to automate the execution of `data_generation_pysnptools.py` in a high-performance computing environment. The script specifies resource requirements, activates a Conda environment, and runs the Python script with user-defined arguments.

---

## Requirements

### Dependencies
Ensure the following Python libraries are installed:
- `numpy`
- `scipy`
- `pysnptools`
- `multiprocessing`
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

#### SLURM Script Parameters
- `--job-name`: Name of the SLURM job.
- `--ntasks`: Number of tasks (always set to 1).
- `--cpus-per-task`: Number of CPUs for parallel processing.
- `--mem`: Memory allocation for the job.
- `--time`: Maximum job runtime.
- `--partition`: SLURM partition to use.

---

## Output Structure
The generated data will be saved in the following directory structure:
```
../Code/Data/N{N}_M{M}_C{C}_P{P}_B{B}/
    ├── Party_1/
    │   ├── X_block_1.npz
    │   ├── ...
    │   └── Z.npy
    │
    ├── Party_2/
    │   ├── X_block_1.npz
    │   ├── ...
    │   └── Z.npy
    │
    └── Party_P/
        ├── X_block_B.npz
        └── y.npy
```
### Files Generated:
- **`X_block_{j}.npz`**: SNP data matrix for each block and party.
- **`Z.npy`**: Covariate matrix for each party.
- **`y.npy`**: Response vector for each party.

---

## Notes
- Ensure sufficient memory and disk space to handle large datasets.
- The SLURM script assumes a cluster environment with Conda installed.
- Modify paths in both scripts as per your directory structure.

---


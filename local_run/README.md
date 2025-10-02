# 🧬🔒 PP-GWAS Local Toy Example

A self-contained, small-scale example of running PP-GWAS (Privacy Preserving Multi-Site Genome-wide Association Studies) locally without MKL. You can either use a Jupyter notebook or a one-click GUI script.

## 📂 What’s inside

- `sample_notebook.ipynb` for end-to-end experimentation
- `ppgwas_oneclick.py` for a GUI “one-click” run that displays results

## ⚙️ Requirements

1. Conda
2. bash

## 🛠 Installation

### 1) Get the code

```bash
git clone https://github.com/mdppml/PP-GWAS.git
cd PP-GWAS/local_run
```

### 2) Create and activate the Conda environment
```bash
conda env create -f environment.yml
conda activate ppgwas_test
```

### 3) Make the helper script executable
```bash
chmod +x ppgwas.sh
```

## ▶️ Choose how to run

### Option A — Jupyter Notebook

1) Launch Jupyter and select the kernel `ppgwas_test`  
   In Jupyter: `Kernel → Change Kernel → ppgwas_test`

2) Open and run: `sample_notebook.ipynb`

This notebook shows how to:
- Generate synthetic data
- Run the PP-GWAS pipeline
- Visualize results (including the Manhattan plot)

### Option B — One-Click App (GUI)

After activating the environment, run:

```bash
python3 ppgwas_oneclick.py
```

A window will open. Enter the parameters in the input boxes and press Run (you may leave BPR empty and it should always be less than or equal to B). The same window will display the Manhattan plot when the computation finishes.

## 🔧 Synthetic Data (defaults used in the notebook)

Recommended minimum/maximum for reproducibility:
- Samples `N`: 1,000/5,000
- Total SNPs `M`: 5,000/20,000
- Blocks `B`: 2/4
- Number of Nodes `P`: 2/6
- Covariates `C`: 2/10
- BPR `Blocks per Run`: <=B (Based on memory availability. Set to B by default.)
  
Notebook defaults:

```python
N (Number of Samples) = 10000
M (Number of SNPs) = 15000
B (Number of Blocks) = 2
C (Number of Covariates) = 10
P (Number of Computational Nodes) = 4
```

You can scale these up for larger experiments depending on available compute. If you encounter memory issues, reduce BPR (this determines how many blocks are processed at a time).

## ❓ Troubleshooting

- Ensure the `ppgwas_test` Conda environment is activated before running either the notebook or the one-click app.
- If `python3` isn’t found, try `python`.
- `Windows` users: please install Ubuntu on WSL to run `ppgwas.sh`.

# PP-GWAS Local Toy Example

A self-contained, small-scale example of running PP-GWAS (Privacy-Preserving Genome-Wide Association Studies) locally without MKL. You can either use a Jupyter notebook or a one-click GUI script.

## ✨ What’s inside

- `sample_notebook.ipynb` for end-to-end experimentation
- `ppgwas_oneclick.py` for a GUI “one-click” run that displays results
- Reproducible synthetic-data generation and visualization (including a Manhattan plot)

## ⚙️ Requirements

1. Conda
2. bash

## 🛠 Installation

Create and activate the environment, then make the helper script executable:

```bash
conda env create -f environment.yml
conda activate ppgwas_test
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

A window will open. Enter the parameters in the input boxes and press Run. The same window will display the Manhattan plot when the computation finishes.

## 🔧 Synthetic Data (defaults used in the notebook)

Recommended minimum for reproducibility:
- Samples `N`: 1,000
- Total SNPs `B`: 10,000
- SNPs per block: 5,000

Notebook defaults:

```python
N = 10000
B = 15000
Blocks = 2
C = 10
P = 4
```

You can scale these up for larger experiments depending on available compute. If you encounter memory issues, reduce batch or block sizes (for example, lower any BPR-like batching parameter in your setup).

## ✅ Outputs

- Notebook: cells will generate figures and outputs inline.
- One-Click App: the Manhattan plot is rendered in the application window after the run completes.

## ❓ Troubleshooting

- Ensure the `ppgwas_test` Conda environment is activated before running either the notebook or the one-click app.
- If `python3` isn’t found, try `python`.
- If Jupyter doesn’t show the `ppgwas_test` kernel, run `python -m ipykernel install --user --name ppgwas_test` inside the environment, then restart Jupyter.

## PP-GWAS Local Toy Example

This repository provides a self-contained, small-scale example of running PP-GWAS (Privacy-Preserving Genome-Wide Association Studies) locally without MKL dependencies. The included `sample_notebook.ipynb` demonstrates how to generate synthetic data, execute the analysis, and visualize the results.


### ⚙️ Requirements

1. **Conda** 
2. **bash** 

---

### 🛠 Installation

1. **Create the Conda environment**

   ```bash
   conda env create -f environment.yml
   ```

2. **Activate the environment**

   ```bash
   conda activate ppgwas_test
   ```

3. **Make the script executable**

   ```bash
   chmod +x ppgwas.sh
   ```

After creating the environment, select ppgwas_test as the Jupyter kernel `(Kernel → Change Kernel → ppgwas_test)` or set it as your IDE interpreter before running `sample_notebook.ipynb`. 

---


### 🔧 Synthetic Data Generation

For reproducibility, we recommend generating a dataset with at least:

* **Samples (N):** 1,000
* **Total SNPs (B):** 10,000
* **SNPs per blocks:** 5,000

In the notebook, the values are set as follows:

```python
N = 10000      
B = 15000     
Blocks = 2    
C = 10        
P = 4         
```

Feel free to scale these parameters up for larger experiments, keeping in mind available compute resources. Set BPR to something smaller if you run into memory issues.

---

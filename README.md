# 🧬🔒 PP-GWAS: Privacy Preserving Multi-Site Genome-wide Association Studies

This repository is associated with the paper titled **["PP-GWAS: Privacy Preserving Multi-Site Genome-wide Association Studies"](https://www.nature.com/articles/s41467-025-66771-z)**. It provides all the necessary code to reproduce the experiments described in the paper, including synthetic dataset generation, dataset loading, distributed computations, and experiment results.

---

## 📂 Repository Structure

- [**Code/**](Code/)  
  PP-GWAS implementation for HPC clusters (SLURM-ready).  
  Default outputs are written to [**test_site/**](test_site/).

- [**Data_Generation/**](Data_Generation/)  
  Scripts for simulating synthetic genotype/phenotype data with `pysnptools`.  
  Generated datasets are saved to [**test_site/**](test_site/) by default.

- [**REGENIE/**](REGENIE/)  
  Instructions to restructure synthetic data into REGENIE-supported formats and run REGENIE.

- [**Results/**](Results/)  
  `.txt` outputs from the experiments reported in the paper.

- [**local_run/**](local_run/)  
  Self-contained setup (Jupyter notebook + simple GUI) for running PP-GWAS locally on small datasets.

---

## 🧭 Pick a Workflow

- **Prepare data and run PP-GWAS locally (small data):** use [**local_run/**](local_run/).  
- **Prepare data and run PP-GWAS on a cluster:** use [**Code/**](Code/). Outputs go to [**test_site/**](test_site/).

- **Run REGENIE:** go to [**REGENIE/**](REGENIE/) and follow the steps there.  
- **Prepare data:** use [**Data_Generation/**](Data_Generation/); outputs appear in [**test_site/**](test_site/).

---

## 📝 Citation

Please consider citing our work if it is beneficial to your research. 

- **Paper**
```bibtex
@article{swaminathan2025pp,
  title={Pp-gwas: Privacy preserving multi-site genome-wide association studies},
  author={Swaminathan, Arjhun and Hannemann, Anika and {\"U}nal, Ali Burak and Pfeifer, Nico and Akg{\"u}n, Mete},
  journal={Nature Communications},
  year={2025},
  publisher={Nature Publishing Group UK London}
}
```
- **Code**
```bibtex
@software{swaminathan2025code,
  author   = {Swaminathan, Arjhun and Hannemann, Anika and {\"U}nal, Ali Burak and Pfeifer, Nico and Akg{\"u}n, Mete},
  title    = {PP-GWAS: Privacy Preserving Multi-Site Genome-wide Association Studies — code},
  version  = {v1.0},
  publisher= {Zenodo},
  year     = {2025},
  doi      = {10.5281/zenodo.17580283},
  url      = {https://doi.org/10.5281/zenodo.17580283}
}
```

---

## 📜 License

This project is released under the MIT License. See [LICENSE](LICENSE) for details.

---

## 📧 Contact for Questions
`arjhun.swaminathan@uni-tuebingen.de`, `mete.akguen@uni-tuebingen.de`

# PP-GWAS: Privacy Preserving Multi-Site Genome-wide Association Studies

This repository is associated with the paper titled **["PP-GWAS: Privacy Preserving Multi-Site Genome-wide Association Studies"](https://arxiv.org/abs/2410.08122)**. It provides all the necessary code to reproduce the experiments described in the paper, including synthetic dataset generation, dataset loading, distributed computations, and experiment results.

## Repository Structure

- **Datasets Folders**: 
  - **Code**: Contains all relevant code to implementing PP-GWAS on the cluster using SLURM. 
  - **Data_Generation**: This folder includes the code used for generating synthetic data using the pysnptools library. The generated data is stored in **test_site**.
  - **REGENIE**: Contains details on how one can restructure synthetic data to fit the formats supported by REGENIE. It also includes instructions on how to run REGENIE on the data.
  - **Results**: This folder holds all output .txt files generated from the experiments reported in the paper.

  - **local_run**: A fully self-contained folder with a jupyter notebook and GUI suitable for local testing of PP-GWAS. 

## Citation

Please consider citing our work if it is beneficial to your research. 

```bibtex
@article{swaminathan2024pp,
  title={PP-GWAS: Privacy Preserving Multi-Site Genome-wide Association Studies},
  author={Swaminathan, Arjhun and Hannemann, Anika and {\"U}nal, Ali Burak and Pfeifer, Nico and Akg{\"u}n, Mete},
  journal={arXiv preprint arXiv:2410.08122},
  year={2024}
}
```

## License

This project is released under the MIT License. See [LICENSE](LICENSE) for details.

## Contact for Questions
`arjhun.swaminathan@uni-tuebingen.de`, `mete.akguen@uni-tuebingen.de`

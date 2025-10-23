This folder contains scripts and commands to prepare data for and run REGENIE.

## Files

- **`plink_file_generation.py`**  
  A Python script that generates essential input files for PLINK (`.map`, `.ped`). These include a phenotypes.txt file, a covairates.txt file, a .map file and a .ped file.
  
- **`plink_bash_commands.sh`**  
  A bash script that uses PLINK commands to create `.bed` and `.bgen` files from `.map` and `.ped` files. 

- **`regenie_bash_commands.sh`**  
  A bash script that runs quantitative analysis in REGENIE using the PLINK-generated files. 

---

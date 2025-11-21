# miRNA_Picea_abies
Analysis of miRNA dynamics in Norway spruce (Picea abies) embryogenesis:


This repository contains scripts and data for analyzing microRNA (miRNA) dynamics during Norway spruce embryogenesis, including differential gene expression (DGE), variance-stabilized transformation (VST), principal component analysis (PCA), and other analyses for somatic embryogenesis (SE) and zygotic embryogenesis (ZE).


## Prerequisites
R: For running .R scripts (e.g., DESeq2, GO_enrich2.R). Install required packages such as DESeq2 and others specified in the scripts.

Python: For running .py scripts. Ensure dependencies like pandas, matplotlib, and seaborn are installed.

Input data files must be placed in the input/ directory before running scripts.

## Usage

1. Clone the repository:
```
git clone https://github.com/ahsan-rizvi/miRNA_Picea_abies.git
```
2. Place input data files in the input/ directory.

3. Run the scripts:
```
bash pip_spruce.sh
```
4. Output files will be generated in the output/ directory.
   
## Notes
Ensure input files are formatted as expected by each script (refer to script documentation or comments).

Scripts assume specific data structures; verify compatibility with your data before execution.

## Contact

Ahsan Z Rizvi
ahsan.zaigam@gmail.com


# miRNA_Picea_abies
Analysis of miRNA dynamics in Norway spruce (Picea abies) embryogenesis:
This repository contains scripts and data for analyzing microRNA (miRNA) dynamics during Norway spruce embryogenesis, including differential gene expression (DGE), variance-stabilized transformation (VST), principal component analysis (PCA), and other analyses for somatic embryogenesis (SE) and zygotic embryogenesis (ZE).
## File and folder Structure
File and folder structure of the repository is as follows:
```.
├── input/  # This folder contains all input files required for scripts
├── output/ # This folder contains all output files generated from the exicutions of the scripts
|── script/ # This folder contains scripts developed for miRNA analysis
|   ├── DESeq_SE.R   # R script for DGE, VST, PCA analysis for SE
|   ├── DESeq_ZE.R   # R script for DGE, VST, PCA analysis for ZE
|   ├── GO_enrich2.R # R script for GO enrichment analysis
|   ├── Heatmap_group_SE_M3.py # Python script for Heatmap generation for miRNA groups
|   ├── RPM_calculator.py # Script for counts to Reads Per Million (RPM) calculator
|   ├── SE_UP_DN_analysis.py  # This script will check miRNA mRNA opposite expressions of SE
|   ├── ZE_UP_DN_analysis.py  # This script will check miRNA mRNA opposite expressions of ZE
|   ├── VST_correlation_check_SE.py # This script will correlate VST based expressions of miRNAs and mRNAs
|   ├── ZE_UP_DN_analysis.py  # This script will check miRNA mRNA opposite expressions of ZE
|   ├── fasta_seq_size_checker.py # This script will check the distributions of miRNAs sequence size
|   ├── mirDP2_mature_miRNA_analysis.py # mirDP2 output Picea_miRNA analysis
|   ├── mirDP2_result_analysis.py # mirDP2 result comparision
|   ├── pre_processing_for_mirDP2.py # Input file maker for mirDP2 
|── secondry_files/ # This folder contains secondry files generated from the exicutions of the scripts
├── README.md       # Project overview and instructions
├── pip_spruce.sh   # pipeline for miRNA analysis         
```
## Prerequisites
R: For running .R scripts (e.g., DESeq2, GO_enrich2.R). Install required packages such as DESeq2 and others specified in the scripts.

Python: For running .py scripts. Ensure dependencies like pandas, matplotlib, and seaborn are installed (use pip install -r requirements.txt if provided).

Input data files must be placed in the input/ directory before running scripts.


## Notes
Ensure input files are formatted as expected by each script (refer to script documentation or comments).

Scripts assume specific data structures; verify compatibility with your data before execution.

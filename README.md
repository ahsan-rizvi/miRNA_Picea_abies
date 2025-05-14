miRNA_Picea_abies
Analysis of miRNA dynamics in Norway spruce (Picea abies) embryogenesis
This repository contains scripts and data for analyzing microRNA (miRNA) dynamics during Norway spruce embryogenesis, including differential gene expression (DGE), variance-stabilized transformation (VST), principal component analysis (PCA), and other analyses for somatic embryogenesis (SE) and zygotic embryogenesis (ZE).
File and Folder Structure
The repository is organized as follows:
├── input/                        # Contains all input files required for scripts
├── output/                       # Contains all output files generated from script executions
├── script/                       # Contains scripts developed for miRNA analysis
│   ├── DESeq_SE.R                # R script for DGE, VST, and PCA analysis for SE
│   ├── DESeq_ZE.R                # R script for DGE, VST, and PCA analysis for ZE
│   ├── GO_enrich2.R              # R script for GO enrichment analysis
│   ├── Heatmap_group_SE_M3.py    # Python script for heatmap generation for miRNA groups
│   ├── RPM_calculator.py         # Python script for converting counts to Reads Per Million (RPM)
│   ├── SE_UP_DN_analysis.py      # Python script to check miRNA-mRNA opposite expressions for SE
│   ├── ZE_UP_DN_analysis.py      # Python script to check miRNA-mRNA opposite expressions for ZE
│   ├── VST_correlation_check_SE.py # Python script to correlate VST-based expressions of miRNAs and mRNAs for SE
│   ├── fasta_seq_size_checker.py # Python script to check distributions of miRNA sequence sizes
│   ├── mirDP2_mature_miRNA_analysis.py # Python script for analyzing mirDP2 output for Picea abies miRNAs

Prerequisites

R: For running .R scripts (e.g., DESeq2, GO_enrich2.R). Install required packages such as DESeq2 and others specified in the scripts.
Python: For running .py scripts. Ensure dependencies like pandas, matplotlib, and seaborn are installed (use pip install -r requirements.txt if provided).
Input data files must be placed in the input/ directory before running scripts.

Usage

Clone the repository:git clone https://github.com/<your-username>/miRNA_Picea_abies.git


Place input data files in the input/ directory.
Run the scripts in the script/ directory as needed. For example:Rscript script/DESeq_SE.R
python script/RPM_calculator.py


Output files will be generated in the output/ directory.

Notes

Ensure input files are formatted as expected by each script (refer to script documentation or comments).
Scripts assume specific data structures; verify compatibility with your data before execution.


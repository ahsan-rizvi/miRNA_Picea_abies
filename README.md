# miRNA_Picea_abies
miRNA dynamics in Norway spruce embryogenesis
## File Structure
The file structure of the repository is as follows:
```.
├── input/  # This folder contains all input files required for scripts
├── output/ # This folder contains all output files generated from the exicutions of the scripts
|── script/ # This folder contains scripts developed for miRNA analysis
|   ├── DESeq_SE.R   # R script for DGE, VST, PCA analysis for SE
|   ├── DESeq_ZE.R   # R script for DGE, VST, PCA analysis for ZE
|   ├── GO_enrich2.R # R script for GO enrichment analysis
|   ├── Heatmap_group_SE_M3.py # Python script for Heatmap generation for miRNA groups
|   ├── RPM_calculator.py # Script for counts to Reads Per Million (RPM) calculator
|   ├── SE_UP_DN_analysis.py  # This script will check miRNA mRNA opposite expressions
|── secondry_files/ # This folder contains secondry files generated from the exicutions of the scripts
├── README.md       # Project overview and instructions
├── pip_spruce.sh   # pipeline for miRNA analysis         
```

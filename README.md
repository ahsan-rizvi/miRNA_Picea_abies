# miRNA_Picea_abies
miRNA dynamics in Norway spruce embryogenesis
## File Structure
The file structure of the repository is as follows:
```.
├── data/
│   ├── lncRNA_feature.pkl   # Pickle file containing all lncRNA features
│   ├── miRNA_feature.pkl    # Pickle file containing all miRNA features
│   ├── lncRNA_idx.csv       # CSV file with all lncRNA names
│   ├── miRNA_idx.csv        # CSV file with all miRNA names
│   └── splits.pkl           # Pickle file containing train/test split data
├── code/
│   ├── model.py             # SGAT-TM model architecture
│   ├── dataset.py           # Data loading and preprocessing
│   ├── funcs.py             # Metric and utility functions
│   ├── main.py              # Training and evaluation script
│   └── layer.py             # Custom layers for the model
├── README.md                # Project overview and instructions
└── requirements.txt         # Python dependencies
```

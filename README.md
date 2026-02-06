# Physics-Informed pKa Predictor

## Introduction
A Physics-Informed Machine Learning project to predict pKa values from molecular SMILES.
Combines semi-empirical quantum chemistry (xtb) with chemoinformatics (RDKit).

## Setup
1. Install dependencies:
   ```bash
   conda env create -f environment.yml
   conda activate pka_predictor
   ```

## Usage
Run the main pipeline:
```bash
python main.py --data data/raw/data.csv
```

## Structure
- `data/`: Data storage (Raw, Interim, Processed)
- `src/`: Source code
- `models/`: Trained models
- `temp/`: Temporary files for xtb calculations

import numpy as np
import pandas as pd
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import xgboost as xgb
from sklearn.model_selection import KFold, cross_validate
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
from sklearn.preprocessing import StandardScaler
import joblib
import logging
from typing import List, Dict, Any, Tuple

logger = logging.getLogger(__name__)

class ModelTrainer:
    """
    Handles feature generation, model training, and evaluation.
    """
    
    def __init__(self, fp_radius: int = 2, fp_bits: int = 2048):
        self.fp_radius = fp_radius
        self.fp_bits = fp_bits
        # Updated feature keys including gsolv
        self.feature_keys = [
            'total_energy', 'homo', 'lumo', 'fermi_level', 
            'gap', 'max_pos_charge', 'min_neg_charge', 'max_h_charge', 'gsolv'
        ]
        self.model = xgb.XGBRegressor(
            n_estimators=1000,
            learning_rate=0.05,
            max_depth=6,
            n_jobs=-1,
            random_state=42
        )
        self.scaler = StandardScaler()

    def generate_fingerprints(self, smiles_list: List[str]) -> np.ndarray:
        """
        Generate Morgan fingerprints for a list of SMILES using optimized numpy conversion.
        """
        fps = []
        valid_indices = []
        
        for i, smi in enumerate(smiles_list):
            mol = Chem.MolFromSmiles(smi)
            if mol:
                # Optimized conversion to numpy array
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, self.fp_radius, nBits=self.fp_bits)
                arr = np.zeros((0,), dtype=np.int8)
                DataStructs.ConvertToNumpyArray(fp, arr)
                fps.append(arr)
                valid_indices.append(i)
            else:
                fps.append(np.zeros((self.fp_bits,), dtype=np.int8))
                
        return np.vstack(fps)

    def prepare_data(self, data: List[Dict[str, Any]]) -> Tuple[np.ndarray, np.ndarray]:
        """
        Prepare X and y matrices from a list of data dictionaries.
        Each dict should contain: 'smiles', 'pka', and physical feature keys.
        """
        smiles = []
        y = []
        phys_feats = []
        
        for entry in data:
            # Check if essential keys are present (physics keys are checked below)
            if 'smiles' not in entry or 'pka' not in entry:
                continue
                
            # Extract physical features safely
            # If a feature is missing, we might want to skip or impute.
            # For now, let's skip to ensure data quality.
            if not all(k in entry for k in self.feature_keys):
                # Check which keys are missing for debugging if needed
                # missing = [k for k in self.feature_keys if k not in entry]
                continue
            
            smiles.append(entry['smiles'])
            y.append(entry['pka'])
            
            feat_vec = [float(entry[k]) for k in self.feature_keys]
            phys_feats.append(feat_vec)
            
        if not smiles:
            raise ValueError("No valid data points found for training.")
            
        # Generate fingerprints
        X_fp = self.generate_fingerprints(smiles)
        
        # Scale physical features
        X_phys = np.array(phys_feats)
        X_phys_scaled = self.scaler.fit_transform(X_phys)
        
        # Concatenate: [Fingerprints, Scaled Physical Features]
        X = np.hstack([X_fp, X_phys_scaled])
        y = np.array(y)
        
        return X, y

    def evaluate(self, X: np.ndarray, y: np.ndarray, n_splits: int = 5):
        """
        Perform K-Fold Cross Validation and print metrics.
        """
        kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
        
        scoring = {
            'r2': 'r2',
            'neg_mae': 'neg_mean_absolute_error',
            'neg_rmse': 'neg_root_mean_squared_error'
        }
        
        results = cross_validate(self.model, X, y, cv=kf, scoring=scoring)
        
        r2 = results['test_r2'].mean()
        mae = -results['test_neg_mae'].mean()
        rmse = -results['test_neg_rmse'].mean()
        
        logger.info(f"CV Results ({n_splits}-fold):")
        logger.info(f"R2: {r2:.4f}")
        logger.info(f"MAE: {mae:.4f}")
        logger.info(f"RMSE: {rmse:.4f}")
        
        return {'r2': r2, 'mae': mae, 'rmse': rmse}

    def analyze_feature_importance(self) -> pd.DataFrame:
        """
        Analyze and return feature importance.
        """
        if not hasattr(self.model, 'feature_importances_'):
            logger.warning("Model has not been trained yet.")
            return pd.DataFrame()
            
        # Generate feature names
        fp_names = [f"FP_{i}" for i in range(self.fp_bits)]
        phys_names = self.feature_keys
        all_names = fp_names + phys_names
        
        importances = self.model.feature_importances_
        
        df_imp = pd.DataFrame({
            'Feature': all_names,
            'Importance': importances
        }).sort_values(by='Importance', ascending=False)
        
        # Log top 10 features
        logger.info("Top 10 Important Features:")
        logger.info(df_imp.head(10).to_string(index=False))
        
        return df_imp

    def train_final_model(self, X: np.ndarray, y: np.ndarray, save_path: str = "model_pkg.joblib"):
        """
        Train on full dataset and save complete model package.
        """
        self.model.fit(X, y)
        
        # Analyze importance
        self.analyze_feature_importance()
        
        # Create model package
        model_pkg = {
            'model': self.model,
            'scaler': self.scaler,
            'fp_radius': self.fp_radius,
            'fp_bits': self.fp_bits,
            'feature_keys': self.feature_keys
        }
        
        joblib.dump(model_pkg, save_path)
        logger.info(f"Full model package saved to {save_path}")

    def load_model(self, path: str):
        """
        Load a trained model package.
        """
        if not os.path.exists(path):
            raise FileNotFoundError(f"Model file not found at {path}")
            
        model_pkg = joblib.load(path)
        self.model = model_pkg['model']
        self.scaler = model_pkg['scaler']
        self.fp_radius = model_pkg.get('fp_radius', self.fp_radius)
        self.fp_bits = model_pkg.get('fp_bits', self.fp_bits)
        self.feature_keys = model_pkg.get('feature_keys', self.feature_keys)
        
        logger.info(f"Model loaded from {path}")

    def predict(self, data: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Generate predictions for new data.
        Returns a list of dicts with 'id', 'smiles', 'predicted_pka'.
        """
        smiles = []
        ids = []
        phys_feats = []
        valid_entries = []
        
        for entry in data:
            if 'smiles' not in entry:
                continue
            
            # Check features
            if not all(k in entry for k in self.feature_keys):
                continue
                
            smiles.append(entry['smiles'])
            ids.append(entry.get('id', 'unknown'))
            
            feat_vec = [float(entry[k]) for k in self.feature_keys]
            phys_feats.append(feat_vec)
            valid_entries.append(entry)
            
        if not smiles:
            return []
            
        # Generate fingerprints
        X_fp = self.generate_fingerprints(smiles)
        
        # Scale physical features (transform only)
        X_phys = np.array(phys_feats)
        X_phys_scaled = self.scaler.transform(X_phys)
        
        # Concatenate
        X = np.hstack([X_fp, X_phys_scaled])
        
        # Predict
        preds = self.model.predict(X)
        
        results = []
        for i, pred in enumerate(preds):
            res = valid_entries[i].copy()
            res['predicted_pka'] = float(pred)
            results.append(res)
            
        return results

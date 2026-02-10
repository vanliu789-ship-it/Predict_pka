import os
import argparse
import pandas as pd
import logging
from multiprocessing import Pool, cpu_count
from functools import partial
from tqdm import tqdm
from typing import List, Dict, Any, cast

from src.preprocessor import MoleculePreprocessor
from src.physics_engine import PhysicsEngine
from src.model_trainer import ModelTrainer

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("project.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Global variables for worker processes (Lazy Loading)
_preprocessor = None
_physics_engine = None

def init_worker(output_dir, xtb_work_dir):
    """
    Initialize global objects for each worker process.
    """
    global _preprocessor, _physics_engine
    _preprocessor = MoleculePreprocessor(output_dir=output_dir)
    _physics_engine = PhysicsEngine(work_dir=xtb_work_dir)

def process_single_molecule(row):
    """
    Worker function for parallel processing.
    Uses globally initialized objects to avoid re-instantiation overhead.
    
    Args:
        row: Dictionary containing (smiles, id, initial_charge, uhf).
        
    Returns:
        Dictionary with physical features or None if failed.
    """
    try:
        global _preprocessor, _physics_engine
        
        # Fallback if init_worker wasn't called (e.g. single threaded debug)
        if _preprocessor is None:
            _preprocessor = MoleculePreprocessor()
        if _physics_engine is None:
            _physics_engine = PhysicsEngine()
            
        smiles = row['smiles']
        mol_id = str(row['id'])
        
        # Dynamic Charge & UHF
        # Default to 0 if not provided
        charge = int(row.get('initial_charge', 0))
        uhf = int(row.get('uhf', 0))
        
        # 1. Preprocess (Generate 3D)
        xyz_path = _preprocessor.process_smiles(smiles, mol_id)
        if not xyz_path:
            return None
        
        # 2. Physics Calculation
        props = _physics_engine.run_calculation(xyz_path, mol_id, charge=charge, uhf=uhf)
        
        if not props:
            return None
            
        # Merge results
        result = row.copy()
        result.update(props)
        return result
        
    except Exception as e:
        logger.error(f"Worker failed for {row.get('id', 'unknown')}: {e}")
        return None

def main():
    parser = argparse.ArgumentParser(description="Physics-Informed pKa Predictor")
    parser.add_argument("--data", type=str, default="data/raw/data.csv", help="Path to input CSV")
    parser.add_argument("--n_jobs", type=int, default=max(1, cpu_count() - 1), help="Number of parallel jobs")
    parser.add_argument("--output_dir", type=str, default="data/interim/structures", help="Directory for 3D structures")
    parser.add_argument("--xtb_dir", type=str, default="temp", help="Directory for xtb calculations")
    parser.add_argument("--processed_file", type=str, default="data/processed/features.csv", help="Path to save processed features")
    parser.add_argument("--model_path", type=str, default="models/pka_model_v1.joblib", help="Path to save trained model")
    parser.add_argument("--chunk_size", type=int, default=50, help="Number of molecules to save at once (checkpointing)")
    
    args = parser.parse_args()
    
    # 1. Load Data
    if not os.path.exists(args.data):
        logger.error(f"Data file {args.data} not found.")
        # Create dummy data if not exists (for first run)
        if not os.path.exists(os.path.dirname(args.data)):
             os.makedirs(os.path.dirname(args.data), exist_ok=True)
             
        logger.info("Creating a dummy dataset for demonstration at data/raw/data.csv...")
        df = pd.DataFrame({
            'id': ['mol1', 'mol2', 'mol3'],
            'smiles': ['CC(=O)O', 'c1ccccc1O', 'C1CCCCC1C(=O)O'], 
            'pka': [4.76, 9.95, 4.90],
            'initial_charge': [0, 0, 0],
            'uhf': [0, 0, 0]
        })
        df.to_csv(args.data, index=False)
    else:
        df = pd.read_csv(args.data)
    
    # Ensure ID column exists
    if 'id' not in df.columns:
        df['id'] = [f"mol_{i}" for i in range(len(df))]
        
    # Determine Mode
    mode = 'train'
    if 'pka' not in df.columns:
        mode = 'predict'
        logger.info("No 'pka' column found. Switching to PREDICTION mode.")
    else:
        logger.info("Found 'pka' column. Running in TRAINING mode.")

    required_cols = ['smiles']
    if mode == 'train':
        required_cols.append('pka')
        
    if not all(col in df.columns for col in required_cols):
        logger.error(f"Input CSV must contain columns: {required_cols}")
        return

    # Checkpointing Logic: Resume from partial file
    processed_ids = set()
    partial_file = args.processed_file
    df_existing = None
    
    if os.path.exists(partial_file):
        try:
            df_existing = pd.read_csv(partial_file)
            if 'id' in df_existing.columns:
                processed_ids = set(df_existing['id'].astype(str))
            logger.info(f"Found existing processed file with {len(processed_ids)} records. Resuming...")
        except Exception as e:
            logger.warning(f"Failed to read existing processed file: {e}. Starting fresh.")
    
    # Filter out already processed molecules
    # Ensure ID type consistency
    df['id'] = df['id'].astype(str)
    df_to_process = df[~df['id'].isin(processed_ids)]
    
    if len(df_to_process) == 0:
        logger.info("All molecules already processed!")
        if df_existing is not None:
            processed_data = cast(List[Dict[str, Any]], df_existing.to_dict('records')) # Load all for training
        else:
            logger.error("Unexpected state: All molecules processed but no existing data found.")
            return
    else:
        data_records = df_to_process.to_dict('records')
        logger.info(f"Processing {len(data_records)} new molecules (Total: {len(df)}).")

        # 2. Parallel Processing
        logger.info(f"Starting parallel processing with {args.n_jobs} jobs...")
        
        # Ensure output directories exist
        os.makedirs(os.path.dirname(args.processed_file), exist_ok=True)
        
        # We use imap_unordered for efficiency, but we need to handle writing carefully.
        # Simple approach: Collect chunks and append to CSV.
        
        with Pool(processes=args.n_jobs, initializer=init_worker, initargs=(args.output_dir, args.xtb_dir)) as pool:
            
            # Using imap (ordered) or imap_unordered. Ordered is nice for tracking but unordered is faster.
            # Let's use imap to keep it simple with tqdm, but we write incrementally.
            
            results_generator = pool.imap(process_single_molecule, data_records)
            
            batch_results = []
            total_new_processed = 0
            
            # Iterate through results with progress bar
            for result in tqdm(results_generator, total=len(data_records)):
                if result is not None:
                    batch_results.append(result)
                
                # Checkpoint saving
                if len(batch_results) >= args.chunk_size:
                    df_batch = pd.DataFrame(batch_results)
                    
                    # Append to CSV with column alignment
                    if os.path.exists(partial_file):
                        try:
                            # Read existing header to ensure column alignment
                            existing_cols = pd.read_csv(partial_file, nrows=0).columns.tolist()
                            # Reindex to match existing columns (adds NaNs if missing, drops extras if not in file)
                            df_batch = df_batch.reindex(columns=existing_cols)
                            header = False
                        except pd.errors.EmptyDataError:
                            # File exists but empty? Treat as new.
                            header = True
                        except Exception as e:
                            logger.warning(f"Error reading existing file header: {e}. Proceeding with potential risk.")
                            header = False
                    else:
                        header = True

                    df_batch.to_csv(partial_file, mode='a', header=header, index=False)
                    
                    total_new_processed += len(batch_results)
                    batch_results = [] # Clear buffer
            
            # Save remaining results
            if batch_results:
                df_batch = pd.DataFrame(batch_results)
                
                if os.path.exists(partial_file):
                    try:
                        existing_cols = pd.read_csv(partial_file, nrows=0).columns.tolist()
                        df_batch = df_batch.reindex(columns=existing_cols)
                        header = False
                    except pd.errors.EmptyDataError:
                        header = True
                    except Exception as e:
                        logger.warning(f"Error reading existing file header: {e}")
                        header = False
                else:
                    header = True
                    
                df_batch.to_csv(partial_file, mode='a', header=header, index=False)
                total_new_processed += len(batch_results)
                
        logger.info(f"Successfully processed {total_new_processed} new molecules.")
        
        # Reload full dataset for training
        if os.path.exists(partial_file):
            processed_data = pd.read_csv(partial_file).to_dict('records')
        else:
            logger.error("No processed data found. Calculation failed for all molecules.")
            return

    if not processed_data:
        logger.error("No valid data available for training. Exiting.")
        return

    # 3. Model Training or Prediction
    trainer = ModelTrainer()

    if mode == 'train':
        logger.info("Starting model training...")
        
        try:
            X, y = trainer.prepare_data(processed_data)
            
            # Evaluate
            trainer.evaluate(X, y)
            
            # Final Train and Save
            os.makedirs(os.path.dirname(args.model_path), exist_ok=True)
            trainer.train_final_model(X, y, save_path=args.model_path)
            
        except Exception as e:
            logger.error(f"Model training failed: {e}")
            
    else:
        logger.info("Starting prediction...")
        try:
            if not os.path.exists(args.model_path):
                logger.error(f"Model file not found at {args.model_path}. Cannot predict without a trained model.")
                return

            trainer.load_model(args.model_path)
            predictions = trainer.predict(processed_data)
            
            if predictions:
                pred_df = pd.DataFrame(predictions)
                
                # Construct output path
                base, ext = os.path.splitext(args.data)
                output_path = f"{base}_predictions{ext}"
                
                pred_df.to_csv(output_path, index=False)
                logger.info(f"Predictions saved to {output_path}")
            else:
                logger.warning("No predictions generated.")
                
        except Exception as e:
            logger.error(f"Prediction failed: {e}")

if __name__ == "__main__":
    # Windows support for multiprocessing
    import multiprocessing
    multiprocessing.freeze_support()
    main()

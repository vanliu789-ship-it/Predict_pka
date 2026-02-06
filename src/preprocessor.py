import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem.MolStandardize import rdMolStandardize
from typing import List, Optional, Tuple
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Top-level function for picklability on Windows
def _process_single_wrapper(args):
    """
    Wrapper to call instance method, but we can't easily pass instance method.
    So we reconstruct preprocessor or pass necessary data.
    Actually, passing the instance method `self.process_smiles` works if `self` is picklable.
    RDKit objects (SaltRemover, Uncharger) inside `self` might NOT be picklable.
    
    Safe approach: Create preprocessor inside worker or pass only necessary configs.
    """
    smiles, mol_id, output_dir = args
    
    # Re-instantiate locally to avoid pickling RDKit C++ objects if they are problematic
    # But this is expensive. Let's try to see if we can just do the work here.
    
    try:
        # 1. Convert SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        # 2. Standardize
        # Re-create standardizers locally
        remover = SaltRemover()
        # NOTE: StripMol might return a clean mol, or it might modify in place?
        # Documentation says it returns the stripped molecule.
        
        # Check if stripping removed everything
        stripped = remover.StripMol(mol)
        if stripped.GetNumAtoms() > 0:
            mol = stripped
        else:
            # If stripped is empty, it means the whole thing was considered a salt?
            # Or maybe SaltRemover defaults are too aggressive?
            # For now, if empty, keep original (maybe it wasn't a salt, just matched a pattern)
            # Actually, if it's empty, we probably shouldn't use it.
            # But let's log and keep original if stripped is empty but original wasn't.
            pass
            
        uncharger = rdMolStandardize.Uncharger()
        
        # uncharge might return None if it fails? No, usually returns mol.
        # But let's be safe.
        uncharged = uncharger.uncharge(mol)
        if uncharged is not None and uncharged.GetNumAtoms() > 0:
             mol = uncharged

        # 3. 3D & Optimize
        mol = Chem.AddHs(mol)
        
        # Check atoms again
        if mol.GetNumAtoms() == 0:
             logger.error(f"Molecule {mol_id} has no atoms after standardization")
             return None

        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        res = AllChem.EmbedMolecule(mol, params)
        
        if res == -1:
            params.useRandomCoords = True
            res = AllChem.EmbedMolecule(mol, params)
            
        if res == -1:
            return None
            
        try:
            AllChem.MMFFOptimizeMolecule(mol)
        except:
            pass

        # 4. Save
        xyz_block = Chem.MolToXYZBlock(mol)
        file_path = os.path.join(output_dir, f"{mol_id}.xyz")
        with open(file_path, "w") as f:
            f.write(xyz_block)
            
        return file_path
        
    except Exception as e:
        logger.error(f"Error in worker for {mol_id}: {e}")
        return None

class MoleculePreprocessor:
    """
    Handles conversion of SMILES to 3D structures, conformation generation,
    and preliminary energy optimization.
    """
    
    def __init__(self, output_dir: str = "data/structures"):
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)
        # We don't initialize RDKit objects here for the pool if we use the static wrapper approach

    def process_smiles(self, smiles: str, mol_id: str) -> Optional[str]:
        # Keep this for single-threaded usage or debugging
        # But for batch, we use the wrapper
        return _process_single_wrapper((smiles, mol_id, self.output_dir))

    def process_batch(self, smiles_list: List[str], ids_list: List[str], n_jobs: int = 1) -> List[Tuple[str, Optional[str]]]:
        results = []
        
        tasks = [(smiles, mol_id, self.output_dir) for smiles, mol_id in zip(smiles_list, ids_list)]
        
        if n_jobs > 1:
            with ProcessPoolExecutor(max_workers=n_jobs) as executor:
                # Map returns results in order
                paths = list(executor.map(_process_single_wrapper, tasks))
                
                for mol_id, path in zip(ids_list, paths):
                    results.append((mol_id, path))
        else:
            # Serial execution
            for task in tasks:
                path = _process_single_wrapper(task)
                mol_id = task[1]
                results.append((mol_id, path))
                    
        return results

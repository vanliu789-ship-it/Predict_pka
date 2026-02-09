import os
import subprocess
import logging
from typing import Dict, Optional
import shutil
from src.utils import parse_xtb_output

logger = logging.getLogger(__name__)

class PhysicsEngine:
    """
    Wrapper for xtb (semi-empirical quantum chemistry) to extract physical properties.
    """
    
    def __init__(self, xtb_path: str = "xtb", work_dir: str = "temp"):
        """
        Initialize the physics engine.
        
        Args:
            xtb_path: Path to the xtb executable (default: assume in PATH).
            work_dir: Directory for temporary xtb files.
        """
        self.xtb_path = xtb_path
        
        # Auto-detect xtb in Conda environment on Windows if not found in PATH
        if shutil.which(self.xtb_path) is None and os.name == 'nt':
            # Try standard Conda Library/bin path
            conda_prefix = os.environ.get('CONDA_PREFIX')
            if conda_prefix:
                potential_path = os.path.join(conda_prefix, 'Library', 'bin', 'xtb.exe')
                if os.path.exists(potential_path):
                    self.xtb_path = potential_path
                    logger.info(f"Auto-detected xtb at: {self.xtb_path}")

        self.work_dir = work_dir
        os.makedirs(self.work_dir, exist_ok=True)
        
        # Verify xtb availability
        if shutil.which(self.xtb_path) is None and not os.path.exists(self.xtb_path):
            logger.warning(f"xtb executable not found at {self.xtb_path}. Ensure it is installed and in PATH.")

    def run_calculation(self, xyz_file: str, mol_id: str, charge: int = 0, uhf: int = 0) -> Optional[Dict[str, float]]:
        """
        Run xtb geometry optimization and property calculation.
        
        Args:
            xyz_file: Path to the input .xyz file.
            mol_id: Unique ID for the molecule (used for isolation).
            charge: Net charge of the molecule.
            uhf: Number of unpaired electrons (spin).
            
        Returns:
            Dictionary containing extracted properties, or None if failed.
        """
        # Create a specific directory for this molecule in temp to avoid conflicts
        mol_dir = os.path.join(self.work_dir, mol_id)
        os.makedirs(mol_dir, exist_ok=True)
        
        abs_xyz_path = os.path.abspath(xyz_file)
        
        # Command: xtb <file> --opt --alpb water --chrg <charge> --uhf <uhf>
        cmd = [
            self.xtb_path,
            abs_xyz_path,
            "--opt",
            "--alpb", "water",
            "--chrg", str(charge),
            "--uhf", str(uhf)
        ]
        
        # Set environment variables to limit threads per xtb process
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = "1"
        env["MKL_NUM_THREADS"] = "1"
        env["OMP_STACKSIZE"] = "2000M"  # Use 'M' suffix for compatibility
        
        try:
            # Run xtb
            result = subprocess.run(
                cmd,
                cwd=mol_dir,
                capture_output=True,
                text=True,
                encoding='utf-8',
                errors='replace',
                env=env,
                timeout=600  # 10 minutes timeout per molecule
            )
            
            if result.returncode != 0:
                logger.error(f"xtb failed for {mol_id}: {result.stderr}")
                # Don't return None immediately, try to parse what we have if useful?
                # No, if return code != 0, calculation crashed.
                return None
            
            # Read charges file if it exists
            charges_content = None
            charges_file = os.path.join(mol_dir, "charges")
            if os.path.exists(charges_file):
                try:
                    with open(charges_file, 'r') as f:
                        charges_content = f.read()
                except Exception as e:
                    logger.warning(f"Failed to read charges file for {mol_id}: {e}")

            # Parse output using utils
            output_text = result.stdout
            properties = parse_xtb_output(output_text, charges_content)
            
            # Check convergence flag
            if properties.get('converged', 0.0) == 0.0:
                 logger.warning(f"Optimization not converged for {mol_id}")
            
            # Clean up (remove molecule directory to save space)
            try:
                shutil.rmtree(mol_dir)
            except Exception as e:
                logger.warning(f"Failed to cleanup {mol_dir}: {e}")
            
            return properties
            
        except subprocess.TimeoutExpired:
            logger.error(f"xtb timed out for {mol_id}")
            return None
        except Exception as e:
            logger.error(f"Error running xtb for {mol_id}: {e}")
            return None

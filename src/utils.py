"""
Utility functions for pKa Predictor project.
"""
import os
import re
import logging
from typing import Dict, List, Any, Optional

logger = logging.getLogger(__name__)

def parse_xtb_output(output: str, charges_content: Optional[str] = None) -> Dict[str, float]:
    """
    Parse xtb standard output to extract physical properties.
    
    Args:
        output: Standard output content from xtb.
        charges_content: Content of the 'charges' file (optional).
    """
    props = {}
    
    # 0. Check Convergence
    if "GEOMETRY OPTIMIZATION CONVERGED" not in output:
        logger.warning("Geometry optimization did NOT converge!")
        props['converged'] = 0.0
    else:
        props['converged'] = 1.0

    # 1. Total Energy
    energy_match = re.search(r"\|\s+TOTAL ENERGY\s+(-?\d+\.\d+)\s+Eh", output)
    if energy_match:
        props['total_energy'] = float(energy_match.group(1))
        
    # 1.1 Solvation Free Energy (Gsolv)
    # Match various formats:
    # "| Gsolv   -0.123 Eh"
    # ":: -> Gsolv   -0.123 Eh"
    # "delta Gsolv   -0.123 Eh"
    gsolv_pattern = r"(?:\||::\s+->|delta)\s+Gsolv\s+(-?\d+\.\d+)\s+Eh"
    gsolv_match = re.search(gsolv_pattern, output, re.IGNORECASE)
    if gsolv_match:
        props['gsolv'] = float(gsolv_match.group(1))
        
    # 2. HOMO/LUMO/Fermi
    homo_match = re.search(r"\(HOMO\)\s+(-?\d+\.\d+)", output)
    if not homo_match:
            homo_match = re.search(r"\|\s+HOMO\s+(-?\d+\.\d+)\s+eV", output)
            
    lumo_match = re.search(r"\(LUMO\)\s+(-?\d+\.\d+)", output)
    if not lumo_match:
        lumo_match = re.search(r"\|\s+LUMO\s+(-?\d+\.\d+)\s+eV", output)
        
    fermi_match = re.search(r"\|\s+Fermi Level\s+(-?\d+\.\d+)\s+eV", output)
    
    if homo_match: props['homo'] = float(homo_match.group(1))
    if lumo_match: props['lumo'] = float(lumo_match.group(1))
    if fermi_match: props['fermi_level'] = float(fermi_match.group(1))
    
    if 'homo' in props and 'lumo' in props:
        props['gap'] = props['lumo'] - props['homo']

    # 3. Partial Charges (Mapping Logic)
    # Extract elements from stdout (reliable order)
    elements_from_stdout = extract_elements_from_output(output)
    
    charges = []
    
    # Try reading from charges file first
    if charges_content:
        file_charges_values = parse_charges_file_values(charges_content)
        
        # If lengths match, pair them up!
        if len(file_charges_values) == len(elements_from_stdout):
            for el, q in zip(elements_from_stdout, file_charges_values):
                charges.append({'element': el, 'charge': q})
        else:
            logger.warning(f"Charges file count ({len(file_charges_values)}) mismatch with stdout elements ({len(elements_from_stdout)}). Falling back to stdout.")
    
    # Fallback to stdout if file parsing failed or mismatched
    if not charges:
        charges = extract_charges_from_output(output)
        
    if charges:
        props['max_pos_charge'] = max(c['charge'] for c in charges)
        props['min_neg_charge'] = min(c['charge'] for c in charges)
        
        # Max Positive Charge on Hydrogen
        h_charges = [c['charge'] for c in charges if c['element'] == 'H']
        if h_charges:
            props['max_h_charge'] = max(h_charges)
        else:
            # Only set to 0.0 if we actually have charges but no H found (and not empty list)
             props['max_h_charge'] = 0.0
            
    return props

def parse_charges_file_values(content: str) -> List[float]:
    """
    Parse xtb 'charges' file and return list of charge values.
    Expected format: x y z q
    """
    charge_values = []
    lines = content.strip().split('\n')
    for line in lines:
        parts = line.strip().split()
        if len(parts) >= 4:
            try:
                q = float(parts[3]) # 4th column is charge
                charge_values.append(q)
            except ValueError:
                pass
    return charge_values

def extract_elements_from_output(output: str) -> List[str]:
    """
    Extract just the list of elements from the charge table in stdout.
    Used for mapping file charges to elements.
    """
    # Use the same logic as extract_charges but just return elements
    charges_dicts = extract_charges_from_output(output)
    return [c['element'] for c in charges_dicts]

def extract_charges_from_output(output: str) -> List[Dict[str, Any]]:
    """
    Extract atomic charges from xtb output.
    Prioritizes CM5 charges over Mulliken.
    """
    charges = []
    
    # 1. Try to find CM5 charges specifically
    # Pattern: "Mulliken/CM5 charges" OR just "CM5 charges"
    # We want the LAST occurrence of either.
    
    # Find all start indices of charge tables
    # Note: xtb might output "Mulliken charges" then "CM5 charges" later?
    # Or "Mulliken/CM5 charges" as a single header.
    # Let's look for "Mulliken/CM5 charges" first (standard xtb 6.x)
    
    pattern = r"(Mulliken/CM5 charges|Mulliken charges)"
    matches = list(re.finditer(pattern, output))
    
    if not matches:
        return []
    
    # Use the last match (converged geometry)
    last_match_idx = matches[-1].start()
    relevant_part = output[last_match_idx:]
    
    lines = relevant_part.split('\n')
    start_parsing = False
    
    for line in lines:
        if "qt" in line and "qA" in line: # Header
            continue
        if "--------" in line:
            if not start_parsing:
                start_parsing = True
                continue
            else:
                pass
        
        if start_parsing:
            parts = line.strip().split()
            # Standard line: index Symbol charge ...
            # 1 C -0.123
            if len(parts) < 3: 
                if len(charges) > 0: break
                continue
            
            try:
                # part[0] is index, part[1] is Element, part[2] is Charge
                element = parts[1]
                charge = float(parts[2])
                charges.append({'element': element, 'charge': charge})
            except ValueError:
                break
                
    return charges

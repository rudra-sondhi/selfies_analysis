"""
Metric calculation functions for SELFIES analysis.
"""

import re
import numpy as np
import selfies
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdFingerprintGenerator
from rdkit.DataStructs import TanimotoSimilarity
from typing import Optional


def selfies_to_smiles(selfies_str: str) -> Optional[str]:
    """Convert SELFIES to SMILES."""
    try:
        return selfies.decoder(selfies_str)
    except Exception:
        return None


def smiles_to_mol(smiles: Optional[str]) -> Optional[Chem.Mol]:
    """Convert SMILES to RDKit Mol object."""
    if smiles is None:
        return None
    try:
        return Chem.MolFromSmiles(smiles)
    except Exception:
        return None


def calculate_hdi(selfies_str: str) -> float:
    """
    Calculate Hydrogen Deficiency Index (HDI) for a SELFIES string.
    
    Parameters
    ----------
    selfies_str : str
        SELFIES representation of a molecule
    
    Returns
    -------
    float
        HDI value, or np.nan if calculation fails
    """
    try:
        smiles = selfies.decoder(selfies_str)
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return np.nan
        
        num_c = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        num_h = sum(atom.GetTotalNumHs() for atom in mol.GetAtoms())
        num_n = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
        num_x = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in [9, 17, 35, 53])
        
        return (2 * num_c + 2 - num_h + num_n - num_x) / 2
    except Exception:
        return np.nan
    
def calculate_sequence_accuracy(real_sf: str, pred_sf: str) -> float:
    """
    Calculate sequence-level accuracy between two SELFIES strings,
    comparing token sequences exactly.

    Parameters
    ----------
    real_sf : str
        Real/target SELFIES string
    pred_sf : str
        Predicted SELFIES string

    Returns
    -------
    float
        Sequence accuracy (1.0 if exact token match, 0.0 otherwise)
    """
    try:
        real_toks = re.findall(r'\[[^\]]+\]', real_sf)
        pred_toks = re.findall(r'\[[^\]]+\]', pred_sf)
        return float(real_toks == pred_toks)
    except Exception:
        return 0.0
    


def calculate_token_accuracy(real_sf: str, pred_sf: str) -> float:
    """
    Calculate token-level accuracy between two SELFIES strings.
    
    Parameters
    ----------
    real_sf : str
        Real/target SELFIES string
    pred_sf : str
        Predicted SELFIES string
    
    Returns
    -------
    float
        Token accuracy (0.0 to 1.0)
    """
    try:
        real_toks = re.findall(r'\[[^\]]+\]', real_sf)
        pred_toks = re.findall(r'\[[^\]]+\]', pred_sf)
        n = min(len(real_toks), len(pred_toks))
        if n == 0:
            return 0.0
        matches = sum(1 for i in range(n) if real_toks[i] == pred_toks[i])
        return matches / n
    except Exception:
        return 0.0


def calculate_tanimoto_similarity(mol1: Optional[Chem.Mol], mol2: Optional[Chem.Mol]) -> float:
    """
    Calculate Tanimoto similarity between two molecules.
    
    Parameters
    ----------
    mol1, mol2 : RDKit Mol objects or None
        Molecules to compare
    
    Returns
    -------
    float
        Tanimoto similarity (0.0 to 1.0)
    """
    if mol1 is None or mol2 is None:
        return 0.0
    
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fp1 = gen.GetCountFingerprint(mol1)
    fp2 = gen.GetCountFingerprint(mol2)
    return TanimotoSimilarity(fp1, fp2)


def calculate_molecular_weight(mol: Optional[Chem.Mol]) -> Optional[float]:
    """
    Calculate molecular weight.
    
    Parameters
    ----------
    mol : RDKit Mol object or None
    
    Returns
    -------
    float or None
        Molecular weight in Daltons
    """
    if mol is None:
        return None
    return Descriptors.MolWt(mol)
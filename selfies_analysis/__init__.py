"""
SELFIES Analysis Package

A package for analyzing SELFIES molecular representations, computing metrics,
and generating comparison plots.

New in this version:
- Individual metric computation methods (get_hdi, get_tanimoto_similarity, etc.)
- Efficient caching of computed values
- Support for predicted vs real molecule comparisons
"""

from .core import SELFIESAnalyzer
from .metrics import (
    calculate_hdi,
    calculate_token_accuracy,
    calculate_tanimoto_similarity,
    calculate_molecular_weight,
    selfies_to_smiles,
    smiles_to_mol
)
from .plots import (
    plot_token_accuracy_hist,
    plot_tanimoto_similarity_hist,
    plot_mw_comparison,
    plot_hdi_comparison
)

__version__ = "0.2.0"
__all__ = [
    "SELFIESAnalyzer",
    "calculate_hdi",
    "calculate_token_accuracy", 
    "calculate_tanimoto_similarity",
    "calculate_molecular_weight",
    "selfies_to_smiles",
    "smiles_to_mol",
    "plot_token_accuracy_hist",
    "plot_tanimoto_similarity_hist",
    "plot_mw_comparison",
    "plot_hdi_comparison"
]
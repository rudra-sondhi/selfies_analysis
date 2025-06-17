# selfies_analysis/core.py
"""
Core analyzer class for SELFIES analysis.
"""

import os
import pandas as pd
from typing import Union, Tuple, List, Dict, Optional
from rdkit import Chem

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
    plot_hdi_comparison,
    plot_molecule_grid
)


class SELFIESAnalyzer:
    """
    Analyzer for SELFIES molecular representations.
    
    Can analyze single SELFIES strings or pairs of real/predicted SELFIES.
    """
    
    def __init__(self, selfies_input: Union[str, Tuple[str, str], List[Tuple[str, str]], List[str], List[List[str]], pd.DataFrame]):
        """
        Initialize the analyzer.
        
        Parameters
        ----------
        selfies_input : str, tuple, list of tuples, list of strings, list of lists, or DataFrame
            - str: Single SELFIES string for basic analysis
            - tuple: (real_selfies, pred_selfies) pair
            - list of tuples: List of (real_selfies, pred_selfies) pairs
            - list of strings: List of SELFIES strings for basic analysis
            - list of lists: If all sublists have exactly 2 elements, treated as (real, pred) pairs.
                            Otherwise, flattened into a single list of SELFIES strings.
            - DataFrame: Must have 'real_selfies' and optionally 'pred_selfies' columns
        """
        self.df = self._prepare_dataframe(selfies_input)
        self._computed_metrics = {}
        
    def _prepare_dataframe(self, selfies_input) -> pd.DataFrame:
        """Convert various input types to a standardized DataFrame."""
        if isinstance(selfies_input, str):
            # Single SELFIES string
            return pd.DataFrame({'real_selfies': [selfies_input]})
            
        elif isinstance(selfies_input, tuple) and len(selfies_input) == 2:
            # Single pair
            return pd.DataFrame([selfies_input], columns=['real_selfies', 'pred_selfies'])
            
        elif isinstance(selfies_input, list):
            # Check if it's a list of tuples, strings, or lists
            if len(selfies_input) == 0:
                return pd.DataFrame(columns=['real_selfies'])
            
            # Check the first element to determine list type
            first_element = selfies_input[0]
            
            if isinstance(first_element, tuple):
                # List of (real, predicted) pairs
                return pd.DataFrame(selfies_input, columns=['real_selfies', 'pred_selfies'])
            elif isinstance(first_element, str):
                # List of SELFIES strings
                return pd.DataFrame({'real_selfies': selfies_input})
            elif isinstance(first_element, list):
                # List of lists - check if they represent pairs or just groups
                # If sublists have exactly 2 elements, treat as (real, pred) pairs
                # Otherwise, flatten into a single list of SELFIES
                if all(len(sublist) == 2 for sublist in selfies_input if isinstance(sublist, list)):
                    # All sublists have exactly 2 elements - treat as (real, pred) pairs
                    pairs = []
                    for sublist in selfies_input:
                        if isinstance(sublist, list) and len(sublist) == 2:
                            pairs.append((sublist[0], sublist[1]))
                        else:
                            raise ValueError("When using list of lists as pairs, all sublists must have exactly 2 elements")
                    return pd.DataFrame(pairs, columns=['real_selfies', 'pred_selfies'])
                else:
                    # Mixed lengths or not all pairs - flatten into single list
                    flattened = []
                    for sublist in selfies_input:
                        if isinstance(sublist, list):
                            flattened.extend(sublist)
                        else:
                            raise ValueError("When using list of lists, all elements must be lists")
                    return pd.DataFrame({'real_selfies': flattened})
            else:
                raise ValueError("List elements must be strings (SELFIES), tuples of (real, pred) SELFIES, or lists of SELFIES")
                
        elif isinstance(selfies_input, pd.DataFrame):
            # Already a DataFrame
            if 'real_selfies' not in selfies_input.columns:
                raise ValueError("DataFrame must have 'real_selfies' column")
            return selfies_input.copy()
            
        else:
            raise ValueError("Invalid input type. Expected str, tuple, list of tuples/strings/lists, or DataFrame")
    
    def _ensure_smiles_and_mols(self) -> None:
        """Ensure SMILES and Mol objects are computed for real molecules."""
        if 'real_smiles' not in self.df.columns:
            self.df['real_smiles'] = self.df['real_selfies'].apply(selfies_to_smiles)
        if 'real_mol' not in self.df.columns:
            self.df['real_mol'] = self.df['real_smiles'].apply(smiles_to_mol)
    
    def _ensure_pred_smiles_and_mols(self) -> None:
        """Ensure SMILES and Mol objects are computed for predicted molecules."""
        if 'pred_selfies' not in self.df.columns:
            raise ValueError("No predicted SELFIES available for comparison metrics")
        if 'pred_smiles' not in self.df.columns:
            self.df['pred_smiles'] = self.df['pred_selfies'].apply(selfies_to_smiles)
        if 'pred_mol' not in self.df.columns:
            self.df['pred_mol'] = self.df['pred_smiles'].apply(smiles_to_mol)
    
    def get_hdi(self) -> Union[List[float], float, Dict[str, Union[List[float], float]]]:
        """
        Get HDI values for the molecules.
        
        Returns
        -------
        Union[List[float], float, Dict[str, Union[List[float], float]]]
            - If only real SELFIES: HDI values for real molecules
            - If real and predicted SELFIES: Dict with 'real' and 'predicted' keys
            Returns single float if only one molecule, otherwise list.
        """
        # Always ensure real HDI is computed
        if 'real_hdi' not in self.df.columns:
            self.df['real_hdi'] = self.df['real_selfies'].apply(calculate_hdi)
        
        real_hdi = self.df['real_hdi'].dropna().tolist()
        
        # Check if we have predicted SELFIES
        if 'pred_selfies' in self.df.columns and self.df['pred_selfies'].notna().any():
            # Compute predicted HDI
            if 'pred_hdi' not in self.df.columns:
                self.df['pred_hdi'] = self.df['pred_selfies'].apply(
                    lambda x: calculate_hdi(x) if pd.notna(x) else None
                )
            
            pred_hdi = self.df['pred_hdi'].dropna().tolist()
            
            # Return both as a dictionary
            result = {
                'real': real_hdi[0] if len(real_hdi) == 1 else real_hdi,
                'predicted': pred_hdi[0] if len(pred_hdi) == 1 else pred_hdi
            }
            return result
        else:
            # Only real data available
            return real_hdi[0] if len(real_hdi) == 1 else real_hdi
    
    def get_molecular_weight(self) -> Union[List[float], float, Dict[str, Union[List[float], float]]]:
        """
        Get molecular weight values for the molecules.
        
        Returns
        -------
        Union[List[float], float, Dict[str, Union[List[float], float]]]
            - If only real SELFIES: Molecular weight values for real molecules
            - If real and predicted SELFIES: Dict with 'real' and 'predicted' keys
            Returns single float if only one molecule, otherwise list.
        """
        self._ensure_smiles_and_mols()
        
        # Always compute real molecular weights
        if 'real_mol_wt' not in self.df.columns:
            self.df['real_mol_wt'] = self.df['real_mol'].apply(
                lambda m: calculate_molecular_weight(m) if m is not None else None
            )
        
        real_weights = self.df['real_mol_wt'].dropna().tolist()
        
        # Check if we have predicted SELFIES
        if 'pred_selfies' in self.df.columns and self.df['pred_selfies'].notna().any():
            self._ensure_pred_smiles_and_mols()
            
            # Compute predicted molecular weights
            if 'pred_mol_wt' not in self.df.columns:
                self.df['pred_mol_wt'] = self.df['pred_mol'].apply(
                    lambda m: calculate_molecular_weight(m) if m is not None else None
                )
            
            pred_weights = self.df['pred_mol_wt'].dropna().tolist()
            
            # Return both as a dictionary
            result = {
                'real': real_weights[0] if len(real_weights) == 1 else real_weights,
                'predicted': pred_weights[0] if len(pred_weights) == 1 else pred_weights
            }
            return result
        else:
            # Only real data available
            return real_weights[0] if len(real_weights) == 1 else real_weights
    
    def get_token_accuracy(self) -> Union[List[float], float]:
        """
        Get token accuracy values comparing real vs predicted SELFIES.
        
        Returns
        -------
        List[float] or float
            Token accuracy values (0-1). Returns single float if only one pair, 
            otherwise list.
            
        Raises
        ------
        ValueError
            If no predicted SELFIES are available for comparison.
        """
        if 'pred_selfies' not in self.df.columns:
            raise ValueError("No predicted SELFIES available for token accuracy calculation")
        
        if 'token_accuracy' not in self.df.columns:
            self.df['token_accuracy'] = self.df.apply(
                lambda row: calculate_token_accuracy(row['real_selfies'], row['pred_selfies']),
                axis=1
            )
        
        accuracies = self.df['token_accuracy'].tolist()
        return accuracies[0] if len(accuracies) == 1 else accuracies
    
    def get_tanimoto_similarity(self) -> Union[List[float], float]:
        """
        Get Tanimoto similarity values comparing real vs predicted molecules.
        
        Returns
        -------
        List[float] or float
            Tanimoto similarity values (0-1). Returns single float if only one pair, 
            otherwise list.
            
        Raises
        ------
        ValueError
            If no predicted SELFIES are available for comparison.
        """
        self._ensure_smiles_and_mols()
        self._ensure_pred_smiles_and_mols()
        
        if 'tanimoto_similarity' not in self.df.columns:
            self.df['tanimoto_similarity'] = self.df.apply(
                lambda row: calculate_tanimoto_similarity(row['real_mol'], row['pred_mol']),
                axis=1
            )
        
        similarities = self.df['tanimoto_similarity'].tolist()
        return similarities[0] if len(similarities) == 1 else similarities
    
    def plot_molecule_grid(self, n_samples: int = 4, 
                          mol_size: Tuple[int, int] = (300, 300),
                          save_dir: str = 'plots',
                          filename: str = 'molecule_comparison_grid.png',
                          title: Optional[str] = None,
                          show_indices: bool = True,
                          show_tanimoto: bool = False,
                          rows: Optional[int] = None) -> None:
        """
        Plot a grid comparing real vs predicted molecules using RDKit molecular drawings.
        
        Parameters
        ----------
        n_samples : int, default=4
            Number of molecule pairs to display in the grid
        mol_size : Tuple[int, int], default=(300, 300)
            Size of each molecule image in pixels (width, height)
        save_dir : str, default='plots'
            Directory to save the plot
        filename : str, default='molecule_comparison_grid.png'
            Filename for the saved plot
        title : str, optional
            Custom title for the plot. If None, uses default title
        show_indices : bool, default=True
            Whether to show sample indices in the grid
        show_tanimoto : bool, default=False
            Whether to display Tanimoto similarity scores for each pair
        rows : int, optional
            Number of rows to arrange the molecule pairs. If None, uses automatic layout.
            
        Raises
        ------
        ValueError
            If no predicted SELFIES are available for comparison.
        """
        if 'pred_selfies' not in self.df.columns:
            raise ValueError("No predicted SELFIES available for molecule grid visualization")
        
        # Extract pairs from dataframe
        pairs = list(zip(self.df['real_selfies'], self.df['pred_selfies']))
        
        # Use the standalone plotting function
        plot_molecule_grid(
            selfies_pairs=pairs,
            n_samples=n_samples,
            mol_size=mol_size,
            save_dir=save_dir,
            filename=filename,
            title=title,
            show_indices=show_indices,
            show_tanimoto=show_tanimoto,
            rows=rows
        )
    
    def compute_all_metrics(self) -> Dict[str, any]:
        """
        Compute all available metrics.
        
        Returns
        -------
        dict
            Dictionary containing all computed metrics and statistics
        """
        metrics = {}
        
        # Add SMILES conversions
        self._ensure_smiles_and_mols()
        
        # HDI
        hdi_values = self.get_hdi()
        metrics['hdi_values'] = hdi_values if isinstance(hdi_values, list) else [hdi_values]
        metrics['mean_hdi'] = self.df['real_hdi'].mean()
        
        # Molecular weight
        mw_values = self.get_molecular_weight()
        metrics['mol_weights'] = mw_values if isinstance(mw_values, list) else [mw_values]
        metrics['mean_mol_weight'] = self.df['real_mol_wt'].mean()
        
        # If we have predictions, compute comparison metrics
        if 'pred_selfies' in self.df.columns:
            self._ensure_pred_smiles_and_mols()
            
            # Get predicted HDI and molecular weight
            pred_hdi_values = self.get_hdi()
            pred_mw_values = self.get_molecular_weight()
            
            # Token accuracy
            token_acc_values = self.get_token_accuracy()
            metrics['token_accuracy'] = token_acc_values if isinstance(token_acc_values, list) else [token_acc_values]
            metrics['mean_token_accuracy'] = self.df['token_accuracy'].mean()
            
            # Tanimoto similarity
            tanimoto_values = self.get_tanimoto_similarity()
            metrics['tanimoto_similarity'] = tanimoto_values if isinstance(tanimoto_values, list) else [tanimoto_values]
            metrics['mean_tanimoto_similarity'] = self.df['tanimoto_similarity'].mean()
        
        self._computed_metrics = metrics
        return metrics
    
    def plot_all(self, save_dir: str = 'plots', prefix: str = '', 
                include_molecule_grid: bool = False, n_grid_samples: int = 4,
                show_tanimoto_on_grid: bool = False) -> None:
        """
        Generate all available plots.
        
        Parameters
        ----------
        save_dir : str
            Directory to save plots
        prefix : str
            Prefix for plot filenames
        include_molecule_grid : bool, default=False
            Whether to include molecular structure grid visualization
        n_grid_samples : int, default=4
            Number of samples to show in molecule grid (only used if include_molecule_grid=True)
        show_tanimoto_on_grid : bool, default=False
            Whether to show Tanimoto similarity scores on the molecule grid
        """
        os.makedirs(save_dir, exist_ok=True)
        
        # Ensure metrics are computed
        if not self._computed_metrics:
            self.compute_all_metrics()
        
        # Generate plots based on available data
        if 'pred_selfies' in self.df.columns:
            # Comparison plots
            if 'token_accuracy' in self.df.columns:
                plot_token_accuracy_hist(
                    self.df['token_accuracy'],
                    save_dir,
                    f"{prefix}token_accuracy_hist.png" if prefix else "token_accuracy_hist.png"
                )
            
            if 'tanimoto_similarity' in self.df.columns:
                plot_tanimoto_similarity_hist(
                    self.df['tanimoto_similarity'],
                    save_dir,
                    f"{prefix}tanimoto_similarity_hist.png" if prefix else "tanimoto_similarity_hist.png"
                )
            
            if 'real_mol_wt' in self.df.columns and 'pred_mol_wt' in self.df.columns:
                plot_mw_comparison(
                    self.df,
                    save_dir=save_dir,
                    filename=f"{prefix}Pred_MW_vs_Target_MW.png" if prefix else "Pred_MW_vs_Target_MW.png"
                )
            
            if 'real_hdi' in self.df.columns and 'pred_hdi' in self.df.columns:
                plot_hdi_comparison(
                    self.df,
                    save_dir=save_dir,
                    filename=f"{prefix}Pred_HDI_vs_Target_HDI.png" if prefix else "Pred_HDI_vs_Target_HDI.png"
                )
            
            # Optional molecule grid
            if include_molecule_grid:
                self.plot_molecule_grid(
                    n_samples=n_grid_samples,
                    save_dir=save_dir,
                    filename=f"{prefix}molecule_comparison_grid.png" if prefix else "molecule_comparison_grid.png",
                    show_tanimoto=show_tanimoto_on_grid
                )
    
    def get_dataframe(self) -> pd.DataFrame:
        """Return the underlying DataFrame with all computed values."""
        return self.df.copy()
    
    def summary(self) -> str:
        """Generate a text summary of the analysis."""
        if not self._computed_metrics:
            self.compute_all_metrics()
        
        lines = ["SELFIES Analysis Summary", "=" * 50]
        lines.append(f"Number of molecules: {len(self.df)}")
        
        if 'mean_hdi' in self._computed_metrics:
            lines.append(f"Mean HDI: {self._computed_metrics['mean_hdi']:.3f}")
        
        if 'mean_mol_weight' in self._computed_metrics:
            lines.append(f"Mean Molecular Weight: {self._computed_metrics['mean_mol_weight']:.2f}")
        
        if 'mean_token_accuracy' in self._computed_metrics:
            lines.append(f"Mean Token Accuracy: {self._computed_metrics['mean_token_accuracy']:.3f}")
        
        if 'mean_tanimoto_similarity' in self._computed_metrics:
            lines.append(f"Mean Tanimoto Similarity: {self._computed_metrics['mean_tanimoto_similarity']:.3f}")
        
        return "\n".join(lines)
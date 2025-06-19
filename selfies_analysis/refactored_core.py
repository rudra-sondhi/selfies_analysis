# selfies_analysis/refactored_core.py
"""
Refactored core analyzer class for SELFIES analysis using modular architecture.
"""

import os
from typing import Union, Tuple, List, Dict, Optional, Any
import pandas as pd

# Try to import RDKit Mol type
try:
    from rdkit import Chem
    RDKit_Mol = Chem.Mol
    HAS_RDKIT = True
except ImportError:
    # Create a placeholder type if RDKit is not available
    class RDKit_Mol:
        pass
    HAS_RDKIT = False

from .input_handler import SELFIESInputHandler
from .metrics_calculator import SELFIESMetricsCalculator
from .summary_reporter import SELFIESSummaryReporter

# Import visualizer only if RDKit is available
if HAS_RDKIT:
    from .visualizer import SELFIESVisualizer


class SELFIESAnalyzer:
    """
    Refactored analyzer for SELFIES molecular representations using modular architecture.
    
    This class serves as a facade that orchestrates the functionality of specialized
    components for input handling, metrics calculation, visualization, and reporting.
    
    Components:
    - SELFIESInputHandler: Normalizes and validates input data
    - SELFIESMetricsCalculator: Computes molecular and SELFIES metrics
    - SELFIESVisualizer: Generates plots and visualizations
    - SELFIESSummaryReporter: Creates summaries and reports
    """
    
    def __init__(self, selfies_input: Union[str, Tuple[str, str], List[Tuple[str, str]], 
                                          List[str], List[List[str]], pd.DataFrame]):
        """
        Initialize the analyzer with input data.
        
        Parameters
        ----------
        selfies_input : str, tuple, list of tuples, list of strings, list of lists, or DataFrame
            Input data in various supported formats. See SELFIESInputHandler for details.
        """
        # Initialize components
        self.input_handler = SELFIESInputHandler(selfies_input)
        self.metrics_calculator = SELFIESMetricsCalculator()
        
        # Initialize visualizer only if RDKit is available
        if HAS_RDKIT:
            self.visualizer = SELFIESVisualizer()
        else:
            self.visualizer = None
            
        self.summary_reporter = SELFIESSummaryReporter()
        
        # Get normalized pairs from input handler
        self.pairs = self.input_handler.normalized_pairs
        
        # Storage for computed metrics
        self.metrics = {}
        self._metrics_computed = False
    
    @property
    def has_predictions(self) -> bool:
        """Check if the dataset contains prediction data."""
        return self.input_handler.has_predictions()
    
    @property
    def input_format(self) -> str:
        """Get the detected input format."""
        return self.input_handler.get_input_format()
    
    @property
    def input_statistics(self) -> Dict[str, Any]:
        """Get basic statistics about the input data."""
        return self.input_handler.get_statistics()
    
    def compute_all_metrics(self) -> Dict[str, Any]:
        """
        Compute all available metrics for the SELFIES pairs.
        
        Returns
        -------
        Dict[str, Any]
            Dictionary containing all computed metrics and statistics
        """
        self.metrics = self.metrics_calculator.compute_all(self.pairs)
        self._metrics_computed = True
        return self.metrics
    
    def get_hdi(self) -> Union[List[float], float, Dict[str, Union[List[float], float]]]:
        """
        Get HDI values for the molecules.
        
        Returns
        -------
        Union[List[float], float, Dict[str, Union[List[float], float]]]
            HDI values. Format depends on whether predictions are available.
        """
        if not self._metrics_computed:
            self.compute_all_metrics()
        
        if self.has_predictions:
            real_hdi = [x for x in self.metrics.get('real_hdi', []) if not pd.isna(x)]
            pred_hdi = [x for x in self.metrics.get('pred_hdi', []) if not pd.isna(x)]
            
            result = {
                'real': real_hdi[0] if len(real_hdi) == 1 else real_hdi,
                'predicted': pred_hdi[0] if len(pred_hdi) == 1 else pred_hdi
            }
            return result
        else:
            hdi_values = [x for x in self.metrics.get('real_hdi', []) if not pd.isna(x)]
            return hdi_values[0] if len(hdi_values) == 1 else hdi_values
    
    def get_molecular_weight(self) -> Union[List[float], float, Dict[str, Union[List[float], float]]]:
        """
        Get molecular weight values for the molecules.
        
        Returns
        -------
        Union[List[float], float, Dict[str, Union[List[float], float]]]
            Molecular weight values. Format depends on whether predictions are available.
        """
        if not self._metrics_computed:
            self.compute_all_metrics()
        
        if self.has_predictions:
            real_mw = [x for x in self.metrics.get('real_mol_wt', []) if not pd.isna(x)]
            pred_mw = [x for x in self.metrics.get('pred_mol_wt', []) if not pd.isna(x)]
            
            result = {
                'real': real_mw[0] if len(real_mw) == 1 else real_mw,
                'predicted': pred_mw[0] if len(pred_mw) == 1 else pred_mw
            }
            return result
        else:
            mw_values = [x for x in self.metrics.get('real_mol_wt', []) if not pd.isna(x)]
            return mw_values[0] if len(mw_values) == 1 else mw_values
    
    def get_token_accuracy(self) -> Union[List[float], float]:
        """
        Get token accuracy values comparing real vs predicted SELFIES.
        
        Returns
        -------
        Union[List[float], float]
            Token accuracy values (0-1)
            
        Raises
        ------
        ValueError
            If no predicted SELFIES are available for comparison
        """
        if not self.has_predictions:
            raise ValueError("No predicted SELFIES available for token accuracy calculation")
        
        if not self._metrics_computed:
            self.compute_all_metrics()
        
        accuracies = [x for x in self.metrics.get('token_accuracy', []) if not pd.isna(x)]
        return accuracies[0] if len(accuracies) == 1 else accuracies
    
    def get_sequence_accuracy(self) -> Union[List[float], float]:
        """
        Get token accuracy values comparing real vs predicted SELFIES.
        
        Returns
        -------
        Union[List[float], float]
            Token accuracy values (0-1)
            
        Raises
        ------
        ValueError
            If no predicted SELFIES are available for comparison
        """
        if not self.has_predictions:
            raise ValueError("No predicted SELFIES available for token accuracy calculation")
        
        if not self._metrics_computed:
            self.compute_all_metrics()
        
        accuracies = [x for x in self.metrics.get('sequence_accuracy', []) if not pd.isna(x)]
        return accuracies[0] if len(accuracies) == 1 else accuracies
    
    def get_tanimoto_similarity(self) -> Union[List[float], float]:
        """
        Get Tanimoto similarity values comparing real vs predicted molecules.
        
        Returns
        -------
        Union[List[float], float]
            Tanimoto similarity values (0-1)
            
        Raises
        ------
        ValueError
            If no predicted SELFIES are available for comparison
        """
        if not self.has_predictions:
            raise ValueError("No predicted SELFIES available for Tanimoto similarity calculation")
        
        if not self._metrics_computed:
            self.compute_all_metrics()
        
        similarities = [x for x in self.metrics.get('tanimoto_similarity', []) if not pd.isna(x)]
        return similarities[0] if len(similarities) == 1 else similarities
    
    def get_smiles(self) -> Union[List[str], str, Dict[str, Union[List[str], str]]]:
        """
        Get SMILES representations for the molecules.
        
        Returns
        -------
        Union[List[str], str, Dict[str, Union[List[str], str]]]
            SMILES strings. Format depends on whether predictions are available.
        """
        # Get molecular objects dataframe which includes SMILES
        df = self.metrics_calculator._get_molecular_objects(self.pairs)
        
        if self.has_predictions:
            real_smiles = df['real_smiles'].dropna().tolist()
            pred_smiles = df['pred_smiles'].dropna().tolist()
            
            result = {
                'real': real_smiles[0] if len(real_smiles) == 1 else real_smiles,
                'predicted': pred_smiles[0] if len(pred_smiles) == 1 else pred_smiles
            }
            return result
        else:
            smiles = df['real_smiles'].dropna().tolist()
            return smiles[0] if len(smiles) == 1 else smiles
    
    def get_mol(self) -> Union[List[Optional[RDKit_Mol]], Optional[RDKit_Mol], Dict[str, Union[List[Optional[RDKit_Mol]], Optional[RDKit_Mol]]]]:
        """
        Get RDKit Mol objects for the molecules.
        
        Returns
        -------
        Union[List[Optional], Optional, Dict[str, Union[List[Optional], Optional]]]
            RDKit Mol objects. Format depends on whether predictions are available.
            None values included for molecules that couldn't be converted.
        """
        # Get molecular objects dataframe which includes Mol objects
        df = self.metrics_calculator._get_molecular_objects(self.pairs)
        
        if self.has_predictions:
            real_mols = df['real_mol'].tolist()
            pred_mols = df['pred_mol'].tolist()
            
            result = {
                'real': real_mols[0] if len(real_mols) == 1 else real_mols,
                'predicted': pred_mols[0] if len(pred_mols) == 1 else pred_mols
            }
            return result
        else:
            mols = df['real_mol'].tolist()
            return mols[0] if len(mols) == 1 else mols
    
    def plot_all(self, save_dir: str = 'plots', 
                prefix: str = '',
                include_molecule_grid: bool = False,
                n_grid_samples: int = 4,
                show_tanimoto_on_grid: bool = False) -> None:
        """
        Generate all available plots.
        
        Parameters
        ----------
        save_dir : str, default='plots'
            Directory to save plots
        prefix : str, default=''
            Prefix for plot filenames
        include_molecule_grid : bool, default=False
            Whether to include molecular structure grid visualization
        n_grid_samples : int, default=4
            Number of samples to show in molecule grid
        show_tanimoto_on_grid : bool, default=False
            Whether to show Tanimoto similarity scores on the molecule grid
        """
        if not HAS_RDKIT:
            raise ImportError("RDKit required for visualization functionality")
            
        if not self._metrics_computed:
            self.compute_all_metrics()
        
        self.visualizer.plot_all(
            metrics=self.metrics,
            pairs=self.pairs,
            save_dir=save_dir,
            prefix=prefix,
            include_molecule_grid=include_molecule_grid,
            n_grid_samples=n_grid_samples,
            show_tanimoto_on_grid=show_tanimoto_on_grid
        )
    
    def plot_molecule_grid(self, n_samples: int = 4,
                          mol_size: Tuple[int, int] = (300, 300),
                          save_dir: str = 'plots',
                          filename: str = 'molecule_comparison_grid.png',
                          title: Optional[str] = None,
                          show_indices: bool = True,
                          show_tanimoto: bool = False,
                          rows: Optional[int] = None) -> None:
        """
        Plot a grid comparing real vs predicted molecules.
        
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
            Custom title for the plot
        show_indices : bool, default=True
            Whether to show sample indices in the grid
        show_tanimoto : bool, default=False
            Whether to display Tanimoto similarity scores for each pair
        rows : int, optional
            Number of rows to arrange the molecule pairs
        """
        if not HAS_RDKIT:
            raise ImportError("RDKit required for molecular structure grid visualization")
            
        self.visualizer.plot_molecule_grid(
            pairs=self.pairs,
            n_samples=n_samples,
            mol_size=mol_size,
            save_dir=save_dir,
            filename=filename,
            title=title,
            show_indices=show_indices,
            show_tanimoto=show_tanimoto,
            rows=rows
        )
    
    def summary(self) -> str:
        """
        Generate a text summary of the analysis.
        
        Returns
        -------
        str
            Human-readable summary of the analysis results
        """
        if not self._metrics_computed:
            self.compute_all_metrics()
        
        return self.summary_reporter.generate_summary(self.metrics)
    
    def get_dataframe(self) -> pd.DataFrame:
        """
        Return the analysis results as a pandas DataFrame.
        
        Returns
        -------
        pd.DataFrame
            DataFrame containing SELFIES pairs and computed metrics
        """
        if not self._metrics_computed:
            self.compute_all_metrics()
        
        return self.summary_reporter.to_dataframe(self.metrics, self.pairs)
    
    def save_results(self, save_dir: str = 'results',
                    prefix: str = '',
                    include_dataframe: bool = True,
                    include_text_summary: bool = True,
                    include_plots: bool = False,
                    include_molecule_grid: bool = False) -> Dict[str, str]:
        """
        Save analysis results to files.
        
        Parameters
        ----------
        save_dir : str, default='results'
            Directory to save files
        prefix : str, default=''
            Prefix for filenames
        include_dataframe : bool, default=True
            Whether to save results as CSV
        include_text_summary : bool, default=True
            Whether to save text summary
        include_plots : bool, default=False
            Whether to generate and save plots
        include_molecule_grid : bool, default=False
            Whether to include molecule grid in plots (only if include_plots=True)
            
        Returns
        -------
        Dict[str, str]
            Dictionary with paths to saved files
        """
        if not self._metrics_computed:
            self.compute_all_metrics()
        
        # Save summary and dataframe
        saved_files = self.summary_reporter.save_summary(
            metrics=self.metrics,
            pairs=self.pairs,
            save_dir=save_dir,
            prefix=prefix,
            include_dataframe=include_dataframe,
            include_text_summary=include_text_summary
        )
        
        # Optionally save plots
        if include_plots:
            if not HAS_RDKIT:
                print("Warning: RDKit not available. Skipping plot generation.")
            else:
                plots_dir = os.path.join(save_dir, 'plots')
                self.plot_all(
                    save_dir=plots_dir,
                    prefix=prefix,
                    include_molecule_grid=include_molecule_grid
                )
                saved_files['plots_directory'] = plots_dir
        
        return saved_files
    
    def get_available_plots(self) -> List[str]:
        """
        Get a list of plots that can be generated based on available data.
        
        Returns
        -------
        List[str]
            List of plot names that can be generated
        """
        if not HAS_RDKIT:
            return []
            
        if not self._metrics_computed:
            self.compute_all_metrics()
        
        return self.visualizer.get_available_plots(self.metrics, self.pairs)
    
    def get_summary_statistics(self) -> Dict[str, Any]:
        """
        Get key summary statistics.
        
        Returns
        -------
        Dict[str, Any]
            Dictionary with key summary statistics
        """
        if not self._metrics_computed:
            self.compute_all_metrics()
        
        return self.summary_reporter.get_summary_statistics(self.metrics)
    
    def clear_cache(self) -> None:
        """Clear internal caches to free memory."""
        self.metrics_calculator.clear_cache()
        self.metrics = {}
        self._metrics_computed = False
    
    def __repr__(self) -> str:
        """String representation of the analyzer."""
        stats = self.input_statistics
        return (f"SELFIESAnalyzer("
                f"molecules={stats['total_pairs']}, "
                f"with_predictions={stats['pairs_with_predictions']}, "
                f"format='{stats['input_format']}')")
    
    def __len__(self) -> int:
        """Return the number of molecule pairs."""
        return len(self.pairs)


# Backward compatibility - keep the original interface methods
class LegacySELFIESAnalyzer(SELFIESAnalyzer):
    """
    Legacy wrapper that maintains backward compatibility with the original API.
    
    This class provides the same interface as the original SELFIESAnalyzer
    while using the new modular architecture internally.
    """
    
    def compute_all_metrics(self) -> Dict[str, any]:
        """Legacy method that returns metrics in the original format."""
        metrics = super().compute_all_metrics()
        
        # Convert to legacy format for backward compatibility
        legacy_metrics = {}
        
        # HDI values
        if 'real_hdi' in metrics:
            hdi_values = self.get_hdi()
            if isinstance(hdi_values, dict):
                legacy_metrics['hdi_values'] = hdi_values['real']
            else:
                legacy_metrics['hdi_values'] = hdi_values if isinstance(hdi_values, list) else [hdi_values]
            
            if 'real_hdi_mean' in metrics:
                legacy_metrics['mean_hdi'] = metrics['real_hdi_mean']
        
        # Molecular weights
        if 'real_mol_wt' in metrics:
            mw_values = self.get_molecular_weight()
            if isinstance(mw_values, dict):
                legacy_metrics['mol_weights'] = mw_values['real']
            else:
                legacy_metrics['mol_weights'] = mw_values if isinstance(mw_values, list) else [mw_values]
            
            if 'real_mol_wt_mean' in metrics:
                legacy_metrics['mean_mol_weight'] = metrics['real_mol_wt_mean']
        
        # Prediction metrics
        if self.has_predictions:
            if 'token_accuracy' in metrics:
                accuracies = self.get_token_accuracy()
                legacy_metrics['token_accuracy'] = accuracies if isinstance(accuracies, list) else [accuracies]
                if 'token_accuracy_mean' in metrics:
                    legacy_metrics['mean_token_accuracy'] = metrics['token_accuracy_mean']
            
            if 'tanimoto_similarity' in metrics:
                similarities = self.get_tanimoto_similarity()
                legacy_metrics['tanimoto_similarity'] = similarities if isinstance(similarities, list) else [similarities]
                if 'tanimoto_similarity_mean' in metrics:
                    legacy_metrics['mean_tanimoto_similarity'] = metrics['tanimoto_similarity_mean']
        
        return legacy_metrics
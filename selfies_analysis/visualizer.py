# selfies_analysis/visualizer.py
"""
Visualization functions for SELFIES analysis.
"""

import os
import pandas as pd
from typing import List, Tuple, Optional, Dict, Any

from .plots import (
    plot_token_accuracy_hist,
    plot_tanimoto_similarity_hist,
    plot_mw_comparison,
    plot_hdi_comparison,
    plot_molecule_grid
)


class SELFIESVisualizer:
    """
    Handles generation of plots from metrics or molecule pairs.
    
    Provides methods to create various visualizations including histograms,
    comparison plots, and molecular structure grids.
    """
    
    def __init__(self):
        """Initialize the visualizer."""
        pass
    
    def plot_token_accuracy(self, metrics: Dict[str, Any], 
                          save_dir: str = 'plots',
                          filename: str = 'token_accuracy_hist.png') -> None:
        """
        Plot histogram of token accuracy values.
        
        Parameters
        ----------
        metrics : Dict[str, Any]
            Metrics dictionary containing 'token_accuracy' key
        save_dir : str, default='plots'
            Directory to save the plot
        filename : str, default='token_accuracy_hist.png'
            Filename for the saved plot
            
        Raises
        ------
        ValueError
            If token_accuracy data is not available in metrics
        """
        if 'token_accuracy' not in metrics:
            raise ValueError("Token accuracy data not found in metrics")
        
        token_accuracy_data = [x for x in metrics['token_accuracy'] if not pd.isna(x)]
        
        if not token_accuracy_data:
            raise ValueError("No valid token accuracy data available")
        
        plot_token_accuracy_hist(token_accuracy_data, save_dir, filename)
    
    def plot_tanimoto_distribution(self, metrics: Dict[str, Any],
                                 save_dir: str = 'plots',
                                 filename: str = 'tanimoto_similarity_hist.png') -> None:
        """
        Plot histogram of Tanimoto similarity values.
        
        Parameters
        ----------
        metrics : Dict[str, Any]
            Metrics dictionary containing 'tanimoto_similarity' key
        save_dir : str, default='plots'
            Directory to save the plot
        filename : str, default='tanimoto_similarity_hist.png'
            Filename for the saved plot
            
        Raises
        ------
        ValueError
            If tanimoto_similarity data is not available in metrics
        """
        if 'tanimoto_similarity' not in metrics:
            raise ValueError("Tanimoto similarity data not found in metrics")
        
        tanimoto_data = [x for x in metrics['tanimoto_similarity'] if not pd.isna(x)]
        
        if not tanimoto_data:
            raise ValueError("No valid Tanimoto similarity data available")
        
        plot_tanimoto_similarity_hist(tanimoto_data, save_dir, filename)
    
    def plot_molecular_weight_comparison(self, metrics: Dict[str, Any],
                                       save_dir: str = 'plots',
                                       filename: str = 'Pred_MW_vs_Target_MW.png') -> None:
        """
        Plot molecular weight comparison (real vs predicted).
        
        Parameters
        ----------
        metrics : Dict[str, Any]
            Metrics dictionary containing molecular weight data
        save_dir : str, default='plots'
            Directory to save the plot
        filename : str, default='Pred_MW_vs_Target_MW.png'
            Filename for the saved plot
            
        Raises
        ------
        ValueError
            If molecular weight data is not available in metrics
        """
        required_keys = ['real_mol_wt', 'pred_mol_wt']
        if not all(key in metrics for key in required_keys):
            raise ValueError("Molecular weight data (real_mol_wt, pred_mol_wt) not found in metrics")
        
        # Create DataFrame for plotting function
        df_data = {
            'real_mol_wt': metrics['real_mol_wt'],
            'pred_mol_wt': metrics['pred_mol_wt']
        }
        df = pd.DataFrame(df_data)
        
        # Filter out NaN values
        df_clean = df.dropna()
        
        if len(df_clean) == 0:
            raise ValueError("No valid molecular weight pairs available for plotting")
        
        plot_mw_comparison(df_clean, save_dir=save_dir, filename=filename)
    
    def plot_hdi_comparison(self, metrics: Dict[str, Any],
                          save_dir: str = 'plots',
                          filename: str = 'Pred_HDI_vs_Target_HDI.png') -> None:
        """
        Plot HDI comparison (real vs predicted).
        
        Parameters
        ----------
        metrics : Dict[str, Any]
            Metrics dictionary containing HDI data
        save_dir : str, default='plots'
            Directory to save the plot
        filename : str, default='Pred_HDI_vs_Target_HDI.png'
            Filename for the saved plot
            
        Raises
        ------
        ValueError
            If HDI data is not available in metrics
        """
        required_keys = ['real_hdi', 'pred_hdi']
        if not all(key in metrics for key in required_keys):
            raise ValueError("HDI data (real_hdi, pred_hdi) not found in metrics")
        
        # Create DataFrame for plotting function
        df_data = {
            'real_hdi': metrics['real_hdi'],
            'pred_hdi': metrics['pred_hdi']
        }
        df = pd.DataFrame(df_data)
        
        # Filter out NaN values
        df_clean = df.dropna()
        
        if len(df_clean) == 0:
            raise ValueError("No valid HDI pairs available for plotting")
        
        plot_hdi_comparison(df_clean, save_dir=save_dir, filename=filename)
    
    def plot_molecule_grid(self, pairs: List[Tuple[Optional[str], Optional[str]]],
                         n_samples: int = 4,
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
        pairs : List[Tuple[Optional[str], Optional[str]]]
            List of (real_selfies, pred_selfies) pairs
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
            
        Raises
        ------
        ValueError
            If no pairs with predictions are available
        """
        # Filter pairs with both real and predicted SELFIES
        valid_pairs = [(real, pred) for real, pred in pairs 
                      if real is not None and pred is not None]
        
        if not valid_pairs:
            raise ValueError("No pairs with predictions available for molecule grid visualization")
        
        plot_molecule_grid(
            selfies_pairs=valid_pairs,
            n_samples=n_samples,
            mol_size=mol_size,
            save_dir=save_dir,
            filename=filename,
            title=title,
            show_indices=show_indices,
            show_tanimoto=show_tanimoto,
            rows=rows
        )
    
    def plot_all(self, metrics: Dict[str, Any], 
                pairs: List[Tuple[Optional[str], Optional[str]]],
                save_dir: str = 'plots',
                prefix: str = '',
                include_molecule_grid: bool = False,
                n_grid_samples: int = 4,
                show_tanimoto_on_grid: bool = False) -> None:
        """
        Generate all available plots based on metrics and pairs.
        
        Parameters
        ----------
        metrics : Dict[str, Any]
            Computed metrics dictionary
        pairs : List[Tuple[Optional[str], Optional[str]]]
            List of (real_selfies, pred_selfies) pairs
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
        os.makedirs(save_dir, exist_ok=True)
        
        plots_generated = []
        
        # Generate comparison plots if predictions are available
        has_predictions = metrics.get('has_predictions', False)
        
        if has_predictions:
            # Token accuracy histogram
            if 'token_accuracy' in metrics:
                try:
                    filename = f"{prefix}token_accuracy_hist.png" if prefix else "token_accuracy_hist.png"
                    self.plot_token_accuracy(metrics, save_dir, filename)
                    plots_generated.append(filename)
                except Exception as e:
                    print(f"Warning: Could not generate token accuracy plot: {e}")
            
            # Tanimoto similarity histogram
            if 'tanimoto_similarity' in metrics:
                try:
                    filename = f"{prefix}tanimoto_similarity_hist.png" if prefix else "tanimoto_similarity_hist.png"
                    self.plot_tanimoto_distribution(metrics, save_dir, filename)
                    plots_generated.append(filename)
                except Exception as e:
                    print(f"Warning: Could not generate Tanimoto similarity plot: {e}")
            
            # Molecular weight comparison
            if 'real_mol_wt' in metrics and 'pred_mol_wt' in metrics:
                try:
                    filename = f"{prefix}Pred_MW_vs_Target_MW.png" if prefix else "Pred_MW_vs_Target_MW.png"
                    self.plot_molecular_weight_comparison(metrics, save_dir, filename)
                    plots_generated.append(filename)
                except Exception as e:
                    print(f"Warning: Could not generate molecular weight comparison plot: {e}")
            
            # HDI comparison
            if 'real_hdi' in metrics and 'pred_hdi' in metrics:
                try:
                    filename = f"{prefix}Pred_HDI_vs_Target_HDI.png" if prefix else "Pred_HDI_vs_Target_HDI.png"
                    self.plot_hdi_comparison(metrics, save_dir, filename)
                    plots_generated.append(filename)
                except Exception as e:
                    print(f"Warning: Could not generate HDI comparison plot: {e}")
            
            # Optional molecule grid
            if include_molecule_grid:
                try:
                    filename = f"{prefix}molecule_comparison_grid.png" if prefix else "molecule_comparison_grid.png"
                    self.plot_molecule_grid(
                        pairs=pairs,
                        n_samples=n_grid_samples,
                        save_dir=save_dir,
                        filename=filename,
                        show_tanimoto=show_tanimoto_on_grid
                    )
                    plots_generated.append(filename)
                except Exception as e:
                    print(f"Warning: Could not generate molecule grid plot: {e}")
        
        else:
            print("No prediction data available. Skipping comparison plots.")
        
        if plots_generated:
            print(f"Generated {len(plots_generated)} plots in {save_dir}:")
            for plot in plots_generated:
                print(f"  - {plot}")
        else:
            print("No plots were generated.")
    
    def get_available_plots(self, metrics: Dict[str, Any], 
                          pairs: List[Tuple[Optional[str], Optional[str]]]) -> List[str]:
        """
        Get a list of plots that can be generated based on available data.
        
        Parameters
        ----------
        metrics : Dict[str, Any]
            Computed metrics dictionary
        pairs : List[Tuple[Optional[str], Optional[str]]]
            List of (real_selfies, pred_selfies) pairs
            
        Returns
        -------
        List[str]
            List of plot names that can be generated
        """
        available_plots = []
        
        has_predictions = metrics.get('has_predictions', False)
        
        if has_predictions:
            if 'token_accuracy' in metrics:
                available_plots.append('token_accuracy_histogram')
            
            if 'tanimoto_similarity' in metrics:
                available_plots.append('tanimoto_similarity_histogram')
            
            if 'real_mol_wt' in metrics and 'pred_mol_wt' in metrics:
                available_plots.append('molecular_weight_comparison')
            
            if 'real_hdi' in metrics and 'pred_hdi' in metrics:
                available_plots.append('hdi_comparison')
            
            # Check if we have valid pairs for molecule grid
            valid_pairs = [(real, pred) for real, pred in pairs 
                          if real is not None and pred is not None]
            if valid_pairs:
                available_plots.append('molecule_grid')
        
        return available_plots
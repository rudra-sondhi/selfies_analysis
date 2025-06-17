# selfies_analysis/plots.py
"""
Plotting functions for SELFIES analysis.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import random

from typing import Union, Optional, List, Tuple
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D


def plot_molecule_grid(selfies_pairs: Union[List[Tuple[str, str]], pd.DataFrame],
                      n_samples: int = 4,
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
    selfies_pairs : List[Tuple[str, str]] or pandas.DataFrame
        SELFIES pairs as (real, predicted) tuples or DataFrame with 'real_selfies' and 'pred_selfies' columns
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
        Number of rows to arrange the molecule pairs. If None, uses automatic layout:
        - Single pair: 1 row, 2 columns
        - Multiple pairs: n_samples rows, 2 columns
        If specified, arranges pairs in a rows x cols grid where cols = 2 * ceil(n_samples / rows)
        
    Returns
    -------
    None
        Saves the plot to the specified directory
        
    Examples
    --------
    >>> pairs = [("[C][C][O]", "[C][O][C]"), ("[C][=C][C]", "[C][C][=C]")]
    >>> plot_molecule_grid(pairs, n_samples=2, save_dir='results', show_tanimoto=True)
    
    >>> # Custom layout: 2 rows with 4 pairs (2 pairs per row)
    >>> plot_molecule_grid(pairs, n_samples=4, rows=2, save_dir='results')
    
    >>> # Using with DataFrame
    >>> df = pd.DataFrame({'real_selfies': ['[C][C][O]'], 'pred_selfies': ['[C][O][C]']})
    >>> plot_molecule_grid(df, n_samples=1, show_tanimoto=True)
    """
    from .metrics import selfies_to_smiles, smiles_to_mol, calculate_tanimoto_similarity
    
    os.makedirs(save_dir, exist_ok=True)
    
    # Convert input to standardized format
    if isinstance(selfies_pairs, pd.DataFrame):
        if 'real_selfies' not in selfies_pairs.columns or 'pred_selfies' not in selfies_pairs.columns:
            raise ValueError("DataFrame must contain 'real_selfies' and 'pred_selfies' columns")
        pairs = list(zip(selfies_pairs['real_selfies'], selfies_pairs['pred_selfies']))
    else:
        pairs = selfies_pairs
    
    # Validate inputs
    if len(pairs) == 0:
        raise ValueError("No SELFIES pairs provided")
    
    n_samples = min(n_samples, len(pairs))
    if n_samples <= 0:
        raise ValueError("n_samples must be positive")
    
    # Sample pairs (take first n_samples)
    selected_pairs = random.sample(pairs, n_samples)
    
    # Convert SELFIES to molecules
    mol_pairs = []
    valid_indices = []
    tanimoto_scores = []
    
    for i, (real_selfies, pred_selfies) in enumerate(selected_pairs):
        try:
            real_smiles = selfies_to_smiles(real_selfies)
            pred_smiles = selfies_to_smiles(pred_selfies)
            
            real_mol = smiles_to_mol(real_smiles) if real_smiles else None
            pred_mol = smiles_to_mol(pred_smiles) if pred_smiles else None
            
            if real_mol is not None and pred_mol is not None:
                mol_pairs.append((real_mol, pred_mol))
                valid_indices.append(i)
                
                # Calculate Tanimoto similarity if requested
                if show_tanimoto:
                    tanimoto_score = calculate_tanimoto_similarity(real_mol, pred_mol)
                    tanimoto_scores.append(tanimoto_score)
                else:
                    tanimoto_scores.append(None)
            else:
                print(f"Warning: Could not convert pair {i} to valid molecules. Skipping.")
                
        except Exception as e:
            print(f"Warning: Error processing pair {i}: {e}. Skipping.")
    
    if not mol_pairs:
        raise ValueError("No valid molecule pairs could be generated from the input SELFIES")
    
    n_valid = len(mol_pairs)
    
    # Calculate grid dimensions
    if rows is None:
        # Automatic layout (original behavior)
        if n_valid == 1:
            grid_rows, grid_cols = 1, 2  # 1 row, 2 columns (real, pred)
        else:
            # For multiple samples, arrange in rows with 2 columns per sample
            grid_rows = n_valid
            grid_cols = 2
    else:
        # Custom layout with specified number of rows
        if rows <= 0:
            raise ValueError("Number of rows must be positive")
        
        grid_rows = rows
        pairs_per_row = np.ceil(n_valid / rows).astype(int)
        grid_cols = pairs_per_row * 2  # Each pair needs 2 columns (real, pred)
        
        # Ensure we don't exceed available pairs
        if grid_rows * pairs_per_row < n_valid:
            # Add extra row if needed
            grid_rows += 1
    
    # Create figure
    fig_width = grid_cols * (mol_size[0] / 100) + 1  # Convert pixels to inches (rough estimate)
    fig_height = grid_rows * (mol_size[1] / 100) + 2
    
    fig, axes = plt.subplots(grid_rows, grid_cols, figsize=(fig_width, fig_height), dpi=300)
    
    # Handle single subplot case
    if grid_rows == 1 and grid_cols == 1:
        axes = np.array([[axes]])
    elif grid_rows == 1:
        axes = axes.reshape(1, -1)
    elif grid_cols == 1:
        axes = axes.reshape(-1, 1)
    else:
        axes = axes.reshape(grid_rows, grid_cols)
    
    # Calculate pairs per row for custom layout
    if rows is not None:
        pairs_per_row = np.ceil(n_valid / rows).astype(int)
    else:
        pairs_per_row = 1  # Original layout: 1 pair per row
    
    # Generate molecule images
    pair_idx = 0
    for row in range(grid_rows):
        for pair_in_row in range(pairs_per_row):
            if pair_idx >= n_valid:
                # No more pairs to display, hide remaining subplots
                for remaining_col in range(pair_in_row * 2, grid_cols):
                    if remaining_col < grid_cols:
                        axes[row, remaining_col].axis('off')
                break
            
            real_mol, pred_mol = mol_pairs[pair_idx]
            sample_idx = valid_indices[pair_idx]
            tanimoto_score = tanimoto_scores[pair_idx]
            
            # Calculate column positions for this pair
            real_col = pair_in_row * 2
            pred_col = pair_in_row * 2 + 1
            
            # Draw real molecule
            real_img = Draw.MolToImage(real_mol, size=mol_size)
            axes[row, real_col].imshow(real_img)
            axes[row, real_col].axis('off')
            
            # Only show "Real"/"Predicted" labels on top row
            if row == 0:
                if show_indices:
                    axes[row, real_col].set_title(f'Real (Sample {sample_idx})', fontsize=12, fontweight='bold')
                    axes[row, pred_col].set_title(f'Predicted (Sample {sample_idx})', fontsize=12, fontweight='bold')
                else:
                    axes[row, real_col].set_title('Real', fontsize=12, fontweight='bold')
                    axes[row, pred_col].set_title('Predicted', fontsize=12, fontweight='bold')
            else:
                # For other rows, only show sample indices if requested
                if show_indices:
                    axes[row, real_col].set_title(f'Sample {sample_idx}', fontsize=10)
                    # No title for predicted column in non-top rows
            
            # Draw predicted molecule
            pred_img = Draw.MolToImage(pred_mol, size=mol_size)
            axes[row, pred_col].imshow(pred_img)
            axes[row, pred_col].axis('off')
            

            # Add Tanimoto similarity score if requested
            if show_tanimoto and tanimoto_score is not None:
                # Add text box in bottom right corner of the predicted molecule subplot
                axes[row, pred_col].text(0.95, 0.05, f'Tanimoto: {tanimoto_score:.3f}', 
                                transform=axes[row, pred_col].transAxes,
                                fontsize=10, fontweight='bold',
                                verticalalignment='bottom',
                                horizontalalignment='right',
                                bbox=dict(boxstyle='round,pad=0.3', 
                                        facecolor='white', 
                                        edgecolor='black',
                                        alpha=0.8))
            
            # Add border around each pair
            for col_offset in [0, 1]:  # real and pred columns
                col = real_col + col_offset
                rect = Rectangle((0, 0), 1, 1, linewidth=2, edgecolor='black', 
                               facecolor='none', transform=axes[row, col].transAxes)
                axes[row, col].add_patch(rect)
            
            pair_idx += 1
    
    # Set overall title
    # if title is None:
    #     title = f'Molecular Structure Comparison\n({n_valid} Sample{"s" if n_valid != 1 else ""})'
    
    fig.suptitle(title, fontsize=16, fontweight='bold', y=0.95)
    
    # Adjust layout
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)  # Make room for title
    
    # Save plot
    save_path = os.path.join(save_dir, filename)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Molecular grid plot saved to: {save_path}")
    print(f"Displayed {n_valid} valid molecule pairs out of {n_samples} requested")


def plot_token_accuracy_hist(token_accuracy: Union[list, pd.Series], 
                           save_dir: str, 
                           filename: str = 'token_accuracy_hist.png') -> None:
    """
    Plot histogram of token accuracy values.
    
    Parameters
    ----------
    token_accuracy : list or pandas.Series
        Token accuracy values (0 to 1)
    save_dir : str
        Directory to save the plot
    filename : str
        Filename for the saved plot
    """
    os.makedirs(save_dir, exist_ok=True)
    save_path = os.path.join(save_dir, filename)
    
    plt.figure(figsize=(10, 6), dpi=300)
    plt.hist(token_accuracy, bins=20, range=(0, 1), color='skyblue', edgecolor='black')
    
    ticks = [i/10 for i in range(0, 11)]
    plt.xticks(ticks, [f"{t:.1f}" for t in ticks], fontsize=14)
    plt.yticks(fontsize=14)
    
    plt.grid(False)
    plt.xlabel('Token Accuracy', fontsize=18, labelpad=15)
    plt.ylabel('Count', fontsize=18, labelpad=15)
    
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()


def plot_tanimoto_similarity_hist(similarity_series: Union[list, pd.Series], 
                                save_dir: str,
                                filename: str = 'tanimoto_similarity_hist.png') -> None:
    """
    Plot histogram of Tanimoto similarity values.
    
    Parameters
    ----------
    similarity_series : list or pandas.Series
        Tanimoto similarity values (0 to 1)
    save_dir : str
        Directory to save the plot
    filename : str
        Filename for the saved plot
    """
    os.makedirs(save_dir, exist_ok=True)
    save_path = os.path.join(save_dir, filename)
    
    plt.figure(figsize=(10, 6), dpi=300)
    plt.hist(similarity_series, bins=20, range=(0, 1),
             color='skyblue', edgecolor='black')
    
    ticks = [i/10 for i in range(0, 11)]
    plt.xticks(ticks, [f"{t:.1f}" for t in ticks], fontsize=14)
    plt.yticks(fontsize=14)
    
    plt.grid(False)
    plt.xlabel('Tanimoto Similarity', fontsize=18, labelpad=15)
    plt.ylabel('Count', fontsize=18, labelpad=15)
    
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()


def plot_hdi_comparison(df: pd.DataFrame,
                       real_hdi_col: str = 'real_hdi',
                       pred_hdi_col: str = 'pred_hdi',
                       save_dir: str = 'plots',
                       filename: str = 'Pred_HDI_vs_Target_HDI.png') -> None:
    """
    Create hexbin plots comparing real vs predicted HDI values.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing HDI columns
    real_hdi_col : str
        Column name for real HDI values
    pred_hdi_col : str
        Column name for predicted HDI values
    save_dir : str
        Directory to save plots
    filename : str
        Base filename for saved plots
    """
    os.makedirs(save_dir, exist_ok=True)
    sub = df[[real_hdi_col, pred_hdi_col]].dropna()
    
    if len(sub) == 0:
        return
    
    X = sub[real_hdi_col].values.reshape(-1, 1)
    y = sub[pred_hdi_col].values
    
    model = LinearRegression().fit(X, y)
    trend = model.predict(X)
    r2 = r2_score(y, trend)
    
    for suffix, annotate in [('', True), ('_bare', False)]:
        fig = plt.figure(figsize=(8, 8), dpi=300)
        grid = plt.GridSpec(4, 4, hspace=0.0, wspace=0.0)
        
        ax_main = fig.add_subplot(grid[1:4, 0:3])
        ax_main.hexbin(sub[real_hdi_col], sub[pred_hdi_col],
                      gridsize=30, cmap='Reds', mincnt=1)
        
        if annotate:
            ax_main.plot(sub[real_hdi_col], trend, 'k-', lw=2)
            ax_main.text(0.05, 0.95, f'$R^2 = {r2:.2f}$',
                        transform=ax_main.transAxes,
                        fontsize=18, va='top', fontweight='bold')
        
        ax_top = fig.add_subplot(grid[0, 0:3], sharex=ax_main)
        ax_top.hist(sub[real_hdi_col], bins=20, color='red', edgecolor='black')
        ax_top.axis('off')
        
        ax_right = fig.add_subplot(grid[1:4, 3], sharey=ax_main)
        ax_right.hist(sub[pred_hdi_col], bins=20, orientation='horizontal',
                     color='red', edgecolor='black')
        ax_right.axis('off')
        
        ax_main.set_xlabel('Target HDI', fontsize=18, fontweight='bold', labelpad=12)
        ax_main.set_ylabel('Predicted HDI', fontsize=18, fontweight='bold', labelpad=12)
        ax_main.tick_params(axis='both', labelsize=16)
        ax_top.spines['bottom'].set_visible(False)
        ax_right.spines['left'].set_visible(False)
        
        plt.tight_layout()
        out_name = filename.replace('.png', f'{suffix}.png')
        plt.savefig(os.path.join(save_dir, out_name))
        plt.close(fig)


def plot_mw_comparison(df: pd.DataFrame,
                      real_wt_col: str = 'real_mol_wt',
                      pred_wt_col: str = 'pred_mol_wt',
                      save_dir: str = 'plots',
                      filename: str = 'Pred_MW_vs_Target_MW.png') -> None:
    """
    Create hexbin plots comparing real vs predicted molecular weights.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing molecular weight columns
    real_wt_col : str
        Column name for real molecular weights
    pred_wt_col : str
        Column name for predicted molecular weights
    save_dir : str
        Directory to save plots
    filename : str
        Base filename for saved plots
    """
    os.makedirs(save_dir, exist_ok=True)
    sub = df[[real_wt_col, pred_wt_col]].dropna()
    
    if len(sub) == 0:
        return
    
    X = sub[real_wt_col].values.reshape(-1, 1)
    y = sub[pred_wt_col].values
    model = LinearRegression().fit(X, y)
    trendline = model.predict(X)
    r2 = r2_score(y, trendline)
    
    for suffix, annotate in [('', True), ('_bare', False)]:
        fig = plt.figure(figsize=(8, 8), dpi=300)
        grid = plt.GridSpec(4, 4, hspace=0.0, wspace=0.0)
        
        ax_main = fig.add_subplot(grid[1:4, 0:3])
        ax_main.hexbin(sub[real_wt_col], sub[pred_wt_col],
                      gridsize=30, cmap='Reds', mincnt=1)
        if annotate:
            ax_main.plot(sub[real_wt_col], trendline, 'k-', lw=2)
            ax_main.text(0.05, 0.95, f'$R^2 = {r2:.2f}$',
                        transform=ax_main.transAxes,
                        fontsize=18, va='top', fontweight='bold')
        
        ax_top = fig.add_subplot(grid[0, 0:3], sharex=ax_main)
        ax_top.hist(sub[real_wt_col], bins=20, color='red', edgecolor='black')
        ax_top.axis('off')
        
        ax_right = fig.add_subplot(grid[1:4, 3], sharey=ax_main)
        ax_right.hist(sub[pred_wt_col], bins=20, orientation='horizontal',
                     color='red', edgecolor='black')
        ax_right.axis('off')
        
        ax_main.set_xlabel('Target Molecular Weight', fontsize=18,
                          fontweight='bold', labelpad=12)
        ax_main.set_ylabel('Predicted Molecular Weight', fontsize=18,
                          fontweight='bold', labelpad=12)
        ax_main.tick_params(axis='both', labelsize=16)
        ax_top.spines['bottom'].set_visible(False)
        ax_right.spines['left'].set_visible(False)
        
        out_name = filename.replace('.png', f'{suffix}.png')
        plt.savefig(os.path.join(save_dir, out_name))
        plt.close(fig)
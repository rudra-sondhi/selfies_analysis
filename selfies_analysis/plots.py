# selfies_analysis/plots.py
"""
Plotting functions for SELFIES analysis.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from typing import Union, Optional


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
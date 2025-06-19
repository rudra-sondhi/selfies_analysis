# selfies_analysis/summary_reporter.py
"""
Summary and reporting functions for SELFIES analysis.
"""

import os
import pandas as pd
from typing import List, Tuple, Optional, Dict, Any
from datetime import datetime


class SELFIESSummaryReporter:
    """
    Handles generation of human-readable summaries and structured results.
    
    Provides methods to create text summaries, convert results to DataFrames,
    and save analysis reports.
    """
    
    def __init__(self):
        """Initialize the summary reporter."""
        pass
    
    def generate_summary(self, metrics: Dict[str, Any]) -> str:
        """
        Generate a human-readable summary of the analysis.
        
        Parameters
        ----------
        metrics : Dict[str, Any]
            Dictionary containing computed metrics and statistics
            
        Returns
        -------
        str
            Formatted text summary of the analysis
        """
        if not metrics:
            return "No metrics available for summary."
        
        lines = []
        lines.append("SELFIES Analysis Summary")
        lines.append("=" * 50)
        lines.append(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        lines.append("")
        
        # Basic information
        total_pairs = metrics.get('total_pairs', 0)
        pairs_with_predictions = metrics.get('pairs_with_predictions', 0)
        pairs_without_predictions = total_pairs - pairs_with_predictions
        
        lines.append("Dataset Information:")
        lines.append(f"  Total molecules: {total_pairs}")
        lines.append(f"  Pairs with predictions: {pairs_with_predictions}")
        lines.append(f"  Single molecules (no predictions): {pairs_without_predictions}")
        lines.append("")
        
        # HDI Analysis
        if 'real_hdi_mean' in metrics and metrics['real_hdi_mean'] is not None:
            lines.append("Hydrogen Deficiency Index (HDI) Analysis:")
            lines.append(f"  Real molecules - Mean HDI: {metrics['real_hdi_mean']:.3f}")
            if 'real_hdi_std' in metrics and metrics['real_hdi_std'] is not None:
                lines.append(f"  Real molecules - Std HDI: {metrics['real_hdi_std']:.3f}")
            
            if 'pred_hdi_mean' in metrics and metrics['pred_hdi_mean'] is not None:
                lines.append(f"  Predicted molecules - Mean HDI: {metrics['pred_hdi_mean']:.3f}")
                if 'pred_hdi_std' in metrics and metrics['pred_hdi_std'] is not None:
                    lines.append(f"  Predicted molecules - Std HDI: {metrics['pred_hdi_std']:.3f}")
                
                # HDI difference
                hdi_diff = metrics['pred_hdi_mean'] - metrics['real_hdi_mean']
                lines.append(f"  Mean HDI difference (Pred - Real): {hdi_diff:+.3f}")
            
            lines.append("")
        
        # Molecular Weight Analysis
        if 'real_mol_wt_mean' in metrics and metrics['real_mol_wt_mean'] is not None:
            lines.append("Molecular Weight Analysis:")
            lines.append(f"  Real molecules - Mean MW: {metrics['real_mol_wt_mean']:.2f} Da")
            if 'real_mol_wt_std' in metrics and metrics['real_mol_wt_std'] is not None:
                lines.append(f"  Real molecules - Std MW: {metrics['real_mol_wt_std']:.2f} Da")
            
            if 'pred_mol_wt_mean' in metrics and metrics['pred_mol_wt_mean'] is not None:
                lines.append(f"  Predicted molecules - Mean MW: {metrics['pred_mol_wt_mean']:.2f} Da")
                if 'pred_mol_wt_std' in metrics and metrics['pred_mol_wt_std'] is not None:
                    lines.append(f"  Predicted molecules - Std MW: {metrics['pred_mol_wt_std']:.2f} Da")
                
                # MW difference
                mw_diff = metrics['pred_mol_wt_mean'] - metrics['real_mol_wt_mean']
                lines.append(f"  Mean MW difference (Pred - Real): {mw_diff:+.2f} Da")
            
            lines.append("")
        
        # Prediction Quality Metrics
        if pairs_with_predictions > 0:
            lines.append("Prediction Quality Metrics:")
            
            # Token Accuracy
            if 'token_accuracy_mean' in metrics and metrics['token_accuracy_mean'] is not None:
                lines.append(f"  Mean Token Accuracy: {metrics['token_accuracy_mean']:.3f}")
                if 'token_accuracy_std' in metrics and metrics['token_accuracy_std'] is not None:
                    lines.append(f"  Token Accuracy Std: {metrics['token_accuracy_std']:.3f}")
                
                # Interpretation
                if metrics['token_accuracy_mean'] >= 0.9:
                    interpretation = "Excellent"
                elif metrics['token_accuracy_mean'] >= 0.8:
                    interpretation = "Good"
                elif metrics['token_accuracy_mean'] >= 0.7:
                    interpretation = "Fair"
                else:
                    interpretation = "Poor"
                lines.append(f"  Token Accuracy Quality: {interpretation}")

            if 'sequence_accuracy_mean' in metrics and metrics['sequence_accuracy_mean'] is not None:
                lines.append(f"  Mean Sequence Accuracy: {metrics['sequence_accuracy_mean']:.3f}")
                if 'sequence_accuracy_std' in metrics and metrics['sequence_accuracy_std'] is not None:
                    lines.append(f"  Sequence Accuracy Std: {metrics['sequence_accuracy_std']:.3f}")
                
            
            # Tanimoto Similarity
            if 'tanimoto_similarity_mean' in metrics and metrics['tanimoto_similarity_mean'] is not None:
                lines.append(f"  Mean Tanimoto Similarity: {metrics['tanimoto_similarity_mean']:.3f}")
                if 'tanimoto_similarity_std' in metrics and metrics['tanimoto_similarity_std'] is not None:
                    lines.append(f"  Tanimoto Similarity Std: {metrics['tanimoto_similarity_std']:.3f}")
                
                # Interpretation
                if metrics['tanimoto_similarity_mean'] >= 0.8:
                    interpretation = "Very Similar"
                elif metrics['tanimoto_similarity_mean'] >= 0.6:
                    interpretation = "Moderately Similar"
                elif metrics['tanimoto_similarity_mean'] >= 0.4:
                    interpretation = "Somewhat Similar"
                else:
                    interpretation = "Dissimilar"
                lines.append(f"  Molecular Similarity: {interpretation}")
            
            lines.append("")
        
        # Data Quality Assessment
        lines.append("Data Quality Assessment:")
        
        # Count valid vs invalid molecules
        if 'real_mol_wt' in metrics:
            valid_real = sum(1 for x in metrics['real_mol_wt'] if not pd.isna(x))
            invalid_real = len(metrics['real_mol_wt']) - valid_real
            lines.append(f"  Valid real molecules: {valid_real}/{len(metrics['real_mol_wt'])} ({100*valid_real/len(metrics['real_mol_wt']):.1f}%)")
            if invalid_real > 0:
                lines.append(f"  Invalid real molecules: {invalid_real}")
        
        if 'pred_mol_wt' in metrics and pairs_with_predictions > 0:
            valid_pred = sum(1 for x in metrics['pred_mol_wt'] if not pd.isna(x))
            invalid_pred = len(metrics['pred_mol_wt']) - valid_pred
            lines.append(f"  Valid predicted molecules: {valid_pred}/{len(metrics['pred_mol_wt'])} ({100*valid_pred/len(metrics['pred_mol_wt']):.1f}%)")
            if invalid_pred > 0:
                lines.append(f"  Invalid predicted molecules: {invalid_pred}")
        
        lines.append("")
        
        # Recommendations
        lines.append("Recommendations:")
        if pairs_with_predictions > 0:
            if 'token_accuracy_mean' in metrics and metrics['token_accuracy_mean'] is not None:
                if metrics['token_accuracy_mean'] < 0.8:
                    lines.append("  - Consider improving token-level prediction accuracy")
            
            if 'tanimoto_similarity_mean' in metrics and metrics['tanimoto_similarity_mean'] is not None:
                if metrics['tanimoto_similarity_mean'] < 0.6:
                    lines.append("  - Focus on improving molecular similarity of predictions")
        
        if 'real_mol_wt' in metrics:
            valid_real = sum(1 for x in metrics['real_mol_wt'] if not pd.isna(x))
            if valid_real / len(metrics['real_mol_wt']) < 0.95:
                lines.append("  - Review and clean input SELFIES data (some molecules are invalid)")
        
        return "\n".join(lines)
    
    def to_dataframe(self, metrics: Dict[str, Any], 
                    pairs: List[Tuple[Optional[str], Optional[str]]]) -> pd.DataFrame:
        """
        Convert metrics and pairs to a pandas DataFrame.
        
        Parameters
        ----------
        metrics : Dict[str, Any]
            Dictionary containing computed metrics
        pairs : List[Tuple[Optional[str], Optional[str]]]
            List of (real_selfies, pred_selfies) pairs
            
        Returns
        -------
        pd.DataFrame
            DataFrame with SELFIES pairs and computed metrics
        """
        data = []
        
        for i, (real_selfies, pred_selfies) in enumerate(pairs):
            row = {
                'index': i,
                'real_selfies': real_selfies,
                'pred_selfies': pred_selfies
            }
            
            # Add metric values if available
            metric_columns = [
                'real_hdi', 'pred_hdi', 'real_mol_wt', 'pred_mol_wt',
                'token_accuracy', 'tanimoto_similarity'
            ]
            
            for col in metric_columns:
                if col in metrics and i < len(metrics[col]):
                    row[col] = metrics[col][i]
                else:
                    row[col] = None
            
            data.append(row)
        
        df = pd.DataFrame(data)
        
        # Add metadata as attributes (not columns)
        for key, value in metrics.items():
            if key not in ['real_hdi', 'pred_hdi', 'real_mol_wt', 'pred_mol_wt',
                          'token_accuracy', 'tanimoto_similarity']:
                setattr(df, f'meta_{key}', value)
        
        return df
    
    def save_summary(self, metrics: Dict[str, Any], 
                    pairs: List[Tuple[Optional[str], Optional[str]]],
                    save_dir: str = 'results',
                    prefix: str = '',
                    include_dataframe: bool = True,
                    include_text_summary: bool = True) -> Dict[str, str]:
        """
        Save analysis summary to files.
        
        Parameters
        ----------
        metrics : Dict[str, Any]
            Dictionary containing computed metrics
        pairs : List[Tuple[Optional[str], Optional[str]]]
            List of (real_selfies, pred_selfies) pairs
        save_dir : str, default='results'
            Directory to save files
        prefix : str, default=''
            Prefix for filenames
        include_dataframe : bool, default=True
            Whether to save results as CSV
        include_text_summary : bool, default=True
            Whether to save text summary
            
        Returns
        -------
        Dict[str, str]
            Dictionary with saved file paths
        """
        os.makedirs(save_dir, exist_ok=True)
        saved_files = {}
        
        # Save text summary
        if include_text_summary:
            summary_text = self.generate_summary(metrics)
            summary_filename = f"{prefix}analysis_summary.txt" if prefix else "analysis_summary.txt"
            summary_path = os.path.join(save_dir, summary_filename)
            
            with open(summary_path, 'w', encoding='utf-8') as f:
                f.write(summary_text)
            
            saved_files['summary'] = summary_path
        
        # Save DataFrame as CSV
        if include_dataframe:
            df = self.to_dataframe(metrics, pairs)
            csv_filename = f"{prefix}analysis_results.csv" if prefix else "analysis_results.csv"
            csv_path = os.path.join(save_dir, csv_filename)
            
            # Save main dataframe
            df.to_csv(csv_path, index=False)
            saved_files['dataframe'] = csv_path
            
            # Save metadata separately
            metadata_filename = f"{prefix}analysis_metadata.txt" if prefix else "analysis_metadata.txt"
            metadata_path = os.path.join(save_dir, metadata_filename)
            
            with open(metadata_path, 'w', encoding='utf-8') as f:
                f.write("Analysis Metadata\n")
                f.write("=================\n\n")
                for attr in dir(df):
                    if attr.startswith('meta_'):
                        key = attr[5:]  # Remove 'meta_' prefix
                        value = getattr(df, attr)
                        f.write(f"{key}: {value}\n")
            
            saved_files['metadata'] = metadata_path
        
        return saved_files
    
    def get_summary_statistics(self, metrics: Dict[str, Any]) -> Dict[str, Any]:
        """
        Extract key summary statistics from metrics.
        
        Parameters
        ----------
        metrics : Dict[str, Any]
            Dictionary containing computed metrics
            
        Returns
        -------
        Dict[str, Any]
            Dictionary with key summary statistics
        """
        summary_stats = {}
        
        # Basic counts
        summary_stats['total_molecules'] = metrics.get('total_pairs', 0)
        summary_stats['molecules_with_predictions'] = metrics.get('pairs_with_predictions', 0)
        summary_stats['has_predictions'] = metrics.get('has_predictions', False)
        
        # HDI statistics
        if 'real_hdi_mean' in metrics:
            summary_stats['mean_real_hdi'] = metrics['real_hdi_mean']
        if 'pred_hdi_mean' in metrics:
            summary_stats['mean_pred_hdi'] = metrics['pred_hdi_mean']
        
        # Molecular weight statistics
        if 'real_mol_wt_mean' in metrics:
            summary_stats['mean_real_mol_wt'] = metrics['real_mol_wt_mean']
        if 'pred_mol_wt_mean' in metrics:
            summary_stats['mean_pred_mol_wt'] = metrics['pred_mol_wt_mean']
        
        # Prediction quality
        if 'token_accuracy_mean' in metrics:
            summary_stats['mean_token_accuracy'] = metrics['token_accuracy_mean']
        if 'tanimoto_similarity_mean' in metrics:
            summary_stats['mean_tanimoto_similarity'] = metrics['tanimoto_similarity_mean']
        
        # Data quality
        if 'real_mol_wt' in metrics:
            valid_real = sum(1 for x in metrics['real_mol_wt'] if not pd.isna(x))
            summary_stats['valid_real_molecules_fraction'] = valid_real / len(metrics['real_mol_wt'])
        
        if 'pred_mol_wt' in metrics and len(metrics['pred_mol_wt']) > 0:
            valid_pred = sum(1 for x in metrics['pred_mol_wt'] if not pd.isna(x))
            summary_stats['valid_pred_molecules_fraction'] = valid_pred / len(metrics['pred_mol_wt'])
        
        return summary_stats
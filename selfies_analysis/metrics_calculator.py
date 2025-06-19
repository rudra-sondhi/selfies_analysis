# selfies_analysis/metrics_calculator.py
"""
Batch metrics calculation for SELFIES analysis.
"""

import pandas as pd
import numpy as np
from typing import List, Tuple, Optional, Dict, Any

from .metrics import (
    calculate_hdi,
    calculate_token_accuracy,
    calculate_tanimoto_similarity,
    calculate_molecular_weight,
    selfies_to_smiles,
    smiles_to_mol,
    calculate_sequence_accuracy
)


class SELFIESMetricsCalculator:
    """
    Handles batch computation of metrics from SELFIES pairs.
    
    Computes various molecular and SELFIES-specific metrics including HDI,
    molecular weight, token accuracy, and Tanimoto similarity.
    """
    
    def __init__(self):
        """Initialize the metrics calculator."""
        self._cache = {}
    
    def _get_molecular_objects(self, pairs: List[Tuple[Optional[str], Optional[str]]]) -> pd.DataFrame:
        """
        Convert SELFIES pairs to molecular objects and cache results.
        
        Parameters
        ----------
        pairs : List[Tuple[Optional[str], Optional[str]]]
            List of (real_selfies, pred_selfies) pairs
            
        Returns
        -------
        pd.DataFrame
            DataFrame with SELFIES, SMILES, and Mol columns
        """
        cache_key = "molecular_objects"
        if cache_key in self._cache:
            return self._cache[cache_key].copy()
        
        data = []
        for real_selfies, pred_selfies in pairs:
            row = {
                'real_selfies': real_selfies,
                'pred_selfies': pred_selfies,
                'real_smiles': None,
                'pred_smiles': None,
                'real_mol': None,
                'pred_mol': None
            }
            
            # Convert real SELFIES
            if real_selfies is not None:
                row['real_smiles'] = selfies_to_smiles(real_selfies)
                if row['real_smiles'] is not None:
                    row['real_mol'] = smiles_to_mol(row['real_smiles'])
            
            # Convert predicted SELFIES
            if pred_selfies is not None:
                row['pred_smiles'] = selfies_to_smiles(pred_selfies)
                if row['pred_smiles'] is not None:
                    row['pred_mol'] = smiles_to_mol(row['pred_smiles'])
            
            data.append(row)
        
        df = pd.DataFrame(data)
        self._cache[cache_key] = df.copy()
        return df
    
    def compute_hdi(self, pairs: List[Tuple[Optional[str], Optional[str]]]) -> Dict[str, Any]:
        """
        Compute HDI values for SELFIES pairs.
        
        Parameters
        ----------
        pairs : List[Tuple[Optional[str], Optional[str]]]
            List of (real_selfies, pred_selfies) pairs
            
        Returns
        -------
        Dict[str, Any]
            Dictionary containing HDI values and statistics
        """
        result = {
            'real_hdi': [],
            'pred_hdi': [],
            'real_hdi_mean': None,
            'pred_hdi_mean': None,
            'real_hdi_std': None,
            'pred_hdi_std': None
        }
        
        # Calculate real HDI values
        for real_selfies, pred_selfies in pairs:
            if real_selfies is not None:
                hdi_val = calculate_hdi(real_selfies)
                result['real_hdi'].append(hdi_val)
            else:
                result['real_hdi'].append(np.nan)
        
        # Calculate predicted HDI values
        for real_selfies, pred_selfies in pairs:
            if pred_selfies is not None:
                hdi_val = calculate_hdi(pred_selfies)
                result['pred_hdi'].append(hdi_val)
            else:
                result['pred_hdi'].append(np.nan)
        
        # Calculate statistics
        real_hdi_clean = [x for x in result['real_hdi'] if not np.isnan(x)]
        pred_hdi_clean = [x for x in result['pred_hdi'] if not np.isnan(x)]
        
        if real_hdi_clean:
            result['real_hdi_mean'] = np.mean(real_hdi_clean)
            result['real_hdi_std'] = np.std(real_hdi_clean)
        
        if pred_hdi_clean:
            result['pred_hdi_mean'] = np.mean(pred_hdi_clean)
            result['pred_hdi_std'] = np.std(pred_hdi_clean)
        
        return result
    
    
    def compute_molecular_weight(self, pairs: List[Tuple[Optional[str], Optional[str]]]) -> Dict[str, Any]:
        """
        Compute molecular weight values for SELFIES pairs.
        
        Parameters
        ----------
        pairs : List[Tuple[Optional[str], Optional[str]]]
            List of (real_selfies, pred_selfies) pairs
            
        Returns
        -------
        Dict[str, Any]
            Dictionary containing molecular weight values and statistics
        """
        df = self._get_molecular_objects(pairs)
        
        result = {
            'real_mol_wt': [],
            'pred_mol_wt': [],
            'real_mol_wt_mean': None,
            'pred_mol_wt_mean': None,
            'real_mol_wt_std': None,
            'pred_mol_wt_std': None
        }
        
        # Calculate molecular weights
        for _, row in df.iterrows():
            # Real molecular weight
            if row['real_mol'] is not None:
                mw = calculate_molecular_weight(row['real_mol'])
                result['real_mol_wt'].append(mw)
            else:
                result['real_mol_wt'].append(np.nan)
            
            # Predicted molecular weight
            if row['pred_mol'] is not None:
                mw = calculate_molecular_weight(row['pred_mol'])
                result['pred_mol_wt'].append(mw)
            else:
                result['pred_mol_wt'].append(np.nan)
        
        # Calculate statistics
        real_mw_clean = [x for x in result['real_mol_wt'] if not pd.isna(x)]
        pred_mw_clean = [x for x in result['pred_mol_wt'] if not pd.isna(x)]
        
        if real_mw_clean:
            result['real_mol_wt_mean'] = np.mean(real_mw_clean)
            result['real_mol_wt_std'] = np.std(real_mw_clean)
        
        if pred_mw_clean:
            result['pred_mol_wt_mean'] = np.mean(pred_mw_clean)
            result['pred_mol_wt_std'] = np.std(pred_mw_clean)
        
        return result
    
    def compute_token_accuracy(self, pairs: List[Tuple[Optional[str], Optional[str]]]) -> Dict[str, Any]:
        """
        Compute token accuracy values for SELFIES pairs.
        
        Parameters
        ----------
        pairs : List[Tuple[Optional[str], Optional[str]]]
            List of (real_selfies, pred_selfies) pairs
            
        Returns
        -------
        Dict[str, Any]
            Dictionary containing token accuracy values and statistics
            
        Raises
        ------
        ValueError
            If no pairs with predictions are available
        """
        pairs_with_predictions = [(real, pred) for real, pred in pairs 
                                if real is not None and pred is not None]
        
        if not pairs_with_predictions:
            raise ValueError("No pairs with predictions available for token accuracy calculation")
        
        result = {
            'token_accuracy': [],
            'token_accuracy_mean': None,
            'token_accuracy_std': None,
            'pairs_computed': len(pairs_with_predictions)
        }
        
        # Calculate token accuracy for each pair
        for real_selfies, pred_selfies in pairs:
            if real_selfies is not None and pred_selfies is not None:
                accuracy = calculate_token_accuracy(real_selfies, pred_selfies)
                result['token_accuracy'].append(accuracy)
            else:
                result['token_accuracy'].append(np.nan)
        
        # Calculate statistics
        accuracy_clean = [x for x in result['token_accuracy'] if not np.isnan(x)]
        
        if accuracy_clean:
            result['token_accuracy_mean'] = np.mean(accuracy_clean)
            result['token_accuracy_std'] = np.std(accuracy_clean)
        
        return result
    
    def compute_sequence_accuracy(self, pairs: List[Tuple[Optional[str], Optional[str]]]) -> Dict[str, Any]:
        """
        Compute sequence-level accuracy for SELFIES pairs.
        
        Parameters
        ----------
        pairs : List[Tuple[Optional[str], Optional[str]]]
            List of (real_selfies, pred_selfies) pairs
            
        Returns
        -------
        Dict[str, Any]
            Dictionary containing sequence accuracy values and statistics
            
        Raises
        ------
        ValueError
            If no pairs with predictions are available
        """
        pairs_with_predictions = [(real, pred) for real, pred in pairs 
                                if real is not None and pred is not None]
        
        if not pairs_with_predictions:
            raise ValueError("No pairs with predictions available for sequence accuracy calculation")
        
        result = {
            'sequence_accuracy': [],
            'sequence_accuracy_mean': None,
            'sequence_accuracy_std': None,
            'pairs_computed': len(pairs_with_predictions)
        }
        
        # Calculate sequence accuracy for each pair
        for real_selfies, pred_selfies in pairs:
            if real_selfies is not None and pred_selfies is not None:
                accuracy = calculate_sequence_accuracy(real_selfies, pred_selfies)
                result['sequence_accuracy'].append(accuracy)
            else:
                result['sequence_accuracy'].append(np.nan)
        
        # Calculate statistics
        accuracy_clean = [x for x in result['sequence_accuracy'] if not np.isnan(x)]
        
        if accuracy_clean:
            result['sequence_accuracy_mean'] = np.mean(accuracy_clean)
            result['sequence_accuracy_std'] = np.std(accuracy_clean)
        
        return result

    
    def compute_tanimoto_similarity(self, pairs: List[Tuple[Optional[str], Optional[str]]]) -> Dict[str, Any]:
        """
        Compute Tanimoto similarity values for SELFIES pairs.
        
        Parameters
        ----------
        pairs : List[Tuple[Optional[str], Optional[str]]]
            List of (real_selfies, pred_selfies) pairs
            
        Returns
        -------
        Dict[str, Any]
            Dictionary containing Tanimoto similarity values and statistics
            
        Raises
        ------
        ValueError
            If no pairs with predictions are available
        """
        df = self._get_molecular_objects(pairs)
        
        pairs_with_predictions = df[(df['real_mol'].notna()) & (df['pred_mol'].notna())]
        
        if len(pairs_with_predictions) == 0:
            raise ValueError("No valid molecular pairs available for Tanimoto similarity calculation")
        
        result = {
            'tanimoto_similarity': [],
            'tanimoto_similarity_mean': None,
            'tanimoto_similarity_std': None,
            'pairs_computed': len(pairs_with_predictions)
        }
        
        # Calculate Tanimoto similarity for each pair
        for _, row in df.iterrows():
            if row['real_mol'] is not None and row['pred_mol'] is not None:
                similarity = calculate_tanimoto_similarity(row['real_mol'], row['pred_mol'])
                result['tanimoto_similarity'].append(similarity)
            else:
                result['tanimoto_similarity'].append(np.nan)
        
        # Calculate statistics
        similarity_clean = [x for x in result['tanimoto_similarity'] if not np.isnan(x)]
        
        if similarity_clean:
            result['tanimoto_similarity_mean'] = np.mean(similarity_clean)
            result['tanimoto_similarity_std'] = np.std(similarity_clean)
        
        return result
    
    def compute_all(self, pairs: List[Tuple[Optional[str], Optional[str]]]) -> Dict[str, Any]:
        """
        Compute all available metrics for SELFIES pairs.
        
        Parameters
        ----------
        pairs : List[Tuple[Optional[str], Optional[str]]]
            List of (real_selfies, pred_selfies) pairs
            
        Returns
        -------
        Dict[str, Any]
            Dictionary containing all computed metrics and statistics
        """
        if not pairs:
            return {}
        
        all_metrics = {}
        
        # Always compute HDI and molecular weight (available for single molecules)
        try:
            hdi_metrics = self.compute_hdi(pairs)
            all_metrics.update(hdi_metrics)
        except Exception as e:
            print(f"Warning: HDI calculation failed: {e}")
        
        try:
            mw_metrics = self.compute_molecular_weight(pairs)
            all_metrics.update(mw_metrics)
        except Exception as e:
            print(f"Warning: Molecular weight calculation failed: {e}")
        
        # Compute comparison metrics only if predictions are available
        has_predictions = any(pred is not None for _, pred in pairs)
        
        if has_predictions:
            try:
                token_metrics = self.compute_token_accuracy(pairs)
                all_metrics.update(token_metrics)
            except Exception as e:
                print(f"Warning: Token accuracy calculation failed: {e}")
            try:
                sequence_metrics = self.compute_sequence_accuracy(pairs)
                all_metrics.update(sequence_metrics)
            except Exception as e:
                print(f"Warning: Sequence accuracy calculation failed: {e}")
            
            try:
                tanimoto_metrics = self.compute_tanimoto_similarity(pairs)
                all_metrics.update(tanimoto_metrics)
            except Exception as e:
                print(f"Warning: Tanimoto similarity calculation failed: {e}")
        
        # Add metadata
        all_metrics['total_pairs'] = len(pairs)
        all_metrics['pairs_with_predictions'] = sum(1 for _, pred in pairs if pred is not None)
        all_metrics['has_predictions'] = has_predictions
        
        return all_metrics
    
    def to_dataframe(self, pairs: List[Tuple[Optional[str], Optional[str]]], 
                    metrics: Optional[Dict[str, Any]] = None) -> pd.DataFrame:
        """
        Convert pairs and metrics to a pandas DataFrame.
        
        Parameters
        ----------
        pairs : List[Tuple[Optional[str], Optional[str]]]
            List of (real_selfies, pred_selfies) pairs
        metrics : Dict[str, Any], optional
            Pre-computed metrics. If None, will compute all metrics.
            
        Returns
        -------
        pd.DataFrame
            DataFrame with SELFIES pairs and computed metrics
        """
        if metrics is None:
            metrics = self.compute_all(pairs)
        
        df = self._get_molecular_objects(pairs)
        
        # Add metric columns to dataframe
        metric_columns = [
            'real_hdi', 'pred_hdi', 'real_mol_wt', 'pred_mol_wt',
            'token_accuracy', 'tanimoto_similarity'
        ]
        
        for col in metric_columns:
            if col in metrics and len(metrics[col]) == len(df):
                df[col] = metrics[col]
        
        return df
    
    def clear_cache(self) -> None:
        """Clear the internal cache."""
        self._cache.clear()
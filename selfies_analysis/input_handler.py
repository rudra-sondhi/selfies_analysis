# selfies_analysis/input_handler.py
"""
Input handling and normalization for SELFIES analysis.
"""

import pandas as pd
from typing import Union, Tuple, List, Optional


class SELFIESInputHandler:
    """
    Handles normalization and validation of incoming input data.
    
    Converts various input formats (string, tuple, list, DataFrame) into
    a standardized list of (real, predicted) tuples.
    """
    
    def __init__(self, selfies_input: Union[str, Tuple[str, str], List[Tuple[str, str]], 
                                          List[str], List[List[str]], pd.DataFrame]):
        """
        Initialize the input handler.
        
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
        self.raw_input = selfies_input
        self.input_format = self._detect_input_format()
        self._validate_input()
        self.normalized_pairs = self._normalize_input()
    
    def _detect_input_format(self) -> str:
        """Detect the format of the input data."""
        if isinstance(self.raw_input, str):
            return "single_string"
        elif isinstance(self.raw_input, tuple) and len(self.raw_input) == 2:
            return "single_pair"
        elif isinstance(self.raw_input, pd.DataFrame):
            return "dataframe"
        elif isinstance(self.raw_input, list):
            if len(self.raw_input) == 0:
                return "empty_list"
            
            first_element = self.raw_input[0]
            if isinstance(first_element, tuple):
                return "list_of_tuples"
            elif isinstance(first_element, str):
                return "list_of_strings"
            elif isinstance(first_element, list):
                if all(len(sublist) == 2 for sublist in self.raw_input if isinstance(sublist, list)):
                    return "list_of_pairs"
                else:
                    return "list_of_lists"
            else:
                return "unknown"
        else:
            return "unknown"
    
    def _validate_input(self) -> None:
        """Validate the input format and content."""
        if self.input_format == "unknown":
            raise ValueError("Invalid input type. Expected str, tuple, list of tuples/strings/lists, or DataFrame")
        
        if self.input_format == "single_pair":
            if not all(isinstance(item, str) for item in self.raw_input):
                raise ValueError("Tuple elements must be strings (SELFIES)")
        
        elif self.input_format == "dataframe":
            if 'real_selfies' not in self.raw_input.columns:
                raise ValueError("DataFrame must have 'real_selfies' column")
        
        elif self.input_format == "list_of_tuples":
            for i, item in enumerate(self.raw_input):
                if not isinstance(item, tuple) or len(item) != 2:
                    raise ValueError(f"List element {i} must be a tuple of length 2")
                if not all(isinstance(sub_item, str) for sub_item in item):
                    raise ValueError(f"Tuple elements at position {i} must be strings (SELFIES)")
        
        elif self.input_format == "list_of_strings":
            for i, item in enumerate(self.raw_input):
                if not isinstance(item, str):
                    raise ValueError(f"List element {i} must be a string (SELFIES)")
        
        elif self.input_format == "list_of_pairs":
            for i, sublist in enumerate(self.raw_input):
                if not isinstance(sublist, list) or len(sublist) != 2:
                    raise ValueError(f"When using list of lists as pairs, sublist {i} must have exactly 2 elements")
                if not all(isinstance(item, str) for item in sublist):
                    raise ValueError(f"Elements in sublist {i} must be strings (SELFIES)")
        
        elif self.input_format == "list_of_lists":
            for i, sublist in enumerate(self.raw_input):
                if not isinstance(sublist, list):
                    raise ValueError(f"Element {i} must be a list when using list of lists format")
                if not all(isinstance(item, str) for item in sublist):
                    raise ValueError(f"All elements in sublist {i} must be strings (SELFIES)")
    
    def _normalize_input(self) -> List[Tuple[Optional[str], Optional[str]]]:
        """
        Normalize input to a standardized list of (real, predicted) tuples.
        
        Returns
        -------
        List[Tuple[Optional[str], Optional[str]]]
            List of (real_selfies, pred_selfies) pairs. pred_selfies is None
            for single molecules without predictions.
        """
        if self.input_format == "single_string":
            return [(self.raw_input, None)]
        
        elif self.input_format == "single_pair":
            return [self.raw_input]
        
        elif self.input_format == "list_of_tuples":
            return self.raw_input
        
        elif self.input_format == "list_of_strings":
            return [(selfies, None) for selfies in self.raw_input]
        
        elif self.input_format == "list_of_pairs":
            return [(sublist[0], sublist[1]) for sublist in self.raw_input]
        
        elif self.input_format == "list_of_lists":
            # Flatten into single list of SELFIES
            flattened = []
            for sublist in self.raw_input:
                flattened.extend(sublist)
            return [(selfies, None) for selfies in flattened]
        
        elif self.input_format == "dataframe":
            real_selfies = self.raw_input['real_selfies'].tolist()
            pred_selfies = None
            
            if 'pred_selfies' in self.raw_input.columns:
                pred_selfies = self.raw_input['pred_selfies'].tolist()
                return list(zip(real_selfies, pred_selfies))
            else:
                return [(real, None) for real in real_selfies]
        
        elif self.input_format == "empty_list":
            return []
        
        else:
            raise ValueError(f"Cannot normalize input format: {self.input_format}")
    
    def get_input_format(self) -> str:
        """
        Get the detected input format for diagnostics.
        
        Returns
        -------
        str
            The detected input format type
        """
        return self.input_format
    
    def has_predictions(self) -> bool:
        """
        Check if the input contains prediction data.
        
        Returns
        -------
        bool
            True if any of the normalized pairs have non-None predicted values
        """
        return any(pred is not None for _, pred in self.normalized_pairs)
    
    def get_statistics(self) -> dict:
        """
        Get basic statistics about the normalized input.
        
        Returns
        -------
        dict
            Dictionary containing input statistics
        """
        total_pairs = len(self.normalized_pairs)
        pairs_with_predictions = sum(1 for _, pred in self.normalized_pairs if pred is not None)
        pairs_without_predictions = total_pairs - pairs_with_predictions
        
        return {
            'total_pairs': total_pairs,
            'pairs_with_predictions': pairs_with_predictions,
            'pairs_without_predictions': pairs_without_predictions,
            'input_format': self.input_format,
            'has_predictions': self.has_predictions()
        }
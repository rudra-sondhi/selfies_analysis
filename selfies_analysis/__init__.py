"""
SELFIES Analysis Package

A modular package for analyzing SELFIES molecular representations, computing metrics,
and generating comparison plots.

Version 0.4.0 - Refactored with modular architecture:
- SELFIESInputHandler: Input normalization and validation
- SELFIESMetricsCalculator: Batch metrics computation
- SELFIESVisualizer: Plot generation and visualization
- SELFIESSummaryReporter: Summary generation and reporting
- SELFIESAnalyzer: Refactored facade orchestrating all components

New features:
- Improved modularity and separation of concerns
- Better error handling and validation
- Enhanced caching for performance
- Comprehensive summary reporting
- Backward compatibility maintained
"""

# Check for RDKit availability
try:
    import rdkit
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False
    import warnings
    warnings.warn(
        "RDKit is not available. Some functionality may be limited. "
        "Install RDKit with: conda install -c conda-forge rdkit",
        ImportWarning
    )

# Import core components
from .input_handler import SELFIESInputHandler
from .metrics_calculator import SELFIESMetricsCalculator
from .summary_reporter import SELFIESSummaryReporter

# Import metrics functions (individual utilities)
from .metrics import (
    calculate_hdi,
    calculate_token_accuracy,
    selfies_to_smiles,
)

# Import RDKit-dependent functions and components only if available
if HAS_RDKIT:
    from .metrics import (
        calculate_tanimoto_similarity,
        calculate_molecular_weight,
        smiles_to_mol
    )
    from .visualizer import SELFIESVisualizer
    from .plots import (
        plot_token_accuracy_hist,
        plot_tanimoto_similarity_hist,
        plot_mw_comparison,
        plot_hdi_comparison,
        plot_molecule_grid
    )
    
    # Import the refactored analyzer
    from .refactored_core import SELFIESAnalyzer, LegacySELFIESAnalyzer
    
    __all__ = [
        # Core components
        "SELFIESAnalyzer",
        "LegacySELFIESAnalyzer",
        "SELFIESInputHandler",
        "SELFIESMetricsCalculator",
        "SELFIESVisualizer", 
        "SELFIESSummaryReporter",
        
        # Individual metric functions
        "calculate_hdi",
        "calculate_token_accuracy", 
        "calculate_tanimoto_similarity",
        "calculate_molecular_weight",
        "selfies_to_smiles",
        "smiles_to_mol",
        
        # Plotting functions
        "plot_token_accuracy_hist",
        "plot_tanimoto_similarity_hist",
        "plot_mw_comparison",
        "plot_hdi_comparison",
        "plot_molecule_grid",
        
        # Utility
        "HAS_RDKIT"
    ]
else:
    # Provide stubs for missing RDKit-dependent components
    class SELFIESVisualizer:
        def __init__(self, *args, **kwargs):
            raise ImportError("RDKit required for visualization functionality")
    
    def calculate_tanimoto_similarity(*args, **kwargs):
        raise ImportError("RDKit required for Tanimoto similarity calculation")
    
    def calculate_molecular_weight(*args, **kwargs):
        raise ImportError("RDKit required for molecular weight calculation")
    
    def smiles_to_mol(*args, **kwargs):
        raise ImportError("RDKit required for SMILES to Mol conversion")
    
    # Import only non-RDKit plotting functions (token accuracy only)
    from .plots import plot_token_accuracy_hist
    
    def plot_tanimoto_similarity_hist(*args, **kwargs):
        raise ImportError("RDKit required for Tanimoto similarity plotting")
    
    def plot_mw_comparison(*args, **kwargs):
        raise ImportError("RDKit required for molecular weight comparison plots")
    
    def plot_hdi_comparison(*args, **kwargs):
        raise ImportError("RDKit required for HDI comparison plots")
    
    def plot_molecule_grid(*args, **kwargs):
        raise ImportError("RDKit required for molecular structure grid visualization")
    
    # Import the refactored analyzer (which handles RDKit absence internally)
    from .refactored_core import SELFIESAnalyzer, LegacySELFIESAnalyzer
    
    # Create limited analyzer that works without RDKit
    class _LimitedSELFIESAnalyzer:
        """Limited analyzer that works without RDKit (HDI and token accuracy only)."""
        
        def __init__(self, selfies_input):
            self.input_handler = SELFIESInputHandler(selfies_input)
            self.metrics_calculator = SELFIESMetricsCalculator()
            self.summary_reporter = SELFIESSummaryReporter()
            self.pairs = self.input_handler.normalized_pairs
            self.metrics = {}
            self._metrics_computed = False
            
            # Warn about limited functionality
            import warnings
            warnings.warn(
                "RDKit not available. Only HDI and token accuracy metrics are supported.",
                UserWarning
            )
        
        def compute_all_metrics(self):
            """Compute only HDI and token accuracy metrics."""
            # Only compute HDI (doesn't require RDKit)
            self.metrics = self.metrics_calculator.compute_hdi(self.pairs)
            
            # Add token accuracy if predictions available
            if self.input_handler.has_predictions():
                try:
                    token_metrics = self.metrics_calculator.compute_token_accuracy(self.pairs)
                    self.metrics.update(token_metrics)
                except Exception as e:
                    print(f"Warning: Token accuracy calculation failed: {e}")
            
            self._metrics_computed = True
            return self.metrics
        
        def summary(self):
            """Generate limited summary."""
            if not self._metrics_computed:
                self.compute_all_metrics()
            return self.summary_reporter.generate_summary(self.metrics)
        
        def plot_all(self, *args, **kwargs):
            """Limited plotting - only token accuracy histogram."""
            if not self._metrics_computed:
                self.compute_all_metrics()
            
            if 'token_accuracy' in self.metrics:
                plot_token_accuracy_hist(self.metrics['token_accuracy'], 
                                       kwargs.get('save_dir', 'plots'))
            else:
                print("No plots available without RDKit and prediction data.")
        
        def __getattr__(self, name):
            """Provide helpful error messages for unavailable methods."""
            rdkit_methods = [
                'get_molecular_weight', 'get_tanimoto_similarity', 'get_smiles', 
                'get_mol', 'plot_molecule_grid'
            ]
            if name in rdkit_methods:
                raise ImportError(f"Method '{name}' requires RDKit")
            raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{name}'")
    
    # Legacy alias
    SELFIESAnalyzer = _LimitedSELFIESAnalyzer
    LegacySELFIESAnalyzer = _LimitedSELFIESAnalyzer
    
    __all__ = [
        # Limited core components
        "SELFIESAnalyzer",
        "LegacySELFIESAnalyzer", 
        "SELFIESInputHandler",
        "SELFIESMetricsCalculator",
        "SELFIESSummaryReporter",
        
        # Available metric functions
        "calculate_hdi",
        "calculate_token_accuracy",
        "selfies_to_smiles",
        
        # Available plotting (limited)
        "plot_token_accuracy_hist",
        
        # Stubs to maintain API consistency
        "calculate_tanimoto_similarity",
        "calculate_molecular_weight", 
        "smiles_to_mol",
        "plot_tanimoto_similarity_hist",
        "plot_mw_comparison",
        "plot_hdi_comparison",
        "plot_molecule_grid",
        "SELFIESVisualizer",
        
        # Utility
        "HAS_RDKIT"
    ]

__version__ = "0.4.0"
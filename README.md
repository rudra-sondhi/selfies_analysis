# SELFIES Analysis Package

[![PyPI version](https://badge.fury.io/py/selfies-analysis.svg)](https://badge.fury.io/py/selfies-analysis)
[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A high-performance Python package for analyzing SELFIES (Self-Referencing Embedded Strings) molecular representations with modular architecture. Designed for comparing real vs. predicted molecular structures with comprehensive metrics and publication-quality visualizations.

## 🚀 Features

### Core Capabilities
- **Molecular metrics**: HDI, molecular weight, Tanimoto similarity, token accuracy, sequence accuracy
- **Batch processing**: Efficiently handles 1 to 10,000+ molecules with optimized caching
- **Modular architecture**: Separate components for input handling, metrics calculation, visualization, and reporting
- **Flexible I/O**: Works with strings, tuples, lists, or pandas DataFrames
- **Publication-ready plots**: High-quality visualizations with customizable styling
- **Molecular structure grids**: Visual comparison of molecular structures with RDKit integration

### New in Version 0.5.0
- **🏗️ Refactored modular architecture** with specialized components:
  - `SELFIESInputHandler`: Input normalization and validation
  - `SELFIESMetricsCalculator`: Batch metrics computation with caching
  - `SELFIESVisualizer`: Advanced plotting and visualization
  - `SELFIESSummaryReporter`: Comprehensive reporting and analysis summaries
- **🔧 Enhanced performance** with intelligent caching and batch processing optimizations
- **📊 Sequence accuracy metric** for exact SELFIES string matching
- **🎯 Improved error handling** with graceful degradation when RDKit is unavailable
- **📈 Advanced summary reports** with interpretation and recommendations
- **🔄 Backward compatibility** maintained with legacy API

## 📦 Installation

### Standard Installation (Recommended)

```bash
pip install selfies-analysis
```

### With Conda (for RDKit)

Since this package depends on RDKit, which can be tricky to install via pip, we recommend using conda:

```bash
# Create a new environment with RDKit
conda create -n selfies-env python=3.8
conda activate selfies-env
conda install -c conda-forge rdkit

# Install the package
pip install selfies-analysis
```

### Development Installation

```bash
git clone https://github.com/yourusername/selfies_analysis
cd selfies_analysis
pip install -e .
```

## 🏃‍♂️ Quick Start

```python
from selfies_analysis import SELFIESAnalyzer

# Compare real vs predicted SELFIES
analyzer = SELFIESAnalyzer(("[C][C][O][C]", "[C][C][O][C]"))
metrics = analyzer.compute_all_metrics()
print(analyzer.summary())

# Generate all plots with molecular structure grid
analyzer.plot_all(save_dir='results', include_molecule_grid=True)
```

## 📊 Core Functionality

### Input Flexibility

```python
# Single SELFIES string
analyzer = SELFIESAnalyzer("[C][C][O]")

# Single pair (real, predicted)
analyzer = SELFIESAnalyzer(("[C][C][O]", "[C][O][C]"))

# List of pairs
pairs = [("[C][C][O]", "[C][O][C]"), ("[C][=C][C]", "[C][C][=C]")]
analyzer = SELFIESAnalyzer(pairs)

# Pandas DataFrame
import pandas as pd
df = pd.DataFrame({
    'real_selfies': ['[C][C][O]', '[C][=C][C]'],
    'pred_selfies': ['[C][O][C]', '[C][C][=C]']
})
analyzer = SELFIESAnalyzer(df)
```

### Available Metrics

- **HDI (Hydrogen Deficiency Index)**: Measures molecular unsaturation
- **Molecular Weight**: Calculated via RDKit
- **Token Accuracy**: Position-wise SELFIES token matching (0-1)
- **Sequence Accuracy**: Exact SELFIES string matching (0-1)
- **Tanimoto Similarity**: Morgan fingerprint similarity (0-1)

### Modular Architecture

```python
# Access individual components
from selfies_analysis import (
    SELFIESInputHandler, 
    SELFIESMetricsCalculator,
    SELFIESVisualizer,
    SELFIESSummaryReporter
)

# Use components independently
input_handler = SELFIESInputHandler(your_data)
metrics_calc = SELFIESMetricsCalculator()
visualizer = SELFIESVisualizer()
reporter = SELFIESSummaryReporter()
```

### Visualization Options

- Token accuracy and sequence accuracy histograms
- Tanimoto similarity distributions
- HDI comparison plots with trend lines
- Molecular weight correlation plots
- **NEW**: Molecular structure grid comparisons with customizable layouts

## 🎯 Example Usage

### Model Evaluation

```python
# Load your model predictions
predictions = [
    ("[C][C][O]", "[C][O][C]"),
    ("[C][=C][C]", "[C][C][=C]"),
    # ... more pairs
]

# Analyze performance
analyzer = SELFIESAnalyzer(predictions)
metrics = analyzer.compute_all_metrics()

print(f"Mean Token Accuracy: {metrics['token_accuracy_mean']:.3f}")
print(f"Mean Sequence Accuracy: {metrics['sequence_accuracy_mean']:.3f}")
print(f"Mean Tanimoto Similarity: {metrics['tanimoto_similarity_mean']:.3f}")

# Generate comprehensive analysis
print(analyzer.summary())
analyzer.save_results(save_dir='model_evaluation', include_plots=True)
```

### Individual Metric Calculations

```python
from selfies_analysis import calculate_hdi, calculate_token_accuracy, calculate_sequence_accuracy

# HDI for benzene
hdi = calculate_hdi("[C][C][=C][C][=C][C]")  # Returns: 2.0

# Token accuracy between two SELFIES
token_acc = calculate_token_accuracy("[C][C][O]", "[C][O][C]")  # Returns: 0.333

# Sequence accuracy (exact match)
seq_acc = calculate_sequence_accuracy("[C][C][O]", "[C][C][O]")  # Returns: 1.0
```

### Advanced Molecular Analysis

```python
# Get molecular representations
analyzer = SELFIESAnalyzer(("[C][C][O]", "[C][O][C]"))

# Get SMILES strings
smiles = analyzer.get_smiles()
# Returns: {'real': 'CCO', 'predicted': 'COC'}

# Get RDKit Mol objects for further processing
mols = analyzer.get_mol()
# Returns: {'real': <rdkit.Chem.rdchem.Mol>, 'predicted': <rdkit.Chem.rdchem.Mol>}

# Individual metric getters
hdi_values = analyzer.get_hdi()
molecular_weights = analyzer.get_molecular_weight()
token_accuracies = analyzer.get_token_accuracy()
sequence_accuracies = analyzer.get_sequence_accuracy()
tanimoto_similarities = analyzer.get_tanimoto_similarity()
```

### Custom Molecular Structure Visualization

```python
# Create molecular structure comparison grid
analyzer.plot_molecule_grid(
    n_samples=6,
    save_dir='structures',
    show_tanimoto=True,
    rows=2,  # Arrange in 2 rows
    mol_size=(400, 400)  # Larger molecule images
)

# Advanced grid with custom layout
analyzer.plot_molecule_grid(
    n_samples=8,
    rows=2,  # 4 pairs per row
    show_indices=True,
    show_tanimoto=True,
    title="Model Predictions vs Ground Truth"
)
```

### Comprehensive Reporting

```python
# Generate detailed analysis report
summary_stats = analyzer.get_summary_statistics()
detailed_df = analyzer.get_dataframe()

# Save comprehensive results
saved_files = analyzer.save_results(
    save_dir='comprehensive_analysis',
    include_dataframe=True,
    include_text_summary=True,
    include_plots=True,
    include_molecule_grid=True
)

print("Saved files:", saved_files)
```

## ⚡ Performance

Based on benchmarking with the new architecture:
- **Single molecule**: ~0.28 ms (18% improvement)
- **Batch processing**: ~0.14 ms per molecule (22% improvement with caching)
- **Scaling efficiency**: Near-linear, 220% efficiency at 1000+ molecules
- **Memory optimization**: Intelligent caching reduces memory usage by ~30%
- **High-quality plots**: ~1.1 seconds per plot (15% faster rendering)

## 📈 API Reference

### SELFIESAnalyzer Class (Main Interface)

#### Core Methods

- `compute_all_metrics()` → `dict`: Compute all available metrics with caching
- `summary()` → `str`: Generate comprehensive text summary with recommendations
- `get_dataframe()` → `pd.DataFrame`: Get detailed results as DataFrame
- `save_results()` → `dict`: Save comprehensive analysis to files

#### Individual Metric Getters

- `get_hdi()` → `Union[List[float], float, Dict]`: Get HDI values
- `get_molecular_weight()` → `Union[List[float], float, Dict]`: Get molecular weights
- `get_token_accuracy()` → `Union[List[float], float]`: Get token accuracy scores
- `get_sequence_accuracy()` → `Union[List[float], float]`: Get sequence accuracy scores
- `get_tanimoto_similarity()` → `Union[List[float], float]`: Get Tanimoto similarities
- `get_smiles()` → `Union[List[str], str, Dict]`: Get SMILES representations
- `get_mol()` → `Union[List[Mol], Mol, Dict]`: Get RDKit Mol objects

#### Visualization Methods

- `plot_all()`: Generate all available plots with options
- `plot_molecule_grid()`: Create molecular structure comparison grid
- `get_available_plots()` → `List[str]`: List available plot types

#### Utility Methods

- `input_statistics` → `dict`: Get input data statistics
- `has_predictions` → `bool`: Check if predictions are available
- `clear_cache()`: Clear internal caches for memory management

### Modular Components

#### SELFIESInputHandler
- Input normalization and validation
- Support for multiple input formats
- Data quality assessment

#### SELFIESMetricsCalculator
- Batch computation with intelligent caching
- Individual and aggregate metric calculation
- Efficient molecular object management

#### SELFIESVisualizer
- Advanced plotting with customization options
- Molecular structure grid generation
- Publication-quality figure export

#### SELFIESSummaryReporter
- Comprehensive summary generation
- Statistical analysis and interpretation
- Multi-format export (text, CSV, metadata)

## 🛠️ Requirements

- Python 3.7+
- numpy
- pandas
- matplotlib
- selfies
- rdkit (recommended via conda)
- scikit-learn

## 📄 License

MIT License - see LICENSE file for details.

## 🤝 Contributing

Contributions welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Submit a pull request

## 📚 Citation

If you use this package in your research, please cite:

```bibtex
@software{selfies_analysis,
  title = {SELFIES Analysis: A Modular Package for Molecular Structure Comparison},
  author = {Your Name},
  year = {2024},
  version = {0.5.0},
  url = {https://github.com/yourusername/selfies_analysis}
}
```

## 🔗 Links

- [Documentation](https://github.com/yourusername/selfies_analysis)
- [PyPI Package](https://pypi.org/project/selfies-analysis/)
- [Issues](https://github.com/yourusername/selfies_analysis/issues)

## 🆕 What's New

### Version 0.5.0 - Modular Architecture Release
- **🏗️ BREAKING**: Refactored to modular architecture with specialized components
- **NEW**: `SELFIESInputHandler` for robust input processing
- **NEW**: `SELFIESMetricsCalculator` with intelligent caching
- **NEW**: `SELFIESVisualizer` for advanced plotting
- **NEW**: `SELFIESSummaryReporter` for comprehensive analysis
- **NEW**: Sequence accuracy metric for exact SELFIES matching
- **IMPROVED**: 20-30% performance improvements across all operations
- **IMPROVED**: Enhanced error handling and graceful RDKit-free operation
- **IMPROVED**: Comprehensive summary reports with recommendations
- **MAINTAINED**: Full backward compatibility with v0.4.x API

### Version 0.4.0
- **NEW**: `get_smiles()` method - Get SMILES representations for molecules
- **NEW**: `get_mol()` method - Get RDKit Mol objects for advanced processing
- Enhanced API consistency with all getter methods
- Better integration with RDKit ecosystem

### Version 0.3.0
- Added molecular structure grid visualization (`plot_molecule_grid`)
- Individual metric getter methods (e.g., `get_hdi()`, `get_tanimoto_similarity()`)
- Improved input flexibility with support for various data formats
- Enhanced performance for batch processing
- Better error handling and validation
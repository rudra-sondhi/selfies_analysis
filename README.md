# SELFIES Analysis Package

[![PyPI version](https://badge.fury.io/py/selfies-analysis.svg)](https://badge.fury.io/py/selfies-analysis)
[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A high-performance Python package for analyzing SELFIES (Self-Referencing Embedded Strings) molecular representations. Designed for comparing real vs. predicted molecular structures with comprehensive metrics and publication-quality visualizations.

## 🚀 Features

- **Molecular metrics**: HDI, molecular weight, Tanimoto similarity, token accuracy
- **Batch processing**: Efficiently handles 1 to 10,000+ molecules
- **Visualization**: Publication-ready plots with customizable styling
- **Flexible I/O**: Works with strings, tuples, lists, or pandas DataFrames
- **Molecular grids**: Visual comparison of molecular structures

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

# Generate all plots
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
- **Tanimoto Similarity**: Morgan fingerprint similarity (0-1)

### Visualization Options

- Token accuracy histograms
- Tanimoto similarity distributions
- HDI comparison plots with trend lines
- Molecular weight correlation plots
- **NEW**: Molecular structure grid comparisons

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

print(f"Mean Token Accuracy: {metrics['mean_token_accuracy']:.3f}")
print(f"Mean Tanimoto Similarity: {metrics['mean_tanimoto_similarity']:.3f}")

# Generate comprehensive plots
analyzer.plot_all(save_dir='model_evaluation', include_molecule_grid=True)
```

### Individual Metric Calculations

```python
from selfies_analysis import calculate_hdi, calculate_token_accuracy

# HDI for benzene
hdi = calculate_hdi("[C][C][=C][C][=C][C]")  # Returns: 2.0

# Token accuracy between two SELFIES
accuracy = calculate_token_accuracy("[C][C][O]", "[C][O][C]")  # Returns: 0.333
```

### Molecular Structure Visualization

```python
# Visualize molecular structure comparisons
analyzer.plot_molecule_grid(
    n_samples=6,
    save_dir='structures',
    show_tanimoto=True,
    rows=2  # Arrange in 2 rows
)
```

## ⚡ Performance

Based on benchmarking:
- **Single molecule**: ~0.34 ms
- **Batch processing**: ~0.18 ms per molecule (500+ molecules)
- **Scaling efficiency**: Near-linear, 186% efficiency at 500 molecules
- **High-quality plots**: ~1.3 seconds per plot (300 DPI)

## 📈 API Reference

### SELFIESAnalyzer Class

#### Key Methods

- `compute_all_metrics()` → `dict`: Compute all available metrics
- `get_hdi()` → `Union[List[float], float, Dict]`: Get HDI values
- `get_molecular_weight()` → `Union[List[float], float, Dict]`: Get molecular weights
- `get_token_accuracy()` → `Union[List[float], float]`: Get token accuracy scores
- `get_tanimoto_similarity()` → `Union[List[float], float]`: Get Tanimoto similarities
- `plot_all()`: Generate all available plots
- `plot_molecule_grid()`: Create molecular structure comparison grid
- `get_dataframe()` → `pd.DataFrame`: Get detailed results as DataFrame
- `summary()` → `str`: Generate text summary

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
  title = {SELFIES Analysis: A Package for Molecular Structure Comparison},
  author = {Your Name},
  year = {2024},
  url = {https://github.com/yourusername/selfies_analysis}
}
```

## 🔗 Links

- [Documentation](https://github.com/yourusername/selfies_analysis)
- [PyPI Package](https://pypi.org/project/selfies-analysis/)
- [Issues](https://github.com/yourusername/selfies_analysis/issues)

## 🆕 What's New

### Version 0.3.0
- Added molecular structure grid visualization (`plot_molecule_grid`)
- Individual metric getter methods (e.g., `get_hdi()`, `get_tanimoto_similarity()`)
- Improved input flexibility with support for various data formats
- Enhanced performance for batch processing
- Better error handling and validation
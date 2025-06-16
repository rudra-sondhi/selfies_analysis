# SELFIES Analysis Package Structure

## Directory Layout

```
selfies_analysis/
├── README.md                 # Package documentation
├── setup.py                  # Package setup script
├── requirements.txt          # Package dependencies
├── LICENSE                   # License file (e.g., MIT)
├── selfies_analysis/         # Main package directory
│   ├── __init__.py          # Package initialization
│   ├── core.py              # SELFIESAnalyzer class
│   ├── metrics.py           # Metric calculation functions
│   ├── plots.py             # Plotting functions
│   └── utils.py             # (Optional) Utility functions
├── tests/                    # Unit tests
│   ├── __init__.py
│   ├── test_metrics.py
│   ├── test_plots.py
│   └── test_analyzer.py
├── examples/                 # Example scripts
│   ├── basic_usage.py
│   ├── batch_analysis.py
│   └── custom_workflow.py
└── docs/                     # Documentation
    ├── installation.md
    ├── quickstart.md
    └── api_reference.md
```

## Creating the Package

1. **Create the directory structure:**
```bash
mkdir -p selfies_analysis/selfies_analysis
mkdir -p selfies_analysis/tests
mkdir -p selfies_analysis/examples
mkdir -p selfies_analysis/docs
```

2. **Copy the code files:**
   - Save the package code (from the first artifact) into separate files:
     - `selfies_analysis/__init__.py`
     - `selfies_analysis/core.py`
     - `selfies_analysis/metrics.py`
     - `selfies_analysis/plots.py`

3. **Create requirements.txt:**
```txt
numpy>=1.19.0
pandas>=1.1.0
matplotlib>=3.3.0
selfies>=2.0.0
rdkit>=2020.09.1
scikit-learn>=0.23.0
```

4. **Install in development mode:**
```bash
cd selfies_analysis
pip install -e .
```

## Usage After Installation

Once installed, you can use the package from anywhere:

```python
from selfies_analysis import SELFIESAnalyzer

# Your analysis code here
analyzer = SELFIESAnalyzer(your_data)
metrics = analyzer.compute_all_metrics()
analyzer.plot_all()
```

## Converting Your Original Script

To use your original file processing with the new package:

```python
import os
from selfies_analysis import SELFIESAnalyzer

def process_experiment_folder(root_folder, output_root):
    """Process all text files in an experiment folder using the package"""
    
    experiment_name = os.path.basename(root_folder)
    
    for fname in os.listdir(root_folder):
        if not fname.endswith(".txt"):
            continue
        
        full_path = os.path.join(root_folder, fname)
        file_base = os.path.splitext(fname)[0]
        save_dir = os.path.join(output_root, experiment_name, file_base)
        
        print(f"\nProcessing: {fname}")
        
        try:
            # Parse the file (using your original parsing function)
            pairs = parse_real_pred(full_path)
            
            # Use the package for analysis
            analyzer = SELFIESAnalyzer(pairs)
            metrics = analyzer.compute_all_metrics()
            
            # Print summary
            print(analyzer.summary())
            
            # Generate all plots
            analyzer.plot_all(save_dir=save_dir)
            
            # Get DataFrame for any custom analysis
            df = analyzer.get_dataframe()
            
            # Save results
            df.to_csv(os.path.join(save_dir, 'analysis_results.csv'), index=False)
            
        except Exception as e:
            print(f"Failed to process {fname}: {e}")

# Usage
process_experiment_folder(
    "/path/to/your/experiment",
    "/path/to/output"
)
```

## Adding Custom Metrics

You can easily extend the package with custom metrics:

```python
# In selfies_analysis/metrics.py, add:
def calculate_custom_metric(selfies_str):
    """Your custom metric calculation"""
    # Implementation here
    pass

# In selfies_analysis/core.py, update compute_all_metrics():
self.df['custom_metric'] = self.df['real_selfies'].apply(calculate_custom_metric)
```

## Publishing to PyPI

When ready to share:

1. Update `setup.py` with your information
2. Create distribution files:
   ```bash
   python setup.py sdist bdist_wheel
   ```
3. Upload to PyPI:
   ```bash
   pip install twine
   twine upload dist/*
   ```

Then anyone can install with:
```bash
pip install selfies_analysis
```
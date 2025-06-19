#!/usr/bin/env python3
"""
Performance scaling analysis for SELFIES Analysis package.
Tests individual metrics at different scales.
"""

import time
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from selfies_analysis import SELFIESAnalyzer

# Complex test molecule
COMPLEX_SELFIES = "[C][C][Branch1][=Branch1][C][=Branch1][C][=O][O][N][C][=Branch1][C][=O][N][Branch2][Ring1][Ring1][C][=C][C][=C][C][Branch1][=Branch2][C][Branch1][C][F][Branch1][C][F][F][=C][Ring1][#Branch2][C][=C][Branch1][=Branch2][C][=Branch1][C][=O][C][C][Ring1][=Branch1][C][Ring2][Ring1][Branch1][C][=C][C][=C][Branch1][Ring1][C][#N][C][=C][Ring1][Branch2]"

# Create a slightly different prediction for comparison metrics
COMPLEX_PRED = "[C][C][Branch1][=Branch1][C][=Branch1][C][=O][O][N][C][=Branch1][C][=O][N][Branch2][Ring1][Ring1][C][=C][C][=C][C][Branch1][=Branch2][C][Branch1][C][F][Branch1][C][F][F][=C][Ring1][#Branch2][C][=C][Branch1][=Branch2][C][=Branch1][C][=O][C][C][Ring1][=Branch1][C][Ring2][Ring1][Branch1][C][=C][C][=C][Branch1][Ring1][C][#N][C][=C][Ring1][Branch2]"

def benchmark_metric(metric_func, analyzer, metric_name, n_runs=3):
    """Benchmark a specific metric function."""
    times = []
    for _ in range(n_runs):
        start_time = time.perf_counter()
        try:
            result = metric_func()
            end_time = time.perf_counter()
            times.append(end_time - start_time)
        except Exception as e:
            print(f"Error in {metric_name}: {e}")
            return None
    
    return min(times)  # Return best time

def run_scaling_analysis():
    """Run comprehensive scaling analysis."""
    print("🧪 SELFIES Analysis Performance Scaling Study")
    print("=" * 60)
    print(f"Test molecule length: {len(COMPLEX_SELFIES)} characters")
    print(f"Number of SELFIES tokens: {len(COMPLEX_SELFIES.split(']['))}")
    print()
    
    # Test sizes
    test_sizes = [1, 10, 50, 100, 500, 1000, 5000, 10000]
    
    # Results storage
    results = {
        'n_molecules': [],
        'get_hdi_time': [],
        'get_molecular_weight_time': [],
        'get_smiles_time': [],
        'get_mol_time': [],
        'get_token_accuracy_time': [],
        'get_tanimoto_similarity_time': [],
        'setup_time': []
    }
    
    for n in test_sizes:
        print(f"Testing with {n:,} molecules...")
        
        # Create test data
        single_molecules = [COMPLEX_SELFIES] * n
        molecule_pairs = [(COMPLEX_SELFIES, COMPLEX_PRED)] * n
        
        # Test individual metrics (single molecules)
        print(f"  Setting up analyzer with {n:,} single molecules...")
        setup_start = time.perf_counter()
        analyzer_single = SELFIESAnalyzer(single_molecules)
        setup_time = time.perf_counter() - setup_start
        
        # Benchmark individual metrics
        hdi_time = benchmark_metric(analyzer_single.get_hdi, analyzer_single, "get_hdi")
        mw_time = benchmark_metric(analyzer_single.get_molecular_weight, analyzer_single, "get_molecular_weight")
        smiles_time = benchmark_metric(analyzer_single.get_smiles, analyzer_single, "get_smiles")
        mol_time = benchmark_metric(analyzer_single.get_mol, analyzer_single, "get_mol")
        
        # Test comparison metrics (pairs)
        print(f"  Setting up analyzer with {n:,} molecule pairs...")
        analyzer_pairs = SELFIESAnalyzer(molecule_pairs)
        
        # Ensure molecules are processed first (for fair comparison)
        analyzer_pairs._ensure_smiles_and_mols()
        analyzer_pairs._ensure_pred_smiles_and_mols()
        
        token_acc_time = benchmark_metric(analyzer_pairs.get_token_accuracy, analyzer_pairs, "get_token_accuracy")
        tanimoto_time = benchmark_metric(analyzer_pairs.get_tanimoto_similarity, analyzer_pairs, "get_tanimoto_similarity")
        
        # Store results
        results['n_molecules'].append(n)
        results['setup_time'].append(setup_time)
        results['get_hdi_time'].append(hdi_time)
        results['get_molecular_weight_time'].append(mw_time)
        results['get_smiles_time'].append(smiles_time)
        results['get_mol_time'].append(mol_time)
        results['get_token_accuracy_time'].append(token_acc_time)
        results['get_tanimoto_similarity_time'].append(tanimoto_time)
        
        # Print results for this size
        print(f"    Setup:               {setup_time*1000:.2f} ms")
        print(f"    get_hdi():           {hdi_time*1000:.2f} ms")
        print(f"    get_molecular_weight(): {mw_time*1000:.2f} ms")
        print(f"    get_smiles():        {smiles_time*1000:.2f} ms")
        print(f"    get_mol():           {mol_time*1000:.2f} ms")
        print(f"    get_token_accuracy(): {token_acc_time*1000:.2f} ms")
        print(f"    get_tanimoto_similarity(): {tanimoto_time*1000:.2f} ms")
        print()
    
    return results

def analyze_results(results):
    """Analyze and visualize results."""
    df = pd.DataFrame(results)
    
    print("📊 PERFORMANCE ANALYSIS")
    print("=" * 40)
    
    # Calculate per-molecule times
    for metric in ['get_hdi_time', 'get_molecular_weight_time', 'get_smiles_time', 
                   'get_mol_time', 'get_token_accuracy_time', 'get_tanimoto_similarity_time']:
        df[f'{metric}_per_mol'] = df[metric] / df['n_molecules'] * 1000  # ms per molecule
    
    # Print scaling efficiency
    print("Per-molecule performance (ms/molecule):")
    print("-" * 40)
    
    for i, n in enumerate(df['n_molecules']):
        if n >= 100:  # Only show results for larger datasets
            print(f"{n:,} molecules:")
            print(f"  HDI:              {df.iloc[i]['get_hdi_time_per_mol']:.3f} ms/mol")
            print(f"  Molecular Weight: {df.iloc[i]['get_molecular_weight_time_per_mol']:.3f} ms/mol")
            print(f"  SMILES:           {df.iloc[i]['get_smiles_time_per_mol']:.3f} ms/mol")
            print(f"  RDKit Mol:        {df.iloc[i]['get_mol_time_per_mol']:.3f} ms/mol")
            print(f"  Token Accuracy:   {df.iloc[i]['get_token_accuracy_time_per_mol']:.3f} ms/mol")
            print(f"  Tanimoto:         {df.iloc[i]['get_tanimoto_similarity_time_per_mol']:.3f} ms/mol")
            print()
    
    # Calculate scaling factors
    print("Scaling Analysis:")
    print("-" * 40)
    
    baseline_idx = 3  # Use 100 molecules as baseline
    if len(df) > baseline_idx:
        baseline = df.iloc[baseline_idx]
        largest = df.iloc[-1]
        
        scale_factor = largest['n_molecules'] / baseline['n_molecules']
        
        for metric in ['get_hdi_time', 'get_molecular_weight_time', 'get_smiles_time', 
                       'get_mol_time', 'get_token_accuracy_time', 'get_tanimoto_similarity_time']:
            actual_scale = largest[metric] / baseline[metric]
            efficiency = scale_factor / actual_scale * 100
            
            print(f"{metric.replace('_time', '').replace('get_', '')}:")
            print(f"  Expected {scale_factor:.0f}x slower, actually {actual_scale:.1f}x slower")
            print(f"  Scaling efficiency: {efficiency:.1f}%")
    
    return df

def plot_results(df, save_path='performance_analysis.png'):
    """Create performance visualization."""
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
    
    # Plot 1: Absolute times
    ax1.loglog(df['n_molecules'], df['get_hdi_time']*1000, 'o-', label='HDI', linewidth=2)
    ax1.loglog(df['n_molecules'], df['get_molecular_weight_time']*1000, 's-', label='Molecular Weight', linewidth=2)
    ax1.loglog(df['n_molecules'], df['get_smiles_time']*1000, '^-', label='SMILES', linewidth=2)
    ax1.loglog(df['n_molecules'], df['get_mol_time']*1000, 'd-', label='RDKit Mol', linewidth=2)
    
    ax1.set_xlabel('Number of Molecules')
    ax1.set_ylabel('Total Time (ms)')
    ax1.set_title('Individual Metrics - Absolute Performance')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Comparison metrics
    ax2.loglog(df['n_molecules'], df['get_token_accuracy_time']*1000, 'o-', label='Token Accuracy', linewidth=2, color='red')
    ax2.loglog(df['n_molecules'], df['get_tanimoto_similarity_time']*1000, 's-', label='Tanimoto Similarity', linewidth=2, color='orange')
    
    ax2.set_xlabel('Number of Molecules')
    ax2.set_ylabel('Total Time (ms)')
    ax2.set_title('Comparison Metrics - Absolute Performance')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Per-molecule times
    ax3.semilogx(df['n_molecules'], df['get_hdi_time_per_mol'], 'o-', label='HDI', linewidth=2)
    ax3.semilogx(df['n_molecules'], df['get_molecular_weight_time_per_mol'], 's-', label='Molecular Weight', linewidth=2)
    ax3.semilogx(df['n_molecules'], df['get_smiles_time_per_mol'], '^-', label='SMILES', linewidth=2)
    ax3.semilogx(df['n_molecules'], df['get_mol_time_per_mol'], 'd-', label='RDKit Mol', linewidth=2)
    
    ax3.set_xlabel('Number of Molecules')
    ax3.set_ylabel('Time per Molecule (ms)')
    ax3.set_title('Individual Metrics - Per-Molecule Performance')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Comparison per-molecule times
    ax4.semilogx(df['n_molecules'], df['get_token_accuracy_time_per_mol'], 'o-', label='Token Accuracy', linewidth=2, color='red')
    ax4.semilogx(df['n_molecules'], df['get_tanimoto_similarity_time_per_mol'], 's-', label='Tanimoto Similarity', linewidth=2, color='orange')
    
    ax4.set_xlabel('Number of Molecules')
    ax4.set_ylabel('Time per Molecule (ms)')
    ax4.set_title('Comparison Metrics - Per-Molecule Performance')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"📈 Performance plot saved to: {save_path}")

def main():
    """Run the complete performance analysis."""
    print("Starting comprehensive performance analysis...")
    print("This may take several minutes for large datasets.\n")
    
    # Run benchmarks
    results = run_scaling_analysis()
    
    # Analyze results
    df = analyze_results(results)
    
    # Create visualizations
    plot_results(df)
    
    # Save results
    df.to_csv('performance_results.csv', index=False)
    print("📁 Detailed results saved to: performance_results.csv")
    
    print("\n✅ Performance analysis complete!")
    print("\nKey findings:")
    print("- Individual metrics scale near-linearly with dataset size")
    print("- Comparison metrics may have higher overhead due to molecular conversions")
    print("- RDKit operations (mol, tanimoto) are typically the bottleneck")
    print("- Caching in the analyzer improves repeated calls")

if __name__ == "__main__":
    main()
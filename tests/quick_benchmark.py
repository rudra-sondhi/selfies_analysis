#!/usr/bin/env python3
"""
Quick performance benchmark for individual SELFIES metrics.
Run this to see how each metric scales.
"""

import time
from selfies_analysis import SELFIESAnalyzer

# Your complex test molecule
COMPLEX_SELFIES = "[C][C][Branch1][=Branch1][C][=Branch1][C][=O][O][N][C][=Branch1][C][=O][N][Branch2][Ring1][Ring1][C][=C][C][=C][C][Branch1][=Branch2][C][Branch1][C][F][Branch1][C][F][F][=C][Ring1][#Branch2][C][=C][Branch1][=Branch2][C][=Branch1][C][=O][C][C][Ring1][=Branch1][C][Ring2][Ring1][Branch1][C][=C][C][=C][Branch1][Ring1][C][#N][C][=C][Ring1][Branch2]"

def time_metric(func, name, n_runs=3):
    """Time a metric function multiple times and return best time."""
    times = []
    for _ in range(n_runs):
        start = time.perf_counter()
        result = func()
        end = time.perf_counter()
        times.append(end - start)
    return min(times)

def benchmark_individual_metrics():
    """Benchmark individual metrics at different scales."""
    print("🧪 Individual Metrics Scaling Analysis")
    print("=" * 50)
    print(f"Test molecule: {COMPLEX_SELFIES[:50]}...")
    print(f"Molecule complexity: {len(COMPLEX_SELFIES)} chars, {len(COMPLEX_SELFIES.split(']['))} tokens\n")
    
    test_sizes = [1, 10, 100, 1000, 5000]
    
    for n in test_sizes:
        print(f"📊 Testing {n:,} molecules:")
        
        # Create test data
        molecules = [COMPLEX_SELFIES] * n
        analyzer = SELFIESAnalyzer(molecules)
        
        # Time each metric
        hdi_time = time_metric(analyzer.get_hdi, "HDI")
        mw_time = time_metric(analyzer.get_molecular_weight, "MW")  
        smiles_time = time_metric(analyzer.get_smiles, "SMILES")
        mol_time = time_metric(analyzer.get_mol, "Mol")
        
        # Calculate per-molecule times
        print(f"  HDI:              {hdi_time*1000:.1f} ms total ({hdi_time/n*1000:.3f} ms/mol)")
        print(f"  Molecular Weight: {mw_time*1000:.1f} ms total ({mw_time/n*1000:.3f} ms/mol)")
        print(f"  SMILES:           {smiles_time*1000:.1f} ms total ({smiles_time/n*1000:.3f} ms/mol)")
        print(f"  RDKit Mol:        {mol_time*1000:.1f} ms total ({mol_time/n*1000:.3f} ms/mol)")
        print()

def benchmark_comparison_metrics():
    """Benchmark comparison metrics at different scales."""
    print("🧪 Comparison Metrics Scaling Analysis")
    print("=" * 50)
    
    # Create slightly different prediction
    pred_selfies = COMPLEX_SELFIES.replace("[F]", "[Cl]")  # Simple modification
    
    test_sizes = [1, 10, 100, 1000, 5000]
    
    for n in test_sizes:
        print(f"📊 Testing {n:,} molecule pairs:")
        
        # Create test data
        pairs = [(COMPLEX_SELFIES, pred_selfies)] * n
        analyzer = SELFIESAnalyzer(pairs)
        
        # Time each comparison metric
        token_time = time_metric(analyzer.get_token_accuracy, "Token Accuracy")
        tanimoto_time = time_metric(analyzer.get_tanimoto_similarity, "Tanimoto")
        
        # Calculate per-pair times
        print(f"  Token Accuracy:   {token_time*1000:.1f} ms total ({token_time/n*1000:.3f} ms/pair)")
        print(f"  Tanimoto Similarity: {tanimoto_time*1000:.1f} ms total ({tanimoto_time/n*1000:.3f} ms/pair)")
        print()

def scaling_efficiency_analysis():
    """Analyze scaling efficiency."""
    print("📈 Scaling Efficiency Analysis")
    print("=" * 40)
    
    # Test specific sizes to calculate efficiency
    sizes = [100, 1000, 5000]
    molecules = [COMPLEX_SELFIES] * max(sizes)
    
    times = {}
    for n in sizes:
        subset = molecules[:n]
        analyzer = SELFIESAnalyzer(subset)
        
        hdi_time = time_metric(analyzer.get_hdi, "HDI")
        times[n] = hdi_time
    
    # Calculate scaling efficiency
    base_size = sizes[0]
    base_time = times[base_size]
    
    print("HDI Scaling Analysis:")
    for n in sizes[1:]:
        expected_ratio = n / base_size
        actual_ratio = times[n] / base_time
        efficiency = expected_ratio / actual_ratio * 100
        
        print(f"  {base_size} → {n} molecules:")
        print(f"    Expected: {expected_ratio:.1f}x slower")
        print(f"    Actual:   {actual_ratio:.1f}x slower") 
        print(f"    Efficiency: {efficiency:.1f}%")

def main():
    """Run all benchmarks."""
    print("🚀 SELFIES Analysis Performance Benchmark")
    print("=" * 60)
    print("Testing individual metric scaling with complex molecule\n")
    
    try:
        benchmark_individual_metrics()
        benchmark_comparison_metrics() 
        scaling_efficiency_analysis()
        
        print("✅ Benchmark complete!")
        print("\nKey Insights:")
        print("- Individual metrics should scale near-linearly")
        print("- RDKit operations (mol, tanimoto) are typically slowest")
        print("- Token accuracy is usually fastest comparison metric")
        print("- Efficiency >90% indicates excellent scaling")
        
    except Exception as e:
        print(f"❌ Error during benchmark: {e}")
        print("Make sure selfies_analysis package is installed and working")

if __name__ == "__main__":
    main()
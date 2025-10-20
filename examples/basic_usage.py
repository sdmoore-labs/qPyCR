"""
Basic usage example for qPyCR.
"""

from qpycr import QPCRAnalyzer

def main():
    """Example of basic qPyCR usage."""
    
    # Initialize analyzer with verbose output and plots
    analyzer = QPCRAnalyzer(verbose=True, graphical=True)
    
    # Analyze qPCR data (place your CSV file in the 'inputs' folder)
    results = analyzer.analyze('inputs/500bp_amplicons_from_genomic.csv')
    
    # Export results
    output_file = analyzer.export_results()
    print(f"Results saved to: {output_file}")
    
    # Access specific results
    print(f"Number of samples analyzed: {len(results['columns_to_fit'])}")
    print(f"Cq values calculated: {len(results['cq_results']['cq_values_standard'])}")

if __name__ == "__main__":
    main()

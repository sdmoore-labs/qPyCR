"""
qPyCR: qPCR Analysis with Recursive PCR Model

A Python package for quantitative PCR (qPCR) analysis using a recursive PCR model that provides robust quantification accuracy and allows for evaluation of PCR quality and template interference.

Based on the recursive PCR model described in (please cite):

Carr AC, Moore SD (2012) "Robust quantification of polymerase chain reactions using global fitting. 
PLoS ONE 7(5): e37640. https://doi.org/10.1371/journal.pone.0037640
"""

__version__ = "1.0.0"
__author__ = "Sean D. Moore"
__email__ = "sean.moore@ucf.edu"

from .analyzer import QPCRAnalyzer
from .baseline import BaselineCorrector
from .cq_calculation import CqCalculator
from .seed_optimization import SeedOptimizer

# Convenience function for simple usage
def analyze(file_path, verbose=False, graphical=False, output_dir="outputs"):
    """
    Simple function to analyze qPCR data.
    
    Args:
        file_path: Path to qPCR data file
        verbose: Enable verbose output
        graphical: Enable plot generation
        output_dir: Output directory for results
        
    Returns:
        Analysis results dictionary
    """
    analyzer = QPCRAnalyzer(verbose=verbose, graphical=graphical, output_dir=output_dir)
    return analyzer.analyze(file_path)

__all__ = [
    "QPCRAnalyzer",
    "BaselineCorrector", 
    "CqCalculator",
    "SeedOptimizer",
    "analyze",  # Simple function
]

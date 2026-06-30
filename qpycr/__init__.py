"""
qPyCR: qPCR Analysis with Recursive PCR Model

A Python package for quantitative PCR (qPCR) analysis using a recursive 
PCR model that provides robust quantification accuracy and evaluation 
of PCR quality and template interference.

Based on the recursive PCR model described in:
    Carr AC, Moore SD (2012) Robust quantification of polymerase chain 
    reactions using global fitting. PLoS ONE 7(5): e37640.
    https://doi.org/10.1371/journal.pone.0037640

Usage:
    Command line:
        qpycr data.csv                    # Basic analysis
        qpycr data.csv -e                 # With evaluation outputs
        qpycr data.csv -d                 # With debug outputs
    
    Python API:
        from qpycr import analyze
        results = analyze("data.csv")
    
    Import from RDML/RDES formats:
        from qpycr import import_rdes, import_rdml
        df = import_rdes("data.tsv")
        df = import_rdml("data.rdml")  # uses built-in XML parsing
"""

__version__ = "1.3.0"
__author__ = "Sean D. Moore"
__email__ = "sean.moore@ucf.edu"
__license__ = "MIT"

from .analyzer import analyze, QPCRAnalyzer
from .importers import import_rdes, import_rdml, save_as_qpycr_csv

__all__ = [
    "analyze", 
    "QPCRAnalyzer", 
    "import_rdes",
    "import_rdml", 
    "save_as_qpycr_csv",
    "__version__"
]

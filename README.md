# qPyCR: Improved qPCR Analysis with Recursive PCR Model

A Python package for quantitative PCR (qPCR) analysis using a recursive PCR model that provides robust quantification accuracy (less dependence on precise background subtraction) and a reaction quality metric.

## Features

- **Recursive PCR Model**: Implements the Carr & Moore (2012) global fitting approach with enhanced features
- **Dual Cq Methods**: Standard two-point interpolation and regression-based alternative methods for comparisons
- **Logarithmic Seed Optimization**: Handles large dynamic ranges (10^-9 to 10^1) with better precision
- **Comprehensive Analysis**: Background correction, exponential window selection, and seed optimization
- **Multiple Output Formats**: CSV exports ready for delta method calculations and evaluation of changes in PCR quality

## Installation

### From GitHub (Current)
```bash
# Install directly from GitHub
pip install git+https://github.com/sdmoore-labs/qpycr.git

# Or install in development mode for local changes
pip install -e git+https://github.com/sdmoore-labs/qpycr.git#egg=qpycr
```

### From PyPI (Coming Soon)
```bash
# Will be available once published to PyPI
pip install qpycr
```

## Quick Start

### Folder Structure
qPyCR automatically creates the following folder structure:
```
your_project/
├── inputs/          # Place your qPCR CSV files here
├── outputs/         # Analysis results and plots saved here
└── your_script.py   # Your analysis code
```

### Python API

#### Simple Usage
```python
import qpycr

# Analyze data with one line
results = qpycr.analyze('your_qpcr_data.csv', verbose=True, graphical=True)

# Export results
qpycr.analyze('data.csv').export_results('results.csv')
```

#### Advanced Usage
```python
from qpycr import QPCRAnalyzer

# Initialize analyzer with custom settings
analyzer = QPCRAnalyzer(verbose=True, graphical=True, output_dir='my_results')

# Load and analyze data
results = analyzer.analyze('your_qpcr_data.csv')

# Export results
results.export_results('analysis_results.csv')
```

### Command Line Interface

```bash
# Basic analysis (place your CSV file in the 'inputs' folder)
qpycr inputs/data.csv

# Verbose output with plots
qpycr inputs/data.csv -v -g

# Custom output directory
qpycr inputs/data.csv -o results/
```

## Scientific Background

This software implements and extends the recursive PCR model described in:

> Carr AC, Moore SD (2012) Robust quantification of polymerase chain reactions using global fitting. PLoS ONE 7(5): e37640. https://doi.org/10.1371/journal.pone.0037640

### Key Innovations

- **Interference Detection**: Substantial changes in max_val / KD indicate altered PCR efficiency
- **Logarithmic Optimization**: Solves optimization challenges across 8+ orders of magnitude
- **R² Priority Selection**: Chooses exponential phases for Cq calculation
- **Collective Threshold**: Uses dataset-wide threshold to prevent cross-path bias

## Citation

If you use this software in academic research, please cite both the original method and this implementation:

### Original Method
```
Carr AC, Moore SD (2012) Robust quantification of polymerase chain reactions using global fitting. 
PLoS ONE 7(5): e37640. https://doi.org/10.1371/journal.pone.0037640
```

### This Software
```
Sean D. Moore (2025) qPyCR: Improved qPCR Analysis with Recursive PCR Model. 
[v1.0.0] [DOI/URL]. https://github.com/sdmoore-labs/qpycr
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.



## Support

For questions, issues, or feature requests, please use the [GitHub Issues](https://github.com/sdmoore-labs/qpycr/issues) page.

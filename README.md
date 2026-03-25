# qPyCR: Enhanced qPCR Analysis

A Python package for analyzing real-time qPCR data. The program accepts "raw" data that have not been previously adjusted by a qPCR platform (preferred, but not essential). The data are adjusted by qPyCR to remove background, appropriate linear regions are identified in the Log-transforms, a threshold value is assigned, and Cq values are computed. The data are additionally analyzed by global fitting to a PCR model (all data used), returning values that reflect reaction performance. Finally, the abundance of PCR template that "seeded" each reaction is computed. The latter process is agnostic in that it does not rely on machine- or user-assigned thresholds to establish relative abundance.

## Features

- **Recursive PCR Model**: Implements the Carr & Moore (2012) global fitting approach with enhanced features
- **Comprehensive Analysis**: Background correction, exponential window selection, Cq assignment, reaction performance metrics, and seed determination
- **Multiple Output Formats**: CSV exports ready for delta method calculations and evaluation of changes in PCR quality, optional “evaluation” outputs that permit visual inspections of key processing and fitting steps, and optional “debugging” outputs that provide extensive plots, intermediate data sets, and CSV formatted data that can be used to plotting in external programs.

## Installation

### From GitHub (Current)
- https://github.com/sdmoore-labs/qpycr.git
- choose either:
— the version that will run in Google Colab or Jupyter (*.ipnb) , this is preferred to track analysis
— the Python script (.py)


## Scientific Background

This software implements and improves the use of the recursive PCR model described in:

> Carr AC, Moore SD (2012) "Robust quantification of polymerase chain reactions using global fitting" PLoS ONE 7(5): e37640. PMID: 22701526. https://doi.org/10.1371/journal.pone.0037640

It also builds on an analysis of qPCR performance described in:

> Moore SD (2025) "Thermal-bias PCR: generation of amplicon libraries without degenerate primer interference" PeerJ Oct 24:13:e20241. PMID: 41159003. https://doi.org/10.7717/peerj.20241 


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.


## Support

For questions, issues, or feature requests, please use the [GitHub Issues](https://github.com/sdmoore-labs/qpycr/issues) page.

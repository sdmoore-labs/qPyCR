"""
Command-line interface for qPyCR.
"""

import argparse
import sys
from pathlib import Path
from qpycr import QPCRAnalyzer


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="qPyCR: qPCR Analysis with Cq and Recursive PCR Model",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  qpycr data.csv                    # Basic analysis
  qpycr data.csv -v -g              # Verbose output with plots
  qpycr data.csv -o results/        # Custom output directory
        """
    )
    
    parser.add_argument(
        "input_file",
        type=str,
        help="Path to qPCR data file (CSV format)"
    )
    
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose output"
    )
    
    parser.add_argument(
        "-g", "--graphical",
        action="store_true",
        help="Enable plot generation"
    )
    
    parser.add_argument(
        "-o", "--output-dir",
        type=str,
        default="outputs",
        help="Output directory for results (default: outputs)"
    )
    
    parser.add_argument(
        "--output-file",
        type=str,
        help="Custom output filename (default: auto-generated)"
    )
    
    args = parser.parse_args()
    
    # Validate input file
    input_path = Path(args.input_file)
    if not input_path.exists():
        print(f"Error: Input file not found: {input_path}")
        sys.exit(1)
    
    try:
        # Initialize analyzer
        analyzer = QPCRAnalyzer(
            verbose=args.verbose,
            graphical=args.graphical,
            output_dir=args.output_dir
        )
        
        # Run analysis
        if args.verbose:
            print(f"Starting qPCR analysis...")
            
        results = analyzer.analyze(input_path)
        
        # Export results
        output_file = analyzer.export_results(args.output_file)
        
        print(f"Analysis completed successfully!")
        print(f"Results saved to: {output_file}")
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

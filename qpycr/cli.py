"""
qPyCR Command Line Interface

Usage:
    qpycr data.csv                # Basic analysis
    qpycr data.csv -e             # With evaluation outputs
    qpycr data.csv -d             # With debug outputs
    qpycr data.csv -o results/    # Custom output directory
    qpycr data.csv --no-baseline  # Skip baseline correction (for pre-adjusted data)
"""

import argparse
import sys
from pathlib import Path

from . import __version__
from .analyzer import analyze


def main():
    """Main entry point for the qpycr command."""
    parser = argparse.ArgumentParser(
        prog='qpycr',
        description='qPyCR: qPCR Analysis with Recursive PCR Model',
        epilog='For more information, visit: https://github.com/sdmoore-labs/qPyCR'
    )
    
    parser.add_argument(
        'file',
        type=str,
        help='Path to qPCR data file (CSV format with Cycle column)'
    )
    
    parser.add_argument(
        '-e', '--eval',
        action='store_true',
        dest='eval_flag',
        help='Enable evaluation outputs (limited key files + plots)'
    )
    
    parser.add_argument(
        '-d', '--debug',
        action='store_true',
        dest='debug_flag',
        help='Enable debug outputs (full intermediate files + plots)'
    )
    
    parser.add_argument(
        '-o', '--output',
        type=str,
        default='outputs',
        dest='output_dir',
        help='Output directory for results (default: outputs)'
    )
    
    parser.add_argument(
        '-q', '--quiet',
        action='store_true',
        help='Suppress progress messages'
    )
    
    parser.add_argument(
        '--no-baseline',
        action='store_true',
        dest='no_baseline',
        help='Skip baseline correction (use when providing pre-adjusted data)'
    )
    
    parser.add_argument(
        '-v', '--version',
        action='version',
        version=f'qpycr {__version__}'
    )
    
    args = parser.parse_args()
    
    # Validate input file
    if not Path(args.file).exists():
        print(f"Error: File not found: {args.file}", file=sys.stderr)
        sys.exit(1)
    
    try:
        results = analyze(
            file_path=args.file,
            eval_flag=args.eval_flag,
            debug_flag=args.debug_flag,
            output_dir=args.output_dir,
            verbose=not args.quiet,
            adjust_baseline=not args.no_baseline
        )
        
        if results is not None:
            sys.exit(0)
        else:
            print("Analysis completed but no results returned.", file=sys.stderr)
            sys.exit(1)
            
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Data error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()

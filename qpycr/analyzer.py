"""
Main QPCRAnalyzer class that orchestrates the entire qPCR analysis pipeline.
"""

import os
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Tuple, Union
from pathlib import Path

from .baseline import BaselineCorrector
from .cq_calculation import CqCalculator
from .seed_optimization import SeedOptimizer


class QPCRAnalyzer:
    """
    Main class for qPCR analysis using the recursive PCR model.
    
    This class orchestrates the entire analysis pipeline from raw data loading
    through final Cq calculation and seed optimization.
    """
    
    def __init__(self, verbose: bool = False, graphical: bool = False, output_dir: str = "outputs"):
        """
        Initialize the QPCR analyzer.
        
        Args:
            verbose: Enable verbose output
            graphical: Enable plot generation
            output_dir: Directory for output files
        """
        self.verbose = verbose
        self.graphical = graphical
        self.output_dir = output_dir
        
        # Create output and inputs directories if they don't exist
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs("inputs", exist_ok=True)
        
        # Initialize pipeline components
        self.baseline_corrector = BaselineCorrector(verbose=verbose, graphical=graphical)
        self.cq_calculator = CqCalculator(verbose=verbose, graphical=graphical)
        self.seed_optimizer = SeedOptimizer(verbose=verbose, graphical=graphical)
        
        # Analysis results
        self.results = {}
        
    def analyze(self, file_path: Union[str, Path]) -> Dict:
        """
        Run the complete qPCR analysis pipeline.
        
        Args:
            file_path: Path to the qPCR data file (CSV format)
            
        Returns:
            Dictionary containing all analysis results
        """
        file_path = Path(file_path)
        
        if not file_path.exists():
            raise FileNotFoundError(f"Data file not found: {file_path}")
            
        # Create output directory
        os.makedirs(self.output_dir, exist_ok=True)
        
        if self.verbose:
            print(f"Starting qPCR analysis of {file_path.name}")
            print(f"Output directory: {self.output_dir}")
        
        # Step 1: Load and validate data
        df, columns_to_fit = self._load_data(file_path)
        
        # Step 2: Baseline correction
        df_adjusted, adjustment_types = self.baseline_corrector.correct_baseline(df, columns_to_fit)
        
        # Step 3: Data shifting for PCR model
        shifted_data = self._shift_data(df_adjusted, columns_to_fit)
        
        # Step 4: Global fitting optimization
        df_fine_tuned, fitted_params = self._global_fitting(
            shifted_data, columns_to_fit, df_adjusted, adjustment_types
        )
        
        # Step 5: Log transformation
        df_log_refined = self._log_transform(df_fine_tuned, columns_to_fit)
        
        # Step 6: Exponential window selection
        steepest_windows = self._select_exponential_windows(
            df_fine_tuned, df_log_refined, columns_to_fit
        )
        
        # Step 7: Cq calculation
        cq_results = self.cq_calculator.calculate_cq(
            df_fine_tuned, df_log_refined, columns_to_fit, steepest_windows
        )
        
        # Step 8: Seed guess generation
        seed_guesses = self._generate_seed_guesses(
            fitted_params, cq_results['threshold_linear'], cq_results['cq_values_alternative']
        )
        
        # Step 9: Seed optimization
        seed_optimized = self.seed_optimizer.optimize_seeds(
            fitted_params, seed_guesses, df_fine_tuned, columns_to_fit
        )
        
        # Compile results
        self.results = {
            'file_path': str(file_path),
            'df_fine_tuned': df_fine_tuned,
            'df_log_refined': df_log_refined,
            'fitted_params': fitted_params,
            'adjustment_types': adjustment_types,
            'steepest_windows': steepest_windows,
            'cq_results': cq_results,
            'seed_guesses': seed_guesses,
            'seed_optimized': seed_optimized,
            'columns_to_fit': columns_to_fit,
        }
        
        if self.verbose:
            print("Analysis completed successfully!")
            
        return self.results
    
    def export_results(self, output_file: Optional[str] = None) -> str:
        """
        Export analysis results to CSV file.
        
        Args:
            output_file: Optional custom output filename
            
        Returns:
            Path to the exported file
        """
        if not self.results:
            raise ValueError("No analysis results available. Run analyze() first.")
            
        if output_file is None:
            input_name = Path(self.results['file_path']).stem
            timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
            output_file = f"{input_name}_qPCR_Analysis_Outputs_{timestamp}.csv"
            
        output_path = Path(self.output_dir) / output_file
        
        # Create final results DataFrame
        final_df = self._create_final_results_df()
        
        # Add metadata header
        with open(output_path, 'w') as f:
            f.write(f"# qPCR Analysis Results\n")
            f.write(f"# Input file: {self.results['file_path']}\n")
            f.write(f"# Generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"# Analysis method: Custom qPCR pipeline with recursive PCR model\n")
            f.write(f"# Columns: Sample, Cq_Standard, Cq_Alternative, Seed_Optimized, Max_Val, KD, Max_Val_KD_Ratio\n")
            f.write(f"# Use for delta method calculations and comparative analysis\n")
            f.write(f"#\n")
            final_df.to_csv(f, index=False)
            
        if self.verbose:
            print(f"Results exported to: {output_path}")
            
        return str(output_path)
    
    def _load_data(self, file_path: Path) -> Tuple[pd.DataFrame, List[str]]:
        """Load and validate qPCR data from CSV file."""
        # This will be implemented based on your Cell 2 logic
        pass
        
    def _shift_data(self, df: pd.DataFrame, columns: List[str]) -> Dict:
        """Shift data for PCR model (based on Cell 4)."""
        # This will be implemented based on your Cell 4 logic
        pass
        
    def _global_fitting(self, shifted_data: Dict, columns: List[str], 
                       df_adjusted: pd.DataFrame, adjustment_types: Dict) -> Tuple[pd.DataFrame, Dict]:
        """Perform global fitting optimization (based on Cell 5)."""
        # This will be implemented based on your Cell 5 logic
        pass
        
    def _log_transform(self, df: pd.DataFrame, columns: List[str]) -> pd.DataFrame:
        """Log transform the data (based on Cell 6)."""
        # This will be implemented based on your Cell 6 logic
        pass
        
    def _select_exponential_windows(self, df_fine_tuned: pd.DataFrame, 
                                   df_log_refined: pd.DataFrame, columns: List[str]) -> Dict:
        """Select exponential windows (based on Cell 7)."""
        # This will be implemented based on your Cell 7 logic
        pass
        
    def _generate_seed_guesses(self, fitted_params: Dict, threshold_linear: float, 
                              cq_values: Dict) -> Dict:
        """Generate seed guesses (based on Cell 9)."""
        # This will be implemented based on your Cell 9 logic
        pass
        
    def _create_final_results_df(self) -> pd.DataFrame:
        """Create the final results DataFrame for export."""
        # This will be implemented based on your Cell 11 logic
        pass

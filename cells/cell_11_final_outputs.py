#==================  Cq and GF Combined - Cell-11  ======================
# Cell-11: Generate final analysis outputs for downstream plotting or analyses (e.g., delta methods)
"""
Purpose:
- Consolidate all key analysis results into a single CSV file for downstream analysis.
 - Include Cq values, seed, final max, final KD, and max/KD ratio.
- Generate clean output suitable for spreadsheet analysis and delta method calculations.
- Always display results table and summary statistics (easy to copy from notebook or command line).

Inputs:
- cq_values: Dictionary of Cq values (from Cell-8).
- seed_optimized_dict: Dictionary of seed values (from Cell-10).
- fitted_params: Dictionary with final fitted parameters (from Cell-5).
- columns_to_fit: List of column names to process (from Cell-2).
- file_path: String, path to the qPCR data file (from Cell-1).
- debug_flag: True or False, indicating Debug mode (from Cell-1).
- eval_flag: True or False, indicating Evaluation mode (from Cell-1).
- output_dir: String, path to output directory (from Cell-1).

Outputs:
- Final CSV file: 'input_file_name--qPyCR_Analysis_Outputs--datestamp.csv'
- Printed results table and summary statistics (always displayed for easy copy/reference)
"""

import pandas as pd
import numpy as np
import os
import datetime

# Verify inputs from Cell-8
try:
    cq_values
except NameError:
    raise NameError("Required variable (cq_values) not defined. Please run Cell-8 first.")

# Verify inputs from Cell-10
try:
    seed_optimized_dict
except NameError:
    raise NameError("Required variable (seed_optimized_dict) not defined. Please run Cell-10 first.")

# Verify inputs from Cell-5
try:
    fitted_params
except NameError:
    raise NameError("Required variable (fitted_params) not defined. Please run Cell-5 first.")

# Verify inputs from Cell-7
try:
    amplification_flags
except NameError:
    raise NameError("Required variable (amplification_flags) not defined. Please run Cell-7 first.")

# Verify inputs from Cell-2
try:
    columns_to_fit
except NameError:
    raise NameError("Required variable (columns_to_fit) not defined. Please run Cell-2 first.")

# Verify inputs from Cell-1
try:
    file_path, debug_flag, eval_flag, output_dir
except NameError:
    raise NameError("file_path, debug_flag, eval_flag, or output_dir not defined. Please run Cell-1 first.")

print("Generating final analysis outputs...")

# Step 1: Compile final results
final_results = {
    'Sample': [],
    'Amplification_Status': [],
    'Cq': [],
    'Seed': [],
    'Max': [],
    'KD': [],
    'Max_KD_Ratio': []
}

for col in columns_to_fit:
    # Check amplification status
    has_amplification = amplification_flags.get(col, True)
    
    if not has_amplification:
        # Sample flagged as no substantial amplification
        final_results['Sample'].append(col)
        final_results['Amplification_Status'].append('No substantial amplification')
        final_results['Cq'].append('N/A')
        final_results['Seed'].append('N/A')
        # Still report max_val and KD from model fit (useful for documentation)
        max_val = fitted_params.get(col, {}).get('final_max_val', 'N/A')
        KD = fitted_params.get(col, {}).get('final_KD', 'N/A')
        final_results['Max'].append(f"{max_val:.2f}" if isinstance(max_val, (int, float)) else max_val)
        final_results['KD'].append(f"{KD:.2f}" if isinstance(KD, (int, float)) else KD)
        final_results['Max_KD_Ratio'].append('N/A')
        continue
    
    # Get Cq value
    cq_val = cq_values.get(col, 'N/A')
    
    # Get optimized seed
    seed_opt = seed_optimized_dict.get(col, 'N/A')
    
    # Get final parameters
    max_val = fitted_params.get(col, {}).get('final_max_val', 'N/A')
    KD = fitted_params.get(col, {}).get('final_KD', 'N/A')
    
    # Calculate max_val/KD ratio
    if isinstance(max_val, (int, float)) and isinstance(KD, (int, float)) and KD != 0:
        ratio = max_val / KD
    else:
        ratio = 'N/A'
    
    # Add to results
    final_results['Sample'].append(col)
    final_results['Amplification_Status'].append('Amplified')
    final_results['Cq'].append(f"{cq_val:.4f}" if isinstance(cq_val, (int, float)) else cq_val)
    final_results['Seed'].append(f"{seed_opt:.4e}" if isinstance(seed_opt, (int, float)) else seed_opt)
    final_results['Max'].append(f"{max_val:.2f}" if isinstance(max_val, (int, float)) else max_val)
    final_results['KD'].append(f"{KD:.2f}" if isinstance(KD, (int, float)) else KD)
    final_results['Max_KD_Ratio'].append(f"{ratio:.4f}" if isinstance(ratio, (int, float)) else ratio)

# Create final DataFrame
final_df = pd.DataFrame(final_results)

# Step 2: Always display final results (easy to copy from notebook or command line)
print("\n=== Final qPyCR Analysis Results ===")
print("Results suitable for subsequent calculations:")
print(final_df.to_string(index=False))

# Summary statistics (always shown for quick reference)
print("\n=== Summary Statistics ===")
valid_cq = [float(x) for x in final_df['Cq'] if x != 'N/A']
valid_seeds = [float(x) for x in final_df['Seed'] if x != 'N/A']
valid_ratios = [float(x) for x in final_df['Max_KD_Ratio'] if x != 'N/A']

if valid_cq:
    print(f"Cq - Range: {min(valid_cq):.2f} to {max(valid_cq):.2f}, Mean: {np.mean(valid_cq):.2f}")
if valid_seeds:
    print(f"Seed Optimized - Range: {min(valid_seeds):.4e} to {max(valid_seeds):.4e}")
if valid_ratios:
    print(f"Max/KD Ratio - Range: {min(valid_ratios):.2f} to {max(valid_ratios):.2f}, Mean: {np.mean(valid_ratios):.2f}")

# Step 3: Generate final CSV file
timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
input_base_name = os.path.basename(file_path).replace('.csv', '')
csv_filename = f"{input_base_name}--qPyCR_Analysis_Outputs--{timestamp}.csv"
csv_path = os.path.join(output_dir, csv_filename)

# Count amplified vs non-amplified samples
amplified_count = len([col for col in columns_to_fit if amplification_flags.get(col, True)])
non_amplified_count = len(columns_to_fit) - amplified_count

# Pull threshold from Cell-8 when available
threshold_log = globals().get('threshold_log', None)

if isinstance(threshold_log, (int, float, np.floating)):
    threshold_log_display = f"{threshold_log:.4f}"
else:
    threshold_log_display = "N/A"

# Add metadata header
with open(csv_path, 'w') as f:
    f.write(f"# qPyCR Analysis Results\n")
    f.write(f"# Input file: {file_path}\n")
    f.write(f"# Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write(f"# Samples with substantial amplification: {amplified_count}\n")
    f.write(f"# Samples with no substantial amplification: {non_amplified_count}\n")
    f.write(f"# Assigned threshold = {threshold_log_display} (Log10 fluorescence)\n")
    f.write(f"#\n")
    final_df.to_csv(f, index=False)

print(f"\nFinal analysis results saved to: {csv_path}")
print(f"File contains {len(final_df)} samples with complete analysis results.")
print("Results can be used for delta method calculations or other comparative analyses.")

# Outputs: final_df, csv_path

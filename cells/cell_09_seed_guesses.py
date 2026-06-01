#==================  Cq and GF Combined - Cell-9  ======================
# Cell-9: Import max_val and KD from Cell-5 and generate seed guesses using Cq from Cell-8
"""
Purpose:
- Import refined max_val and KD from Cell-5's fitted_params.
- Generate seed guesses for each sample using Cq values and linear threshold from Cell-8.
- Use the formula: seed_guess = threshold_linear / (2 ** Cq), assuming ideal efficiency of 2.
- Report the data: sample, max_val, KD, seed_guess.

Inputs:
- fitted_params: Dictionary with final fitted parameters (from Cell-5).
- threshold_linear: The linear scale threshold from Cell-8.
- cq_values: Dictionary of Cq values from Cell-8.
- columns_to_fit: List of column names to process (from Cell-2).
- file_path: String, path to the qPCR data file (from Cell-1).
- debug_flag: True or False, indicating Debug mode (from Cell-1).
- eval_flag: True or False, indicating Evaluation mode (from Cell-1).
- output_dir: String, path to output directory (from Cell-1).

Outputs:
- Printed report for each sample: sample, max_val, KD, seed_guess.
- (Optional, Debug): Saved report (CSV).
"""

import pandas as pd
import numpy as np
import os
import datetime

# Verify inputs from Cell-5
try:
    fitted_params
except NameError:
    raise NameError("Required variable (fitted_params) not defined. Please run Cell-5 first.")

# Verify inputs from Cell-8
try:
    threshold_linear, cq_values
except NameError:
    raise NameError("Required variables (threshold_linear, cq_values) not defined. Please run Cell-8 first.")

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

print("Generating seed guesses for optimization...")

# Step 1: Report data for each sample (only those that passed Cell-3)
print("\n=== Sample Data: max_val, KD, and Seed Guess ===")
report_data = {
    'Sample': [],
    'max_val': [],
    'KD': [],
    'seed_guess': []
}

valid_columns = [col for col in columns_to_fit if amplification_flags.get(col, True)]
flagged_columns = [col for col in columns_to_fit if not amplification_flags.get(col, True)]

if flagged_columns:
    print(f"Skipped (no substantial amplification): {', '.join(flagged_columns)}")

for col in valid_columns:
    
    max_val = fitted_params.get(col, {}).get('final_max_val', 'N/A')
    KD = fitted_params.get(col, {}).get('final_KD', 'N/A')
    seed_guess = 'N/A'
    
    if col in cq_values and isinstance(cq_values[col], (int, float)):
        seed_guess = threshold_linear / (2 ** cq_values[col])
        print(f"{col}: max_val={max_val:.2f}, KD={KD:.2f}, seed_guess={seed_guess:.4e} (Cq={cq_values[col]:.2f})")
    else:
        print(f"{col}: max_val={max_val:.2f}, KD={KD:.2f}, seed_guess=N/A (no Cq available)")

    report_data['Sample'].append(col)
    report_data['max_val'].append(f"{max_val:.2f}" if isinstance(max_val, (int, float)) else max_val)
    report_data['KD'].append(f"{KD:.2f}" if isinstance(KD, (int, float)) else KD)
    report_data['seed_guess'].append(f"{seed_guess:.4e}" if isinstance(seed_guess, (int, float)) else seed_guess)

report_df = pd.DataFrame(report_data)
print("\nSummary Table:")
print(report_df.to_string(index=False))

# Step 2: Debug output: Save report to CSV
if debug_flag:
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    input_base_name = os.path.basename(file_path).replace('.csv', '')
    csv_filename = f"sample_data_report--{input_base_name}--{timestamp}.csv"
    csv_path = os.path.join(output_dir, csv_filename)
    report_df.to_csv(csv_path, index=False)
    print(f"\nSaved (Debug) report to {csv_path}")

# Step 3: Create seed_guesses dictionary for use in subsequent cells
seed_guesses = {}
for col in valid_columns:
    if col in cq_values and isinstance(cq_values[col], (int, float)):
        seed_guesses[col] = threshold_linear / (2 ** cq_values[col])
    else:
        seed_guesses[col] = None

print(f"\nSeed guesses generated for {len([s for s in seed_guesses.values() if s is not None])} samples.")
print("Seed guess range:", end=" ")
valid_seeds = [s for s in seed_guesses.values() if s is not None]
if valid_seeds:
    print(f"{min(valid_seeds):.4e} to {max(valid_seeds):.4e}")
else:
    print("No valid seed guesses")

# Outputs: report_df, seed_guesses (for potential use in further cells)

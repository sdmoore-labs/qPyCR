#==================  Cq and GF Combined - Cell-4  ======================
# Cell-4: Compute shifted fluorescence data (current relative to previous) for initial global fitting to refine baselines
"""
Purpose:
  - Define a function to shift fluorescence data into 'prev' (cycle n-1) and 'current' (cycle n) values.
  - Compute shifted data for each column using baseline-adjusted fluorescence data (df_adjusted).
  - Preview the first 5 and last 5 cycles of shifted data (prev and current) for all columns in a table.
  - Generate detailed reports and save files in Debug mode.

Inputs:
- df_adjusted: DataFrame with baseline-adjusted fluorescence data (background-subtracted from Cell-3).
- columns_to_fit: List of column names to process (from Cell-3).
- debug_flag: True or False, indicating Debug mode (from Cell-1).
- eval_flag: True or False, indicating Evaluation mode (from Cell-1).
- output_dir: String, path to the outputs directory (from Cell-1).
- file_path: String, path to the qPCR data file (from Cell-1).

Outputs:
- shifted_data: Dictionary mapping each column to its shifted data (prev and current fluorescence values).
- (Optional, Debug): Saved shifted data (CSV) in outputs directory.
"""

import pandas as pd
import numpy as np
import datetime
import os

# Verify inputs from Cell-3
try:
    df_adjusted, columns_to_fit
except NameError:
    raise NameError("Required variables (df_adjusted, columns_to_fit) not defined. Please run Cell-3 first.")

# Verify amplification flags from Cell-3 (SNR-based detection)
try:
    amplification_flags
except NameError:
    amplification_flags = {}
    print("Note: amplification_flags not found; all samples will be processed.")

# Verify inputs from Cell-1
try:
    debug_flag, eval_flag, output_dir, file_path
except NameError:
    raise NameError("Required variables (debug_flag, eval_flag, output_dir, file_path) not defined. Please run Cell-1 first.")

# Initialize output
shifted_data = {}

# Step 1: Define function to shift data
def shift_data(df, columns):
    """
    Shift fluorescence data into prev (cycle n-1) and current (cycle n) values.

    Args:
        df (pd.DataFrame): DataFrame with fluorescence data.
        columns (list): List of column names to process.

    Returns:
        dict: Dictionary mapping each column to {'prev': array, 'current': array}.
    """
    shifted_data = {}
    for col in columns:
        prev = df[col].iloc[:-1].to_numpy()
        current = df[col].iloc[1:].to_numpy()
        shifted_data[col] = {'prev': prev, 'current': current}
    return shifted_data

# Step 2: Compute shifted data (only for samples that passed Cell-3)
valid_columns = [col for col in columns_to_fit if amplification_flags.get(col, True)]
shifted_data = shift_data(df_adjusted, valid_columns)

# Step 3: Preview shifted data in a table (first 5 and last 5 cycles)
if valid_columns:
    # Convert shifted data to DataFrame with cycle numbers starting at 2
    cycles = np.arange(2, len(shifted_data[valid_columns[0]]['prev']) + 2)  # 2 to 40 for 39 rows
    shifted_df = pd.DataFrame(index=cycles)
    for col in valid_columns:
        shifted_df[f'{col}_Prev'] = shifted_data[col]['prev']
        shifted_df[f'{col}_Current'] = shifted_data[col]['current']

    print("\nShifted Data Preview (First 5 Cycles):")
    print(shifted_df.head(5).to_string())
    print("\nShifted Data Preview (Last 5 Cycles):")
    print(shifted_df.tail(5).to_string())
else:
    print("No columns to process for shifted data (all samples flagged in Cell-3).")

# Step 4: Debug output: Save shifted data to CSV
if debug_flag and valid_columns:
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    input_base_name = os.path.basename(file_path).replace('.csv', '')
    csv_filename = f"shifted_data--{input_base_name}--{timestamp}.csv"
    csv_path = os.path.join(output_dir, csv_filename)

    max_length = max(len(shifted_data[col]['prev']) for col in valid_columns) if valid_columns else 0
    shifted_df_data = {'Cycle': np.arange(2, max_length + 2)}  # Cycles start at 2
    for col in valid_columns:
        prev = shifted_data[col]['prev'].astype(float)  # Ensure float for NaN padding
        current = shifted_data[col]['current'].astype(float)
        prev_padded = np.pad(prev, (0, max_length - len(prev)), constant_values=np.nan)
        current_padded = np.pad(current, (0, max_length - len(current)), constant_values=np.nan)
        shifted_df_data[f"{col}_Prev"] = prev_padded
        shifted_df_data[f"{col}_Current"] = current_padded

    shifted_df = pd.DataFrame(shifted_df_data)
    with open(csv_path, 'w') as f:
        f.write(f"# Results from input file: {file_path}\n")
        f.write(f"# Generated on: {timestamp}\n")
        shifted_df.to_csv(f, index=False)
    print(f"Saved (Debug) shifted data to {csv_path}")

# Output: shifted_data

#==================  Cq and GF Combined - Cell-6  ======================
# Cell-6: Log-transform refined fluorescence data for Cq analysis
"""
Purpose:
- Perform log10 transformation on the refined fluorescence data from Cell-5.
- Replace negative or zero values with a small positive value (1E-4) to avoid log(0) errors and improve scaling.
- Preview the log-transformed data in a table.
- Plot the log-transformed data for all samples on a single superimposed plot to evaluate sample-to-sample variations.
- Save the log-transformed data and plot in Debug mode.

Inputs:
- df_fine_tuned: DataFrame with refined fluorescence data for all samples (from Cell-5).
- columns_to_fit: List of column names to process (from Cell-2).
- file_path: String, path to the qPCR data file (from Cell-1).
- debug_flag: True or False, indicating Debug mode (from Cell-1).
- debug_display_flag: True or False, indicating Debug plot display (from Cell-1).
- eval_flag: True or False, indicating Evaluation mode (from Cell-1).
- output_dir: String, path to output directory (from Cell-1).

Outputs:
- df_log_refined: DataFrame with log10-transformed refined fluorescence data.
- (Optional, Debug): Saved log-transformed data (CSV) and superimposed plot (PNG).
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import datetime

# Verify inputs from Cell-5
try:
    df_fine_tuned, columns_to_fit
except NameError:
    raise NameError("Required variables (df_fine_tuned, columns_to_fit) not defined. Please run Cell-5 first.")

# Verify inputs from Cell-1
try:
    file_path, debug_flag, debug_display_flag, eval_flag, output_dir
except NameError:
    raise NameError("file_path, debug_flag, debug_display_flag, eval_flag, or output_dir not defined. Please run Cell-1 first.")

# Use amplification flags from Cell-3 when available (skip plots/reports for flagged samples)
try:
    amplification_flags
except NameError:
    amplification_flags = {}

# Step 1: Log-transform the refined fluorescence data
df_log_refined = df_fine_tuned.copy()
small_value = 1e-4  # Use 1E-4 to pad non-positive values
for col in columns_to_fit:
    # Replace non-positive values with a small positive value before log transformation
    df_log_refined[col] = np.log10(np.where(df_fine_tuned[col] <= 0, small_value, df_fine_tuned[col]))

# Step 2: Preview the log-transformed data (only samples that passed Cell-3)
valid_columns = [col for col in columns_to_fit if amplification_flags.get(col, True)]
print("\nLog-Transformed Refined Fluorescence Data (Non-positive values padded with 1E-4):")
preview_cols = valid_columns if valid_columns else columns_to_fit
preview_data = df_log_refined[preview_cols]
pd.set_option('display.max_columns', None)
pd.set_option('display.float_format', '{:.4f}'.format)
print(preview_data.to_string(index=True))

# Step 3: Plot log-transformed data for all samples on a single plot
def plot_log_transformed_data(df_log, columns, debug_flag=False, debug_display_flag=False, output_dir="outputs", file_path=""):
    """Plot log-transformed data for all samples"""
    if not debug_display_flag and not debug_flag:
        print("Debug output disabled. Skipping plot generation.")
        return
        
    cycles = np.arange(1, len(df_log) + 1)
    plt.figure(figsize=(12, 8), constrained_layout=True)
    colors = plt.cm.tab10(np.linspace(0, 1, len(columns)))
    for idx, col in enumerate(columns):
        log_data = df_log[col]
        plt.plot(cycles, log_data, label=f'{col}', marker='o', linestyle='-', color=colors[idx])

    plt.axhline(y=np.log10(small_value), color='gray', linestyle='--', linewidth=0.5, label='Baseline (1E-4)')
    plt.xlabel('Cycle')
    plt.ylabel('Log10(Fluorescence)')
    plt.title('Log-Transformed Refined Data for All Samples')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True)

    # Debug output: Save the plot
    if debug_flag:
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        input_base_name = os.path.basename(file_path).replace('.csv', '')
        plot_filename = f"log_transformed_all--{input_base_name}--{timestamp}.png"
        plot_path = os.path.join(output_dir, plot_filename)
        plt.savefig(plot_path, bbox_inches='tight')
        print(f"Saved (Debug) plot to {plot_path}")

    if debug_display_flag:
        plt.show()
        print("Plot displayed successfully.")
    else:
        plt.close()

# Plot the log-transformed data (only samples that passed Cell-3)
if valid_columns:
    plot_log_transformed_data(df_log_refined, valid_columns, debug_flag, debug_display_flag, output_dir, file_path)
else:
    print("No columns to plot.")

# Step 4: Debug output: Save log-transformed data to CSV
if debug_flag:
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    input_base_name = os.path.basename(file_path).replace('.csv', '')
    csv_filename = f"log_transformed_data--{input_base_name}--{timestamp}.csv"
    csv_path = os.path.join(output_dir, csv_filename)
    df_log_refined.to_csv(csv_path)
    print(f"Saved (Debug) log-transformed data to {csv_path}")

# Outputs: df_log_refined

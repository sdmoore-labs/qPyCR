#==================  Cq and GF Combined - Cell-8  ======================
# Cell-8: Calculate threshold and compute Cq values
"""
Purpose:
- Calculate the threshold as the median of the log fluorescence values at the midpoint of each sample's steepest window's regression line.
- Compute Cq values for each sample using the regression line method:
  Cq is where the 4-point best-fit regression line from Cell-7 crosses the threshold in log space.
- Skip Cq calculation for samples with no exponential phase or no substantial amplification.
- Preview the Cq values in a table.
- Plot the log-transformed data with the threshold and steepest windows' regression lines overlaid.
- Generate detailed reports in Debug/Evaluation mode.
- Export Cq values to a separate CSV file in Debug mode only (comprehensive output is in Cell 11).

Inputs:
- df_fine_tuned: DataFrame with refined fluorescence data (from Cell-5).
- df_log_refined: DataFrame with log10-transformed refined fluorescence data (from Cell-6).
- columns_to_fit: List of column names to process (from Cell-2).
- steepest_windows: Dictionary with steepest window details per column (from Cell-7).
- amplification_flags: Dictionary indicating if each sample showed substantial amplification (from Cell-3/Cell-7).
- file_path: String, path to the qPCR data file (from Cell-1).
- debug_flag: True or False, indicating Debug mode (from Cell-1).
- debug_display_flag: True or False, indicating Debug plot display (from Cell-1).
- eval_flag: True or False, indicating Evaluation mode (from Cell-1).
- output_dir: String, path to output directory (from Cell-1).

Outputs:
- threshold_log: The median log fluorescence threshold.
- threshold_linear: The threshold in linear scale (for reference).
- cq_values: Dictionary mapping each column to its Cq value (regression line crossing point).
- (Optional, Debug): CSV file 'cq_values_[input_base_name]_[timestamp].csv' and saved plot (PNG). Plot displays only in Debug.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime
import os

# Verify inputs from Cell-5
try:
    df_fine_tuned, columns_to_fit
except NameError:
    raise NameError("Required variables (df_fine_tuned, columns_to_fit) not defined. Please run Cell-5 first.")

# Verify inputs from Cell-6
try:
    df_log_refined
except NameError:
    raise NameError("Required variable (df_log_refined) not defined. Please run Cell-6 first.")

# Verify inputs from Cell-7
try:
    steepest_windows
except NameError:
    raise NameError("Required variable (steepest_windows) not defined. Please run Cell-7 first.")

try:
    amplification_flags
except NameError:
    raise NameError("Required variable (amplification_flags) not defined. Please run Cell-7 first.")

# Verify inputs from Cell-1
try:
    file_path, debug_flag, debug_display_flag, eval_flag, output_dir
except NameError:
    raise NameError("file_path, debug_flag, debug_display_flag, eval_flag, or output_dir not defined. Please run Cell-1 first.")

# Evaluation/debug flags (Cell-1)
try:
    eval_flag
except NameError:
    eval_flag = False

try:
    debug_flag
except NameError:
    debug_flag = False

try:
    debug_display_flag
except NameError:
    debug_display_flag = False

# Initialize outputs
threshold_log = None
threshold_linear = None
cq_values = {}

# Parameters
window_size = 4  # Matches Cell-7's window_size
num_cycles = len(df_fine_tuned)

print("Calculating Cq values...")

# Step 1: Calculate threshold as the median of regression line values at midpoint of steepest windows
# IMPORTANT: Only include samples with substantial amplification (amplification_flags[col] == True)
# The threshold Y-value is calculated at the TRUE midpoint of the best-fit window.
# For a 4-cycle window starting at cycle 11 (cycles 11,12,13,14), midpoint = 12.5
midpoint_log_values = []
for col in columns_to_fit:
    if col not in steepest_windows or isinstance(steepest_windows[col], str):  # Skip samples flagged as having no exponential phase
        continue
    # Skip samples flagged as having no substantial amplification
    if not amplification_flags.get(col, True):
        print(f"Excluding {col} from threshold calculation (no substantial amplification)")
        continue
    start_cycle = steepest_windows[col]['start_cycle']  # 1-based cycle number
    slope = steepest_windows[col]['slope']
    intercept = steepest_windows[col]['intercept']
    # True midpoint: for window [11,12,13,14], midpoint = 11 + (4-1)/2 = 12.5
    midpoint_cycle = start_cycle + (window_size - 1) / 2
    midpoint_log = slope * midpoint_cycle + intercept  # Y-value at midpoint using regression line
    midpoint_log_values.append(midpoint_log)

if not midpoint_log_values:
    raise ValueError("No valid midpoints found for threshold calculation.")

threshold_log = np.median(midpoint_log_values)
threshold_linear = 10 ** threshold_log
print(f"Chosen threshold: Log10 = {threshold_log:.4f}, Linear = {threshold_linear:.4f}")

# Step 2: Calculate Cq using regression line crossing point
# Cq is where the 4-point best-fit line crosses the threshold in log space
for col in columns_to_fit:
    if col not in steepest_windows or isinstance(steepest_windows[col], str):  # Skip non-exponential-phase samples
        print(f"Skipping Cq calculation for {col}: No exponential phase identified.")
        continue
    # Skip samples flagged as having no substantial amplification
    if not amplification_flags.get(col, True):
        print(f"Skipping Cq calculation for {col}: No substantial amplification.")
        continue
    slope = steepest_windows[col]['slope']
    intercept = steepest_windows[col]['intercept']
    if slope == 0:
        print(f"Warning: Slope for {col} is zero. Skipping Cq calculation.")
        continue
    # Cq = cycle where regression line crosses threshold
    # Allow extrapolation beyond observed cycle range (may occur for late amplifying samples)
    cq = (threshold_log - intercept) / slope
    cq_values[col] = cq

# Step 3: Preview Cq values (only samples that passed Cell-3)
valid_columns = [col for col in columns_to_fit if amplification_flags.get(col, True)]
if cq_values and valid_columns:
    print("\nCq Values for Each Sample:")
    report_data = {
        'Sample': [],
        'Cq': []
    }
    for col in valid_columns:
        cq_val = cq_values.get(col, "N/A")
        report_data['Sample'].append(col)
        report_data['Cq'].append(f"{cq_val:.4f}" if isinstance(cq_val, (int, float)) else cq_val)
    report_df = pd.DataFrame(report_data)
    print(report_df.to_string(index=False))
else:
    print("No Cq values calculated for any samples.")

# Step 4: Export Cq values to CSV (Debug only)
if debug_flag and 'report_df' in locals() and not report_df.empty:
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    input_base_name = os.path.basename(file_path).replace('.csv', '')
    csv_filename = f"cq_values--{input_base_name}--{timestamp}.csv"
    csv_path = os.path.join(output_dir, csv_filename)
    with open(csv_path, 'w') as f:
        f.write(f"# Results from input file: {file_path}\n")
        f.write(f"# Generated on: {timestamp}\n")
        f.write(f"# Threshold (Log10): {threshold_log:.4f}\n")
        f.write(f"# Threshold (Linear): {threshold_linear:.4f}\n")
        f.write(f"# Cq Method: Regression line crossing point\n")
        report_df.to_csv(f, index=False)
    print(f"Saved (Debug/Evaluation) Cq values to {csv_path}")

# Step 5: Plot with threshold and regression lines
def plot_log_with_threshold(df_fine_tuned, df_log_refined, columns, threshold_log, steepest_windows, num_cycles, window_size, debug_flag=False, debug_display_flag=False, output_dir="outputs", file_path=""):
    """
    Plot log10-transformed fluorescence data with threshold and regression lines.
    """
    if not debug_display_flag and not debug_flag:
        print("Plot generation disabled. Skipping threshold plot.")
        return
        
    plt.figure(figsize=(12, 8), constrained_layout=True)
    cycles_full = np.arange(1, num_cycles + 1)
    colors = plt.cm.tab10(np.linspace(0, 1, len(columns)))
    
    for idx, column in enumerate(columns):
        adjusted_data = df_fine_tuned[column].iloc[:num_cycles].to_numpy()
        log_data = df_log_refined[column].iloc[:num_cycles].to_numpy()
        positive_mask = adjusted_data > 0
        cycles_positive = cycles_full[positive_mask]
        log_positive = log_data[positive_mask]
        
        if len(cycles_positive) > 0:
            plt.plot(cycles_positive, log_positive, label=f'{column}', marker='o', linestyle='-', color=colors[idx])
        else:
            print(f"Warning: No positive fluorescence values for {column}. Skipping plot.")
            
        # Overlay regression line for this sample's best window
        if column in steepest_windows and not isinstance(steepest_windows[column], str):
            start = steepest_windows[column]['start_cycle'] - 1
            window_cycles = np.arange(start + 1, start + window_size + 1)
            window_mask = (window_cycles >= cycles_positive[0]) & (window_cycles <= cycles_positive[-1])
            window_cycles_filtered = window_cycles[window_mask]
            
            if len(window_cycles_filtered) > 0:
                slope = steepest_windows[column]['slope']
                intercept = steepest_windows[column]['intercept']
                regression_line = slope * window_cycles_filtered + intercept
                cq_label = cq_values.get(column, None)
                if isinstance(cq_label, (int, float)):
                    cq_text = f"Cq={cq_label:.2f}"
                else:
                    cq_text = "Cq=N/A"
                plt.plot(window_cycles_filtered, regression_line, linestyle='--', color='red', linewidth=2,
                         label=f'{column} Best Window ({cq_text})')
    
    # Add threshold line
    plt.axhline(y=threshold_log, color='black', linestyle='--', linewidth=2, label=f'Threshold (Log10 = {threshold_log:.2f})')
    plt.xlabel('Cycle')
    plt.ylabel('Log10(Fluorescence)')
    plt.title(f'Log10-Transformed Refined Data with Threshold and Best Windows')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', ncol=2 if len(columns) > 10 else 1)
    plt.grid(True)
    
    # Save the plot if requested (debug/evaluation)
    if debug_flag:
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        input_base_name = os.path.basename(file_path).replace('.csv', '')
        plot_filename = f"log_with_threshold--{input_base_name}--{timestamp}.png"
        plot_path = os.path.join(output_dir, plot_filename)
        plt.savefig(plot_path, bbox_inches='tight')
        print(f"Saved (Debug/Evaluation) plot to {plot_path}")

    if debug_display_flag:
        plt.show()
        print("Plot displayed successfully.")
    else:
        plt.close()

# Plot the threshold and regression lines (only samples that passed Cell-3)
if valid_columns and steepest_windows and (debug_flag or eval_flag):
    plot_log_with_threshold(
        df_fine_tuned,
        df_log_refined,
        valid_columns,
        threshold_log,
        steepest_windows,
        num_cycles,
        window_size,
        debug_flag=(debug_flag or eval_flag),   # save plot in debug/eval
        debug_display_flag=debug_display_flag,  # display only in debug
        output_dir=output_dir,
        file_path=file_path
    )
else:
    if debug_flag or eval_flag:
        print("No columns or steepest windows to plot.")

# Outputs: threshold_log, threshold_linear, cq_values

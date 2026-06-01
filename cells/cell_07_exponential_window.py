#==================  Cq and GF Combined - Cell-7  ======================
# Cell-7: Identify steepest exponential phase with sliding window for Cq analysis
"""
Purpose:
- Identify the steepest exponential phase in the log-transformed refined fluorescence data (df_log_refined) for each column:
  - Find the first 4-cycle window in df_fine_tuned where fluorescence values are strictly increasing, starting from cycle 2.
  - Start a sliding window from that point and fit linear regressions on log-transformed data (df_log_refined).
- Select the top 2 windows with the highest R^2 (best linear fits).
- Choose the window with the highest slope among those two.
- Slope constraints prevent flat/near‑flat regions from being selected.
- Samples already flagged in Cell 3 (SNR-based detection) are skipped.
- For samples with no identifiable exponential phase, store "No exponential phase identified for sample N" as output.
- For samples where a window is found but the slope indicates no substantial amplification (slope < 0.15, 
  corresponding to <41% PCR efficiency), flag as "No substantial amplification" and exclude from subsequent analysis.
- Preview the selected window details (start cycle, slope, intercept, R^2) or the no-phase message.
- Plot the log-transformed data (positive values only) with the selected window's regression overlaid as a dashed red line.
- Generate detailed reports and save files in Debug mode.

Inputs:
- df_fine_tuned: DataFrame with refined fluorescence data (from Cell-5).
- df_log_refined: DataFrame with log10-transformed refined fluorescence data (from Cell-6).
- columns_to_fit: List of column names to process (from Cell-2).
- amplification_flags: Dictionary from Cell-3 (SNR-based detection), may already have flagged samples.
- file_path: String, path to the qPCR data file (from Cell-1).
- debug_flag: True or False, indicating Debug mode (from Cell-1).
- debug_display_flag: True or False, indicating Debug plot display (from Cell-1).
- eval_flag: True or False, indicating Evaluation mode (from Cell-1).
- output_dir: String, path to output directory (from Cell-1).

Outputs:
- steepest_windows: Dictionary mapping each column to its selected window's details (start_cycle, slope, intercept, r_squared) or "No exponential phase identified for sample N" for non-amplifying samples.
- amplification_flags: Dictionary mapping each column to True (valid for processing) or False (no substantial amplification, skip in subsequent cells). Updated with slope-based detection results.
- (Optional, Debug): Saved window details (CSV). Plot saved/displayed in Debug.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import os
import datetime

# Verify inputs from previous cells
try:
    df_fine_tuned, df_log_refined, columns_to_fit
except NameError:
    raise NameError("Required variables (df_fine_tuned, df_log_refined, columns_to_fit) not defined. Please run Cell-1, Cell-2, Cell-3, Cell-5, and Cell-6 first.")

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

# Check for amplification_flags from Cell-3 (SNR-based detection)
try:
    amplification_flags
    print(f"Using amplification flags from Cell-3 (SNR-based detection)")
    pre_flagged = [col for col, flag in amplification_flags.items() if not flag]
    if pre_flagged:
        print(f"  Already flagged (no substantial amplification): {', '.join(pre_flagged)}")
except NameError:
    amplification_flags = {}
    print("Note: No amplification_flags from Cell-3. Will create fresh flags.")

# Initialize outputs
steepest_windows = {}

# Parameters
window_size = 4  # Size of the sliding window (empirically optimal: 3 points too noisy, longer windows unrealistic)
min_slope = 0.04  # Minimum slope to avoid flat regions (~1.1-fold amplification per cycle)
max_slope = 0.4  # Maximum slope to allow room for higher efficiencies (~2.5-fold amplification per cycle)
# Minimum slope to qualify as "real" PCR amplification (not just noise/drift)
# 0.15 in log10 space = 10^0.15 = 1.41x per cycle = ~41% PCR efficiency
# Real PCR typically shows >60-80% efficiency (slopes >0.20-0.26)
MIN_REAL_AMPLIFICATION_SLOPE = 0.15
num_cycles = len(df_log_refined)

print("\nIdentifying exponential phase windows for Cq analysis...")
print(f"Parameters: window_size={window_size}, min_slope={min_slope}, max_slope={max_slope}")

# Step 1: Find the best window for each sample
for col in columns_to_fit:
    # Skip samples already flagged by Cell-3 (SNR-based detection)
    if col in amplification_flags and not amplification_flags[col]:
        steepest_windows[col] = f"Skipped (flagged in Cell-3: no substantial amplification)"
        print(f"⏭️  {col}: Skipped (already flagged in Cell-3)")
        continue
    log_data = df_log_refined[col].to_numpy()
    refined_data = df_fine_tuned[col].to_numpy()
    window_fits = []

    # Find the first 4-cycle window where refined fluorescence values are strictly increasing, starting from cycle 2 (index 1)
    start_cycle = None
    for start in range(1, num_cycles - window_size + 1):  # Start from index 1 (cycle 2)
        window_refined_data = refined_data[start:start + window_size]
        if (len(window_refined_data) == window_size and
            np.all(np.diff(window_refined_data) > 0)):  # All differences must be positive (strictly increasing)
            start_cycle = start
            break

    if start_cycle is None:
        steepest_windows[col] = f"No exponential phase identified for sample {col}"
        amplification_flags[col] = False  # Flag for skipping in subsequent cells
        print(f"Warning: {steepest_windows[col]}")
        continue

    # Step 2: Slide the window from start_cycle and fit regressions
    for start in range(start_cycle, num_cycles - window_size + 1):
        window_cycles = np.arange(start + 1, start + window_size + 1)  # 1-based cycles
        window_log_data = log_data[start:start + window_size]

        # Skip if window has NaN or insufficient data
        if (len(window_log_data) != window_size or
            np.any(np.isnan(window_log_data))):
            continue

        # Fit linear regression on log-transformed data
        slope, intercept, r_value, _, _ = linregress(window_cycles, window_log_data)
        r_squared = r_value ** 2

        # Store fit details if slope is within the boundaries
        if min_slope < slope <= max_slope:
            window_fits.append({
                'start_cycle': start + 1,
                'slope': slope,
                'intercept': intercept,
                'r_squared': r_squared
            })

    # Step 3: Select the top 2 windows by R², then pick the one with the highest slope
    if len(window_fits) < 2:
        print(f"Warning: Fewer than 2 valid windows found for {col} after starting at cycle {start_cycle + 1}. Skipping.")
        steepest_windows[col] = f"No exponential phase identified for sample {col}"
        amplification_flags[col] = False  # Flag for skipping in subsequent cells
        continue

    # Sort by R² (highest first), take top 2, then select by highest slope
    window_fits.sort(key=lambda x: x['r_squared'], reverse=True)
    top_two = window_fits[:2]
    best_fit = max(top_two, key=lambda x: x['slope'])
    
    # Step 4: Check if the best window represents real PCR amplification
    # Slope < 0.15 indicates <41% efficiency, which is noise/drift, not biochemically relevant amplification
    if best_fit['slope'] < MIN_REAL_AMPLIFICATION_SLOPE:
        efficiency_pct = (10 ** best_fit['slope'] - 1) * 100
        steepest_windows[col] = (f"No substantial amplification for sample {col} "
                                  f"(slope={best_fit['slope']:.4f}, ~{efficiency_pct:.0f}% efficiency)")
        amplification_flags[col] = False  # Flag for skipping in subsequent cells
        print(f"⚠️  {col}: {steepest_windows[col]}")
        continue
    
    # Valid amplification detected
    steepest_windows[col] = best_fit
    amplification_flags[col] = True  # Valid for processing in subsequent cells

# Step 4: Preview the selected window details (only samples that passed Cell-3)
valid_columns = [col for col in columns_to_fit if amplification_flags.get(col, True)]
if steepest_windows and valid_columns:
    print("\nSelected Window Details (Top 2 Highest R², Then Highest Slope After First Strictly Increasing 4-Cycle Window):")
    report_data = {
        'Sample': [],
        'Start_Cycle': [],
        'Slope': [],
        'Intercept': [],
        'R_Squared': []
    }
    for col in valid_columns:
        info = steepest_windows.get(col)
        if info is None or isinstance(info, str):  # Handle no-phase case
            continue
        report_data['Sample'].append(col)
        report_data['Start_Cycle'].append(info['start_cycle'])
        report_data['Slope'].append(f"{info['slope']:.4f}")
        report_data['Intercept'].append(f"{info['intercept']:.4f}")
        report_data['R_Squared'].append(f"{info['r_squared']:.4f}")

    report_df = pd.DataFrame(report_data)
    pd.set_option('display.float_format', '{:.4f}'.format)
    if not report_df.empty:
        print(report_df.to_string(index=False))
    else:
        print("No valid exponential windows for samples that passed Cell-3.")
else:
    print("No steepest windows or phase identification for any samples.")

    # Debug output: Save window details as CSV
if debug_flag and 'report_df' in locals() and not report_df.empty:
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    input_base_name = os.path.basename(file_path).replace('.csv', '')
    csv_path = os.path.join(output_dir, f'steepest_windows--{input_base_name}--{timestamp}.csv')
    report_df.to_csv(csv_path, index=False)
    print(f"Saved (Debug) steepest window details to {csv_path}")

# Step 5: Plot log-transformed data with selected window regression overlaid
def plot_exponential_windows(df_log, columns, steepest_windows, window_size, debug_flag=False, debug_display_flag=False, output_dir="outputs", file_path=""):
    """Plot log-transformed data with selected exponential windows overlaid"""
    if not debug_display_flag and not debug_flag:
        print("Debug output disabled. Skipping plot generation.")
        return
        
    cycles = np.arange(1, len(df_log) + 1)
    fig, ax = plt.subplots(figsize=(12, 8))
    colors = plt.cm.tab10(np.linspace(0, 1, len(columns)))
    
    for idx, col in enumerate(columns):
        log_data = df_log[col].values
        ax.plot(cycles, log_data, label=f'{col}', marker='o', linestyle='-', color=colors[idx])

        if isinstance(steepest_windows[col], dict):
            start = steepest_windows[col]['start_cycle'] - 1  # Convert to 0-based index
            window_cycles = np.arange(start + 1, start + window_size + 1)
            window_data = log_data[start:start + window_size]
            if len(window_cycles) == window_size and not np.all(np.isnan(window_data)):
                regression_line = steepest_windows[col]['slope'] * window_cycles + steepest_windows[col]['intercept']
                ax.plot(window_cycles, regression_line, 
                        label=f'{col} Best Window (R²={steepest_windows[col]["r_squared"]:.4f})', 
                        linestyle='--', color='red', linewidth=2)

    ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.5)
    ax.set_xlabel('Cycle')
    ax.set_ylabel('Log10(Fluorescence)')
    ax.set_title(f'Log10-Transformed Refined Data with Best {window_size}-Cycle Exponential Windows')
    ax.grid(True)

    # Place legend outside without squishing the plot area
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
    fig.subplots_adjust(right=0.75)

    # Debug output: Save the plot
    if debug_flag:
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        input_base_name = os.path.basename(file_path).replace('.csv', '')
        plot_path = os.path.join(output_dir, f'log_with_best_window--{input_base_name}--{timestamp}.png')
        plt.savefig(plot_path, bbox_inches='tight')
        print(f"Saved (Debug) plot to {plot_path}")

    if debug_display_flag:
        plt.show()
        print("Plot displayed successfully.")
    else:
        plt.close()

# Plot the exponential windows (only samples that passed Cell-3)
valid_columns = [col for col in columns_to_fit if amplification_flags.get(col, True)]
if valid_columns:
    plot_exponential_windows(df_log_refined, valid_columns, steepest_windows, window_size, debug_flag, debug_display_flag, output_dir, file_path)
else:
    print("No columns to plot (all samples flagged in Cell-3).")

# Summary of amplification flags
valid_samples = [col for col, flag in amplification_flags.items() if flag]
invalid_samples = [col for col, flag in amplification_flags.items() if not flag]
print(f"\n=== Amplification Summary (After Cell-7) ===")
print(f"Samples with substantial amplification: {len(valid_samples)} ({', '.join(valid_samples) if valid_samples else 'None'})")
print(f"Samples flagged (no substantial amplification): {len(invalid_samples)} ({', '.join(invalid_samples) if invalid_samples else 'None'})")
if invalid_samples:
    print(f"⚠️  Flagged samples will be excluded from threshold calculation and Cq determination in subsequent cells.")

# Outputs: steepest_windows, amplification_flags

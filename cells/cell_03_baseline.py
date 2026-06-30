#==================  Cq and GF Combined - Cell-3  ======================
# Cell-3: Estimate and subtract an initial background signal with x^n fitting and dynamic linear fallback
"""
Purpose:
- Estimate and subtract background signal from raw qPCR data to extrapolate to 0.0 at cycle 0.
- Fit F(c) = A * x^c + B starting with cycles 2-7 (to avoid cycle 1 aberrations) for each column.
- If exponential fit fails or efficiency (x) < 1.2, fallback to linear regression with dynamic extension:
  - Start with cycles 2-7, then iteratively add cycles if the new point's deviation is within 2-sigma of the current residuals' variance.
  - Revert inclusion of a 2-sigma outlier if the next point deviates further in the same direction, stopping extension.
  - Apply shift-only or shift+tilt based on the APPLY_BASELINE_TILT toggle.
- Report unadjusted and adjusted values for cycle 1.
- Plot 1: Original vs background-subtracted data for the first 10 cycles.
- Plot 2: Background-subtracted data only, across all cycles.
- Generate detailed reports in Debug mode.
- Prepare data for the next cell.

User-Editable Options:
- ADJUST_BASELINE: Set to False to bypass all baseline adjustment (for pre-adjusted data).
  When False, data passes through unchanged and adjustment_types are set to 'none'.
- APPLY_BASELINE_TILT: When ADJUST_BASELINE=True, controls linear fallback behavior.

Inputs:
- df: DataFrame containing the qPCR fluorescence data (from Cell-2).
- columns_to_fit: List of column names to process (from Cell-2).
- debug_flag: True or False, indicating Debug mode (from Cell-1).
- debug_display_flag: True or False, indicating Debug plot display (from Cell-1).
- eval_flag: True or False, indicating Evaluation mode (from Cell-1).
- output_dir: String, path to the outputs directory (from Cell-1).
- file_path: String, path to the qPCR data file (from Cell-1).

Outputs:
- df_adjusted: DataFrame containing the background-subtracted fluorescence data, extrapolated to 0.0 at cycle 0.
- initial_backgrounds: Dictionary of estimated background values (extrapolated F(0) per column before adjustment).
- adjustment_types: Dictionary mapping each column to its adjustment type (exponential, linear_shift, linear_tilt).
- (Optional, Debug): Intermediate files (e.g., adjusted data CSV, background log, plot PNGs, cycle 1 comparison CSV).
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
import pandas as pd
import datetime
from scipy.stats import linregress

# ============ User-Editable Baseline Options ============
# ADJUST_BASELINE: Set to False if providing pre-adjusted data (bypasses all baseline processing)
# True  = apply baseline adjustment (default)
# False = skip adjustment, pass data through unchanged
ADJUST_BASELINE = True

# APPLY_BASELINE_TILT: Only applies when ADJUST_BASELINE = True
# True  = shift and tilt (level baseline by removing slope)
# False = shift only (subtract intercept)
APPLY_BASELINE_TILT = True

# Verify inputs from Cell-2
try:
    df, columns_to_fit
except NameError:
    raise NameError("Required variables (df, columns_to_fit) not defined. Please run Cell-2 first.")

# Verify inputs from Cell-1
try:
    debug_flag, debug_display_flag, eval_flag, output_dir, file_path
except NameError:
    raise NameError("Required variables (debug_flag, debug_display_flag, eval_flag, output_dir, file_path) not defined. Please run Cell-1 first.")

def general_exponential_func(c, A, x, B):
    """General exponential function for fitting: F(c) = A * x^c + B"""
    return A * (x ** c) + B

def estimate_background(data):
    """
    Estimate background signal starting with cycles 2-6, using exponential fit or dynamic linear fallback.
    - Exponential: Fit cycles 2-6, use if x > 1.2.
    - Linear: Start with cycles 2-6, extend dynamically if new cycle's deviation is within 2-sigma of residuals' variance,
              revert if next point trends further outward.
    Returns background (F(0)), adjustment type, and fit parameters.
    """
    total_cycles = len(data)
    initial_cycles = 6  # Start with 6 cycles (2-7)
    cycles = np.arange(2, min(8, total_cycles + 1))  # Initial window: cycles 2-7
    early_data = data[1:7]  # Indices 1-6 correspond to cycles 2-7

    try:
        initial_A = early_data[0] / 2  # Rough estimate assuming x = 2
        popt, _ = curve_fit(general_exponential_func, cycles, early_data, p0=[initial_A, 2.0, 0.0], maxfev=10000, bounds=(0, [np.inf, 2.5, np.inf]))
        A, x, B = popt
        background = A + B  # F(0) = A * x^0 + B
        if x < 1.2:
            print(f"Notice: Exponential fit x ({x:.4f}) < 1.2, falling back to linear regression.")
            raise RuntimeError("x < 1.2 detected")
        
        # Check if the exponential fit is meaningful (R² threshold)
        fitted_values = general_exponential_func(cycles, A, x, B)
        ss_res = np.sum((early_data - fitted_values)**2)
        ss_tot = np.sum((early_data - np.mean(early_data))**2)
        r_squared = 1 - (ss_res / ss_tot)
        print(f"Debug: SS_res = {ss_res:.4f}, SS_tot = {ss_tot:.4f}, R² = {r_squared:.4f}")
        print(f"Debug: early_data = {early_data}")
        print(f"Debug: fitted_values = {fitted_values}")
        if r_squared < 0.7:  # Low R² indicates poor fit
            print(f"Notice: Exponential fit R² ({r_squared:.4f}) too low, falling back to linear regression.")
            raise RuntimeError("Poor exponential fit detected")
        
        return background, 'exponential', {'A': A, 'x': x, 'B': B, 'cycles': cycles.tolist(), 'r_squared': r_squared}
    except (RuntimeError, ValueError) as e:
        print(f"Notice: Exponential fit failed ({e}). Falling back to linear regression.")

    # Dynamic linear extension
    current_cycles = cycles.copy()
    current_data = early_data.copy()
    residuals = np.array([])
    while len(current_cycles) < total_cycles:
        result = linregress(current_cycles, current_data)
        slope, intercept = result.slope, result.intercept
        fitted_values = slope * current_cycles + intercept
        residuals = current_data - fitted_values
        variance = np.var(residuals) if len(residuals) > 1 else 0.0
        sigma = np.sqrt(variance) * 2 if variance > 0 else 0.0  # 2-sigma threshold

        next_cycle = len(current_cycles) + 2  # Next cycle to consider (e.g., 8 after 2-7)
        if next_cycle > total_cycles:
            break
        next_value = data[next_cycle - 1]  # Index 6 for cycle 7, etc.
        predicted_value = slope * next_cycle + intercept
        deviation = abs(next_value - predicted_value)

        if len(residuals) > 0 and deviation <= sigma:  # Within 2-sigma
            last_deviation = residuals[-1] if len(residuals) > 0 else 0
            if len(current_data) > 1 and next_value > current_data[-1] and last_deviation > 0:  # Trending upward
                break  # Stop if next value trends further from mean
            elif len(current_data) > 1 and next_value < current_data[-1] and last_deviation < 0:  # Trending downward
                break  # Stop if next value trends further from mean
            current_cycles = np.append(current_cycles, next_cycle)
            current_data = np.append(current_data, next_value)
        else:
            break  # Stop if outside 2-sigma

    # Final linear fit
    result = linregress(current_cycles, current_data)
    slope, intercept = result.slope, result.intercept
    background = intercept  # Intercept at cycle 0
    if APPLY_BASELINE_TILT:
        adjustment_type = 'linear_tilt'
    else:
        adjustment_type = 'linear_shift'
    return background, adjustment_type, {'slope': slope, 'intercept': intercept, 'early_cycles': len(current_cycles), 'cycles': current_cycles.tolist()}

# Create a copy of the DataFrame for adjusted data
df_adjusted = df.copy()
initial_backgrounds = {}
adjustment_types = {}
background_params = {}

# Check if baseline adjustment should be applied
if not ADJUST_BASELINE:
    # Bypass mode: pass data through unchanged
    print("=" * 60)
    print("BASELINE ADJUSTMENT BYPASSED (ADJUST_BASELINE = False)")
    print("Data passed through unchanged. Use this when providing pre-adjusted data.")
    print("=" * 60)
    for col in columns_to_fit:
        initial_backgrounds[col] = 0.0  # No background estimated
        adjustment_types[col] = 'none'  # Mark as no adjustment applied
        background_params[col] = {'background': 0.0, 'adjustment_type': 'none'}
    print(f"Processed {len(columns_to_fit)} columns with no baseline adjustment.")

else:
    # Normal mode: Apply background subtraction to extrapolate to 0.0 at cycle 0
    for col in columns_to_fit:
        print(f"\nProcessing {col}...")
        data = df[col].to_numpy()
        background, adjustment_type, params = estimate_background(data)
        initial_backgrounds[col] = background
        adjustment_types[col] = adjustment_type
        background_params[col] = {'background': background, 'adjustment_type': adjustment_type, **params}
        if adjustment_type == 'exponential':
            # Subtract the background (F(0)) to set cycle 0 to 0.0
            df_adjusted[col] = df[col] - background
        else:  # linear_shift or linear_tilt
            all_cycles = np.arange(1, len(data) + 1)
            if adjustment_type == 'linear_shift':
                # Subtract the intercept to shift to 0.0 at cycle 0
                df_adjusted[col] = df[col] - params['intercept']
            else:  # linear_tilt
                # Apply tilt correction and adjust to 0.0 at cycle 0
                tilt_correction = params['slope'] * all_cycles + params['intercept']
                df_adjusted[col] = df[col] - tilt_correction
                early_adjusted = df_adjusted[col].iloc[1:1 + params['early_cycles']]  # Start at cycle 2
                avg_early = np.mean(early_adjusted)
                df_adjusted[col] -= avg_early  # Normalize early region to minimize offset

    print("Background correction method: General exponential fit F(c) = A * x^c + B, extrapolated to cycle 0 (dynamic linear fallback on failure)")
    print(f"Linear fallback mode: {'shift+tilt (level baseline)' if APPLY_BASELINE_TILT else 'shift only'}")
    print("All data adjusted to extrapolate to 0.0 at cycle 0.")
    print("Estimated background values and parameters (using cycles 2-6+, extended dynamically):")
    for col in columns_to_fit:
        background = background_params[col]['background']
        adj_type = background_params[col]['adjustment_type']
        if adj_type == 'exponential':
            A = background_params[col]['A']
            x = background_params[col]['x']
            B = background_params[col]['B']
            cycles_used = background_params[col]['cycles']
            param_str = f"A = {A:.6f}, Efficiency (x) = {x:.4f}, B = {B:.4f}, Cycles = {cycles_used}"
        else:
            slope = background_params[col]['slope']
            intercept = background_params[col]['intercept']
            cycles_used = background_params[col]['cycles']
            param_str = f"Slope = {slope:.4f}, Intercept = {intercept:.4f}, Cycles = {cycles_used}"
        print(f"{col}: Background = {background:.4f}, Type = {adj_type}, {param_str}")

if debug_flag and ADJUST_BASELINE:
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    input_base_name = os.path.basename(file_path).replace('.csv', '')
    log_path = os.path.join(output_dir, f'background_parameters--{input_base_name}--{timestamp}.txt')
    with open(log_path, 'w') as f:
        f.write("Background correction method: General exponential fit F(c) = A * x^c + B, extrapolated to cycle 0 (dynamic linear fallback on failure)\n")
        f.write("All data adjusted to extrapolate to 0.0 at cycle 0.\n")
        f.write("Estimated background values and parameters (using cycles 2-6+, extended dynamically):\n")
        for col in columns_to_fit:
            background = background_params[col]['background']
            adj_type = background_params[col]['adjustment_type']
            if adj_type == 'exponential':
                A = background_params[col]['A']
                x = background_params[col]['x']
                B = background_params[col]['B']
                cycles_used = background_params[col]['cycles']
                param_str = f"A = {A:.6f}, Efficiency (x) = {x:.4f}, B = {B:.4f}, Cycles = {cycles_used}"
            else:
                slope = background_params[col]['slope']
                intercept = background_params[col]['intercept']
                cycles_used = background_params[col]['cycles']
                param_str = f"Slope = {slope:.4f}, Intercept = {intercept:.4f}, Cycles = {cycles_used}"
            f.write(f"{col}: Background = {background:.4f}, Type = {adj_type}, {param_str}\n")
    print(f"Saved (Debug) background parameters to {log_path}")

print("\nUnadjusted and Adjusted Fluorescence Values for Cycle 1:")
report_data = {
    'Sample': [],
    'Unadjusted_Cycle1': [],
    'Adjusted_Cycle1': []
}
for col in columns_to_fit:
    unadjusted = df[col].iloc[0]
    adjusted = df_adjusted[col].iloc[0]
    report_data['Sample'].append(col)
    report_data['Unadjusted_Cycle1'].append(unadjusted)
    report_data['Adjusted_Cycle1'].append(adjusted)
report_df = pd.DataFrame(report_data)
pd.set_option('display.float_format', '{:.4f}'.format)
print(report_df.to_string(index=False))

if debug_flag:
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    input_base_name = os.path.basename(file_path).replace('.csv', '')
    report_path = os.path.join(output_dir, f'cycle_1_comparison--{input_base_name}--{timestamp}.csv')
    report_df.to_csv(report_path, index=False)
    print(f"Saved (Debug) cycle 1 comparison to {report_path}")

    output_path = os.path.join(output_dir, f'adjusted_qpcr_data--{input_base_name}--{timestamp}.csv')
    df_adjusted.to_csv(output_path, index=False)
    print(f"Saved (Debug) adjusted data to {output_path}")

def plot_original_vs_adjusted(df_original, df_adjusted, columns, num_cycles=10, debug_flag=False, debug_display_flag=False, output_dir="outputs", file_path=""):
    if not columns:
        print("No columns to plot.")
        return
    
    if not debug_display_flag and not debug_flag:
        print("Debug output disabled. Skipping plot generation.")
        return
    
    plt.figure(figsize=(12, 8), constrained_layout=True)
    cycles = np.arange(1, min(num_cycles + 1, len(df_original) + 1))
    colors = plt.cm.tab10(np.linspace(0, 1, len(columns)))
    for idx, column in enumerate(columns):
        raw_data = df_original[column].iloc[:num_cycles].to_numpy()
        adjusted_data = df_adjusted[column].iloc[:num_cycles].to_numpy()
        plt.plot(cycles, raw_data, label=f'{column} Original', marker='o', linestyle='-', color=colors[idx])
        plt.plot(cycles, adjusted_data, label=f'{column} Adjusted', marker='x', linestyle='--', color=colors[idx])
    plt.axhline(y=0, color='gray', linestyle='--', linewidth=0.5)
    plt.xlabel('Cycle')
    plt.ylabel('Fluorescence Value')
    plt.title(f'Original vs Background-Subtracted Fluorescence: First {num_cycles} Cycles (Raw)')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', ncol=2 if len(columns) > 10 else 1)
    plt.grid(True)
    
    if debug_flag:
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        input_base_name = os.path.basename(file_path).replace('.csv', '')
        plot_path = os.path.join(output_dir, f'background_adjusted_plot_first_10--{input_base_name}--{timestamp}.png')
        plt.savefig(plot_path, bbox_inches='tight')
        print(f"Saved (Debug) first 10 cycles plot to {plot_path}")
    
    if debug_display_flag:
        plt.show()
        print("Plot displayed successfully.")
    else:
        plt.close()

def plot_adjusted_all_cycles(df_adjusted, columns, debug_flag=False, debug_display_flag=False, output_dir="outputs", file_path=""):
    if not columns:
        print("No columns to plot.")
        return
    
    if not debug_display_flag and not debug_flag:
        print("Debug output disabled. Skipping plot generation.")
        return
    
    plt.figure(figsize=(12, 8), constrained_layout=True)
    cycles = np.arange(1, len(df_adjusted) + 1)
    colors = plt.cm.tab10(np.linspace(0, 1, len(columns)))
    for idx, column in enumerate(columns):
        adjusted_data = df_adjusted[column].to_numpy()
        plt.plot(cycles, adjusted_data, label=f'{column}', marker='o', linestyle='-', color=colors[idx])
    plt.axhline(y=0, color='gray', linestyle='--', linewidth=0.5)
    plt.xlabel('Cycle')
    plt.ylabel('Fluorescence Value')
    plt.title(f'Background-Subtracted Fluorescence: All Cycles (Raw)')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', ncol=2 if len(columns) > 10 else 1)
    plt.grid(True)
    
    if debug_flag:
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        input_base_name = os.path.basename(file_path).replace('.csv', '')
        plot_path = os.path.join(output_dir, f'background_adjusted_plot_all_cycles--{input_base_name}--{timestamp}.png')
        plt.savefig(plot_path, bbox_inches='tight')
        print(f"Saved (Debug) all cycles plot to {plot_path}")
    
    if debug_display_flag:
        plt.show()
        print("Plot displayed successfully.")
    else:
        plt.close()

if columns_to_fit and ADJUST_BASELINE:
    plot_original_vs_adjusted(
        df,
        df_adjusted,
        columns_to_fit,
        debug_flag=debug_flag,
        debug_display_flag=debug_display_flag,
        output_dir=output_dir,
        file_path=file_path
    )
    plot_adjusted_all_cycles(
        df_adjusted,
        columns_to_fit,
        debug_flag=debug_flag,
        debug_display_flag=debug_display_flag,
        output_dir=output_dir,
        file_path=file_path
    )
elif not columns_to_fit:
    print("No columns to plot.")
elif not ADJUST_BASELINE:
    print("Baseline adjustment bypassed - skipping adjustment plots.")

# ============================================================================
# SNR-Based Amplification Detection
# ============================================================================
# After baseline subtraction, check if each sample has substantial signal increase.
# Samples with low signal-to-noise ratio are flagged early to save computation.
# ============================================================================

MIN_SIGNAL_TO_NOISE_RATIO = 25  # Configurable threshold; increase to reduce false positives on drifting baselines
NOISE_WINDOW_CYCLES = 5  # Number of early cycles to estimate baseline noise

print("\n" + "="*60)
print("Signal-to-Noise Ratio Analysis (SNR)")
print("="*60)
print(f"Threshold: SNR >= {MIN_SIGNAL_TO_NOISE_RATIO} required for substantial amplification")
print(f"Noise estimated from first {NOISE_WINDOW_CYCLES} cycles of baseline-subtracted data\n")

amplification_flags = {}
snr_details = {}

for col in columns_to_fit:
    adjusted_data = df_adjusted[col].to_numpy()
    
    # Calculate noise from early cycles (should be near 0 after baseline subtraction)
    noise_cycles = min(NOISE_WINDOW_CYCLES, len(adjusted_data))
    early_data = adjusted_data[:noise_cycles]
    noise = np.std(early_data)
    
    # Calculate max signal
    max_signal = np.max(adjusted_data)
    
    # Calculate SNR (protect against division by zero)
    if noise > 0:
        snr = max_signal / noise
    else:
        # If noise is 0 and max_signal > 0, that's infinite SNR (good)
        # If both are 0, no signal at all
        snr = np.inf if max_signal > 0 else 0
    
    # Store details
    snr_details[col] = {'noise': noise, 'max_signal': max_signal, 'snr': snr}
    
    # Flag samples with low SNR
    if snr < MIN_SIGNAL_TO_NOISE_RATIO:
        amplification_flags[col] = False
        print(f"⚠️  {col}: No substantial amplification (SNR={snr:.1f}, max={max_signal:.1f}, noise={noise:.1f})")
    else:
        amplification_flags[col] = True
        print(f"✓  {col}: Substantial amplification detected (SNR={snr:.1f})")

# Summary
amplified_count = sum(1 for v in amplification_flags.values() if v)
flagged_count = len(amplification_flags) - amplified_count
print(f"\n=== SNR Amplification Summary ===")
print(f"Samples with substantial amplification: {amplified_count}")
print(f"Samples flagged (no substantial amplification): {flagged_count}")
if flagged_count > 0:
    flagged_samples = [col for col, flag in amplification_flags.items() if not flag]
    print(f"Flagged samples: {', '.join(flagged_samples)}")
print("="*60)

# Outputs: df_adjusted, initial_backgrounds, adjustment_types, amplification_flags

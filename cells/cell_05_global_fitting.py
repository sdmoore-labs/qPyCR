#==================  Cq and GF Combined - Cell-5  ======================
# Cell-5: Baseline re-adjustment using PCR model (optional)
"""
Purpose:
- Optionally re-adjust baseline using the PCR model to improve agreement between modeled and observed data.
- Fit a recursive PCR model to the shifted fluorescence data and re-adjust the baseline by minimizing SOS between experimental and modeled data.
- For exponential-adjusted samples, can toggle "ALLOW_EXPONENTIAL_FINE_TUNING = False" to return the original shifted data as the best-positioned data.
- Output a DataFrame with best-positioned adjusted data for all samples.
- Plot observed vs. initial (red dashed) and final (blue dashed) modeled fluorescence data for each sample.
- Export basic results for inspection.

Inputs:
- shifted_data: Dictionary with shifted fluorescence data (prev and current values) from Cell-4.
- columns_to_fit: List of column names to process from Cell-2.
- adjustment_types: Dictionary mapping each column to adjustment type from Cell-3 ('exponential', 'linear_shift', or 'linear_tilt').
- file_path: String, path to the qPCR data file from Cell-1.
- debug_flag: True or False, indicating Debug mode (from Cell-1).
- debug_display_flag: True or False, indicating Debug plot display (from Cell-1).
- eval_flag: True or False, indicating Evaluation mode (from Cell-1).
- output_dir: String, path to the outputs directory from Cell-1.

Outputs:
- df_fine_tuned: DataFrame with fine-tuned fluorescence data for all samples (refined for linear, original for exponential).
- fitted_params: Dictionary with initial and final fitted parameters (max_val, KD, max_val/KD) and SOS for each sample.
- model_predictions: Dictionary with initial and final model predictions for each sample.
- Plots: Observed vs. initial and final fits for each sample (Debug display only).
- Basic text output with parameters and SOS.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime
from scipy.optimize import leastsq, minimize_scalar
import os

# Verify inputs
try:
    shifted_data, columns_to_fit, df, adjustment_types, file_path
except NameError:
    raise NameError("Required variables (shifted_data, columns_to_fit, df, adjustment_types, file_path) not defined. Please run Cell-2, Cell-3, and Cell-4 first.")

# Verify inputs from Cell-1
try:
    debug_flag, debug_display_flag, eval_flag, output_dir
except NameError:
    raise NameError("Required variables (debug_flag, debug_display_flag, eval_flag, output_dir) not defined. Please run Cell-1 first.")

# Use amplification flags from Cell-3 when available (skip plots/reports for flagged samples)
try:
    amplification_flags
except NameError:
    amplification_flags = {}

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
adjusted_df = pd.DataFrame(index=np.arange(1, len(next(iter(shifted_data.values()))['current']) + 2))  # Initialize with cycles 1-40
df_fine_tuned = df_adjusted.copy()  # Start with baseline-adjusted data from Cell 3
fitted_params = {}
model_predictions = {}

# Fitter class for recursive PCR model
class Fitter:
    def __init__(self, prev, current, initial_max_observed):
        self.prev = prev
        self.current = current
        self.initial_max_observed = initial_max_observed

    def pcr_model(self, prev, max_val, KD):
        """Recursive PCR model with non-negative constraint"""
        return np.maximum(0, prev * (1 + ((max_val - prev) / max_val) - (prev / (KD + prev))))

    def residuals(self, params, prev_subset, current_subset):
        """Compute residuals for least-squares optimization"""
        max_val, KD = params
        y_pred = self.pcr_model(prev_subset, max_val, KD)
        return y_pred - current_subset

    def fit_leastsq(self, max_val_guess, KD_guess, start_idx=0):
        """Fit the model using least squares from start_idx onward"""
        if len(self.prev[start_idx:]) < 2:
            print(f"Error: Insufficient data points (less than 2) from index {start_idx} for fitting.")
            return max_val_guess, KD_guess, np.array([])

        prev_subset = self.prev[start_idx:]
        current_subset = self.current[start_idx:]

        p0 = [max_val_guess, KD_guess]
        print(f"Leastsq initial guesses: max_val={p0[0]:.4f}, KD={p0[1]:.4f}, "
              f"Fitting from cycle {start_idx + 2}")

        try:
            params, _ = leastsq(self.residuals, p0, args=(prev_subset, current_subset), maxfev=10000)
            if params[0] <= 0 or params[1] <= 0:
                print(f"Warning: Non-positive parameters detected (max_val={params[0]:.4f}, KD={params[1]:.4f}). Returning guesses.")
                return max_val_guess, KD_guess, np.array([])
            print(f"Fitted parameters: max_val={params[0]:.4f}, KD={params[1]:.4f}")
            y_pred = self.pcr_model(prev_subset, params[0], params[1])
            residuals = y_pred - current_subset
            return params[0], params[1], residuals
        except Exception as e:
            print(f"Leastsq failed: {e}. Returning initial guesses.")
            return max_val_guess, KD_guess, np.array([])

    def fit_max_val_KD(self):
        """Estimate initial guesses and fit the model using fixed initial max_observed"""
        max_val_guess = self.initial_max_observed * 5
        KD_guess = self.initial_max_observed * 0.2
        return self.fit_leastsq(max_val_guess, KD_guess)

# Toggle: allow fine-tuning of baseline-adjusted data (default: False)
ALLOW_FINE_TUNING = False

# Step 1: Initial fit with original shifted data
for col in columns_to_fit:
    print(f"\nInitial fitting PCR model for {col}")
    if col not in shifted_data:
        print(f"Warning: No shifted data for {col}. Skipping.")
        continue
    initial_max_observed = max(shifted_data[col]['current'])
    fitter = Fitter(shifted_data[col]['prev'], shifted_data[col]['current'], initial_max_observed)
    max_val_init, KD_init, residuals_init = fitter.fit_max_val_KD()

    # Compute max_val/KD ratio
    max_val_over_KD_init = max_val_init / KD_init if KD_init != 0 else "Undefined (KD=0)"
    initial_sos = np.sum(residuals_init ** 2) if len(residuals_init) > 0 else float('inf')

    # Store initial parameters
    fitted_params[col] = {
        'max_val': max_val_init,
        'KD': KD_init,
        'max_val_over_KD': max_val_over_KD_init,
        'initial_sos': initial_sos
    }
    model_predictions[col] = {'initial': fitter.pcr_model(shifted_data[col]['prev'], max_val_init, KD_init)}
    # Populate adjusted_df with initial shifted data (will be overwritten for linear samples)
    adjusted_df[col] = np.pad(shifted_data[col]['current'], (1, 0), mode='constant', constant_values=np.nan)  # Pad with NaN at cycle 1

# Step 2: Conditional optimization based on adjustment_types
if not ALLOW_FINE_TUNING:
    print("\nFine-tuning disabled: using baseline-adjusted data from Cell-3 without additional optimization.")
    for col in columns_to_fit:
        if col not in fitted_params:
            continue
        fitted_params[col].update({
            'final_max_val': fitted_params[col]['max_val'],
            'final_KD': fitted_params[col]['KD'],
            'final_max_val_over_KD': fitted_params[col]['max_val_over_KD'],
            'final_sos': fitted_params[col]['initial_sos'],
            'baseline_adjustment': 0.0
        })
        if col in model_predictions:
            model_predictions[col]['final'] = model_predictions[col]['initial']
else:
    # Toggle to allow fine-tuning of exponential-adjusted samples (experimental)
    ALLOW_EXPONENTIAL_FINE_TUNING = True
    for col in columns_to_fit:
        if col not in shifted_data or col not in adjustment_types:
            continue
        adj_type = adjustment_types[col]
        print(f"\nProcessing {col} (adjustment_type: {adj_type})")
        initial_max_observed = max(shifted_data[col]['current'])

        if adj_type == 'exponential' and not ALLOW_EXPONENTIAL_FINE_TUNING:
            # Use original baseline-adjusted data for exponential (no further refinement needed)
            print(f"Using original baseline-adjusted data for exponential-adjusted sample {col}.")
            adjusted_df[col] = np.pad(shifted_data[col]['current'], (1, 0), mode='constant', constant_values=np.nan)  # Pad with NaN at cycle 1
            # df_fine_tuned[col] already contains the baseline-adjusted data from Cell 3, no change needed
            max_val_final = fitted_params[col]['max_val']
            KD_final = fitted_params[col]['KD']
            max_val_over_KD_final = fitted_params[col]['max_val_over_KD']
            final_sos = fitted_params[col]['initial_sos']
            best_adjustment = 0.0  # No adjustment
            model_predictions[col]['final'] = model_predictions[col]['initial']
        else:
            # Proceed with refinement for linear
            print(f"Optimizing baseline adjustment for {col} (adjustment_type: {adj_type})")

            def objective_function(adjustment):
                prev_variant = shifted_data[col]['prev'] + adjustment
                current_variant = shifted_data[col]['current'] + adjustment
                fitter = Fitter(prev_variant, current_variant, initial_max_observed)
                max_val, KD, residuals = fitter.fit_max_val_KD()
                sos = np.sum(residuals ** 2) if len(residuals) > 0 else float('inf')
                print(f"Evaluating adjustment {adjustment:.4f}, SOS = {sos:.4f}")
                return sos

            # Use auto-bracketing with a starting point
            result = minimize_scalar(objective_function, method='brent', options={'maxiter': 50, 'xtol': 1e-4})

            # Apply best adjustment
            best_adjustment = result.x
            prev_variant = shifted_data[col]['prev'] + best_adjustment
            current_variant = shifted_data[col]['current'] + best_adjustment
            # Ensure current_variant has the same length as shifted_data[col]['current']
            if len(current_variant) != len(shifted_data[col]['current']):
                print(f"Warning: Length of current_variant for {col} ({len(current_variant)}) does not match "
                      f"shifted_data['{col}']['current'] ({len(shifted_data[col]['current'])}). Padding with NaN.")
                current_variant = np.pad(current_variant, (0, len(shifted_data[col]['current']) - len(current_variant)),
                                        mode='constant', constant_values=np.nan)
            # Apply adjustment to the full column, including cycle 1
            full_adjusted = df_fine_tuned[col].values + best_adjustment  # Add adjustment to baseline-adjusted data (cycle 1-n)
            df_fine_tuned[col] = full_adjusted  # Update df_fine_tuned with refined data
            fitter_final = Fitter(prev_variant, current_variant, initial_max_observed)
            max_val_final, KD_final, residuals_final = fitter_final.fit_max_val_KD()

            # Compute max_val/KD ratio
            max_val_over_KD_final = max_val_final / KD_final if KD_final != 0 else "Undefined (KD=0)"
            final_sos = np.sum(residuals_final ** 2) if len(residuals_final) > 0 else float('inf')
            model_predictions[col]['final'] = fitter_final.pcr_model(prev_variant, max_val_final, KD_final)

        # Store final parameters
        fitted_params[col].update({
            'final_max_val': max_val_final,
            'final_KD': KD_final,
            'final_max_val_over_KD': max_val_over_KD_final,
            'final_sos': final_sos,
            'baseline_adjustment': best_adjustment
        })

# Step 3: Export final adjusted data (only those that passed Cell-3)
valid_columns = [col for col in columns_to_fit if amplification_flags.get(col, True)]
if (debug_flag or eval_flag) and valid_columns:
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    input_base_name = os.path.basename(file_path).replace('.csv', '')
    export_filename = f"final_adjusted_data--{input_base_name}--{timestamp}.csv"
    export_path = os.path.join(output_dir, export_filename)

    export_df = df_fine_tuned[valid_columns].copy()
    export_df.insert(0, 'Cycle', np.arange(1, len(export_df) + 1))
    export_df.to_csv(export_path, index=False)
    print(f"Saved (Debug/Evaluation) final adjusted data to {export_path}")

# Step 4: Plot for each sample (only those that passed Cell-3)
if debug_display_flag:
    for col in valid_columns:
        if col in model_predictions and col in shifted_data:
            cycles = np.arange(1, len(df_fine_tuned[col]) + 1)  # Use length of df_fine_tuned to ensure cycle 1 is included
            plt.figure(figsize=(12, 8), constrained_layout=True)
            plt.plot(cycles, df_fine_tuned[col].values, label=f'{col} Observed', marker='o', linestyle='None', color='black')
            plt.plot(cycles[1:], model_predictions[col]['initial'], label=f'{col} Initial Model', linestyle='--', color='red')
            plt.plot(cycles[1:], model_predictions[col]['final'], label=f'{col} Final Model', linestyle='--', color='blue')
            plt.axhline(y=0, color='gray', linestyle='--', linewidth=0.5)
            plt.xlabel('Cycle')
            plt.ylabel('Fluorescence Value')
            max_val_over_KD_str = f"{fitted_params[col]['final_max_val_over_KD']:.2f}" if isinstance(fitted_params[col]['final_max_val_over_KD'], (int, float)) else fitted_params[col]['final_max_val_over_KD']
            plt.title(f'{col}: Observed (Black) from 1-n vs Initial (Red) and Final (Blue) Model from 2-n\n'
                      f'(Initial: max_val={fitted_params[col]["max_val"]:.2f}, KD={fitted_params[col]["KD"]:.2f})\n'
                      f'(Final: max_val={fitted_params[col]["final_max_val"]:.2f}, KD={fitted_params[col]["final_KD"]:.2f}, max_val/KD={max_val_over_KD_str})')
            plt.legend()
            plt.grid(True)

            if debug_flag:
                timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
                input_base_name = os.path.basename(file_path).replace('.csv', '')
                plot_filename = f"pcr_model_fit--{col}--{input_base_name}--{timestamp}.png"
                plot_path = os.path.join(output_dir, plot_filename)
                plt.savefig(plot_path, bbox_inches='tight')
                print(f"Saved (Debug) plot to {plot_path}")

            plt.show()
            print("Plot displayed successfully.")
else:
    print("Debug display disabled. Skipping plot generation.")

# Outputs: df_fine_tuned, fitted_params, model_predictions

#==================  Cq and GF Combined - Cell-10  ======================
# Cell-10: Optimize seed values using SOS minimization with PCR model propagation
"""
Purpose:
- Initialize seed_optimize values from seed_guess in Cell-9 for each sample.
- Compute modeled value for cycle 1 using seed and propagate values recursively up to the length of df_fine_tuned using the PCR model with refined max_val and KD from Cell-5.
- Optimize seed_optimized to minimize SOS between modeled and observed data.
- Plot observed data (black dots), initial model (blue dashed), and optimized model (green dashed) for visual comparison.
- Report max_val, KD, seed_guess, seed_optimized, initial SOS, and final SOS.

Inputs:
- fitted_params: Dictionary with final fitted parameters (from Cell-5).
- seed_guesses: Dictionary with seed guesses from Cell-9.
- df_fine_tuned: DataFrame with refined fluorescence data (from Cell-5).
- columns_to_fit: List of column names to process (from Cell-2).
- file_path: String, path to the qPCR data file (from Cell-1).
- debug_flag: True or False, indicating Debug mode (from Cell-1).
- debug_display_flag: True or False, indicating Debug plot display (from Cell-1).
- eval_flag: True or False, indicating Evaluation mode (from Cell-1).
- output_dir: String, path to output directory (from Cell-1).

Outputs:
- Printed report of max_val, KD, seed_guess, seed_optimized, initial SOS, and final SOS.
- Plots of observed data (black dots), initial model (blue dashed), and optimized model (green dashed) for each sample (Debug display only).
- (Optional, Debug): Saved report (CSV) and saved plots.
- (Optional, Debug/Evaluation): Global fitting plot (optimized models with observed data dots, single axes).
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import datetime
from scipy.optimize import minimize_scalar

# Verify inputs from Cell-5
try:
    fitted_params, df_fine_tuned
except NameError:
    raise NameError("Required variables (fitted_params, df_fine_tuned) not defined. Please run Cell-5 first.")

# Verify inputs from Cell-9
try:
    seed_guesses
except NameError:
    raise NameError("Required variable (seed_guesses) not defined. Please run Cell-9 first.")

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

print("Optimizing seed values using SOS minimization...")

# PCR model function
def pcr_model(prev, max_val, KD):
    """Single-step PCR model with non-negative constraint"""
    return np.maximum(0, prev * (1 + ((max_val - prev) / max_val) - (prev / (KD + prev))))

# Objective function for SOS minimization
def sos_objective(seed, max_val, KD, observed, n_cycles):
    """Compute SOS between modeled and observed data for a given seed."""
    # Generate full sequence starting from seed
    sequence = []
    prev = seed
    for i in range(n_cycles):
        if i == 0:
            # Cycle 1: use seed as previous value
            current = pcr_model(seed, max_val, KD)
        else:
            # Cycles 2-n: use previous cycle's result
            current = pcr_model(prev, max_val, KD)
        sequence.append(current)
        prev = current
    modeled = np.array(sequence)
    return np.sum((observed - modeled) ** 2)

# Logarithmic objective function for better scaling
def log_sos_objective(log_seed, max_val, KD, observed, n_cycles):
    """Compute SOS using logarithmic seed scaling for better optimization."""
    seed = 10 ** log_seed  # Convert from log space
    return sos_objective(seed, max_val, KD, observed, n_cycles)

# Step 1: Initialize and optimize seed values
print("\n=== Seed Optimization Results ===")
report_data = {
    'Sample': [],
    'max_val': [],
    'KD': [],
    'seed_guess': [],
    'initial_SOS': [],
    'seed_optimized': [],
    'final_SOS': [],
    'SOS_improvement': []
}

valid_columns = [col for col in columns_to_fit if amplification_flags.get(col, True)]
flagged_columns = [col for col in columns_to_fit if not amplification_flags.get(col, True)]

if flagged_columns:
    print(f"Skipped (no substantial amplification): {', '.join(flagged_columns)}")

# Include flagged samples in the report (as N/A), but do not optimize them
for col in columns_to_fit:
    if not amplification_flags.get(col, True):
        report_data['Sample'].append(col)
        report_data['max_val'].append('N/A')
        report_data['KD'].append('N/A')
        report_data['seed_guess'].append('N/A')
        report_data['initial_SOS'].append('N/A')
        report_data['seed_optimized'].append('N/A (no substantial amplification)')
        report_data['final_SOS'].append('N/A')
        report_data['SOS_improvement'].append('N/A')

for col in valid_columns:
    print(f"\nProcessing {col}...")
    
    # Get parameters
    max_val = fitted_params.get(col, {}).get('final_max_val', None)
    KD = fitted_params.get(col, {}).get('final_KD', None)
    seed_guess = seed_guesses.get(col, None)
    
    if max_val is None or KD is None or seed_guess is None:
        print(f"Warning: Missing parameters for {col}. Skipping.")
        report_data['Sample'].append(col)
        report_data['max_val'].append('N/A')
        report_data['KD'].append('N/A')
        report_data['seed_guess'].append('N/A')
        report_data['initial_SOS'].append('N/A')
        report_data['seed_optimized'].append('N/A')
        report_data['final_SOS'].append('N/A')
        report_data['SOS_improvement'].append('N/A')
        continue
    
    # Get observed data
    observed = df_fine_tuned[col].values.copy()
    n_cycles = len(observed)
    
    # Calculate initial SOS
    initial_sos = sos_objective(seed_guess, max_val, KD, observed, n_cycles)
    print(f"  Initial SOS: {initial_sos:.4e}")
    
    # Optimize seed using logarithmic scaling for better handling of large dynamic ranges
    try:
        # Convert seed guess to log space for optimization
        log_seed_guess = np.log10(seed_guess)
        
        # Set bounds in log space (factor of 3 in each direction = 1000x range)
        log_bounds_factor = 3.0
        log_lower_bound = log_seed_guess - log_bounds_factor
        log_upper_bound = log_seed_guess + log_bounds_factor
        
        # Try logarithmic optimization
        result = minimize_scalar(
            log_sos_objective,
            args=(max_val, KD, observed, n_cycles),
            method='bounded',
            bounds=(log_lower_bound, log_upper_bound),
            options={'xatol': 1e-4, 'maxiter': 100}
        )
        
        if result.success:
            log_seed_optimized = result.x
            seed_optimized = 10 ** log_seed_optimized  # Convert back to linear space
            final_sos = result.fun
            sos_improvement = ((initial_sos - final_sos) / initial_sos) * 100
            
            # Check if optimization actually improved the fit
            if final_sos < initial_sos:
                print(f"  Optimized seed: {seed_optimized:.4e} (log: {log_seed_optimized:.4f})")
                print(f"  Final SOS: {final_sos:.4e}")
                print(f"  SOS improvement: {sos_improvement:.2f}%")
            else:
                print(f"  Optimization found worse fit, using seed guess")
                seed_optimized = seed_guess
                final_sos = initial_sos
                sos_improvement = 0.0
        else:
            print(f"  Logarithmic optimization failed: {result.message}")
            # Fall back to seed guess
            seed_optimized = seed_guess
            final_sos = initial_sos
            sos_improvement = 0.0
            
    except Exception as e:
        print(f"  Optimization error: {e}")
        seed_optimized = seed_guess
        final_sos = initial_sos
        sos_improvement = 0.0
    
    # Store results
    report_data['Sample'].append(col)
    report_data['max_val'].append(f"{max_val:.2f}")
    report_data['KD'].append(f"{KD:.2f}")
    report_data['seed_guess'].append(f"{seed_guess:.4e}")
    report_data['initial_SOS'].append(f"{initial_sos:.4e}")
    report_data['seed_optimized'].append(f"{seed_optimized:.4e}")
    report_data['final_SOS'].append(f"{final_sos:.4e}")
    report_data['SOS_improvement'].append(f"{sos_improvement:.2f}%")

# Create report DataFrame
report_df = pd.DataFrame(report_data)
print("\n=== Optimization Summary ===")
print(report_df.to_string(index=False))

# Step 2: Generate plots for each sample
def plot_seed_optimization(col, max_val, KD, seed_guess, seed_optimized, observed, cycles, debug_flag=False, debug_display_flag=False, output_dir="outputs", file_path=""):
    """Plot observed vs modeled data for seed optimization.
    
    Model line starts at cycle 0 (seed position) through max observed cycle.
    Observed data is plotted at actual cycles (e.g., 3-40).
    """
    if not debug_display_flag and not debug_flag:
        print(f"Debug output disabled. Skipping plot for {col}.")
        return
    
    # Model starts at cycle 0, extends to max observed cycle
    max_cycle = int(max(cycles))
    model_cycles = np.arange(0, max_cycle + 1)  # 0, 1, 2, ..., max_cycle
    
    # Generate initial model sequence from cycle 0
    initial_sequence = [seed_guess]  # Cycle 0 = seed
    prev = seed_guess
    for i in range(max_cycle):  # Generate cycles 1 through max_cycle
        current = pcr_model(prev, max_val, KD)
        initial_sequence.append(current)
        prev = current
    
    # Generate optimized model sequence from cycle 0
    optimized_sequence = [seed_optimized]  # Cycle 0 = seed
    prev = seed_optimized
    for i in range(max_cycle):  # Generate cycles 1 through max_cycle
        current = pcr_model(prev, max_val, KD)
        optimized_sequence.append(current)
        prev = current
    
    # Create plot
    fig, ax = plt.subplots(figsize=(12, 8))
    plt.plot(cycles, observed, 'ko', label='Observed Data', markersize=4)
    plt.plot(model_cycles, initial_sequence, 'b--', label=f'Initial Model (seed={seed_guess:.2e})', linewidth=2)
    plt.plot(model_cycles, optimized_sequence, 'g--', label=f'Optimized Model (seed={seed_optimized:.2e})', linewidth=2)
    plt.axhline(y=0, color='gray', linestyle='--', linewidth=0.5)
    plt.xlabel('Cycle')
    plt.ylabel('Fluorescence Value')
    plt.title(f'{col}: Seed Optimization Results\nmax_val={max_val:.2f}, KD={KD:.2f}')
    plt.legend()
    plt.grid(True)
    
    # Save plot in Debug mode
    if debug_flag:
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        input_base_name = os.path.basename(file_path).replace('.csv', '')
        plot_filename = f"seed_optimization--{col}--{input_base_name}--{timestamp}.png"
        plot_path = os.path.join(output_dir, plot_filename)
        plt.savefig(plot_path, bbox_inches='tight')
        print(f"Saved (Debug) plot to {plot_path}")
    
    if debug_display_flag:
        plt.show()
        print(f"Plot displayed for {col}.")
    else:
        plt.close()

# Generate plots for each sample (only those that passed Cell-3)
for col in valid_columns:
    
    max_val = fitted_params.get(col, {}).get('final_max_val', None)
    KD = fitted_params.get(col, {}).get('final_KD', None)
    seed_guess = seed_guesses.get(col, None)
    seed_optimized = report_df.loc[report_df['Sample'] == col, 'seed_optimized'].iloc[0]
    
    if max_val is not None and KD is not None and seed_guess is not None:
        # Convert seed_optimized back to float
        try:
            seed_optimized = float(seed_optimized)
        except:
            seed_optimized = seed_guess
            
        observed = df_fine_tuned[col].values
        cycles = df_fine_tuned.index.values  # Use actual cycle numbers
        plot_seed_optimization(col, max_val, KD, seed_guess, seed_optimized, observed, cycles, debug_flag, debug_display_flag, output_dir, file_path)

# Step 3: Debug output: Save report to CSV (only valid samples)
if debug_flag and report_data['Sample']:
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    input_base_name = os.path.basename(file_path).replace('.csv', '')
    csv_filename = f"seed_optimization_results--{input_base_name}--{timestamp}.csv"
    csv_path = os.path.join(output_dir, csv_filename)
    report_df.to_csv(csv_path, index=False)
    print(f"\nSaved (Debug) optimization results to {csv_path}")

# Step 4: Create final seed_optimized dictionary for use in subsequent cells
seed_optimized_dict = {}
for col in valid_columns:
    if col in report_df['Sample'].values:
        seed_opt = report_df.loc[report_df['Sample'] == col, 'seed_optimized'].iloc[0]
        try:
            seed_optimized_dict[col] = float(seed_opt)
        except:
            seed_optimized_dict[col] = seed_guesses.get(col, None)

print(f"\nSeed optimization completed for {len([s for s in seed_optimized_dict.values() if s is not None])} samples.")

# Step 5 (Evaluation/Debug): Plot global fitting models (optimized only)
def plot_global_fitting_models(seed_optimized_dict, fitted_params, df_fine_tuned, output_dir="outputs", file_path="", save_plot=False, show_plot=False):
    """Plot all optimized model curves on a single axes with observed data dots."""
    if not seed_optimized_dict:
        print("No optimized seeds available for global fitting plot.")
        return

    fig, ax = plt.subplots(figsize=(12, 8))
    for col, seed_opt in seed_optimized_dict.items():
        if seed_opt is None:
            continue
        max_val = fitted_params.get(col, {}).get('final_max_val', None)
        KD = fitted_params.get(col, {}).get('final_KD', None)
        if max_val is None or KD is None:
            continue

        observed_cycles = df_fine_tuned.index.values  # Actual observed cycle numbers (e.g., 3-40)
        observed = df_fine_tuned[col].values
        
        # Model starts at cycle 0, extends to max observed cycle
        max_cycle = int(max(observed_cycles))
        model_cycles = np.arange(0, max_cycle + 1)  # 0, 1, 2, ..., max_cycle
        
        # Generate model from cycle 0
        modeled = [seed_opt]  # Cycle 0 = seed
        prev = seed_opt
        for i in range(max_cycle):  # Generate cycles 1 through max_cycle
            current = pcr_model(prev, max_val, KD)
            modeled.append(current)
            prev = current
        
        ax.plot(observed_cycles, observed, 'ko', markersize=3, alpha=0.6, label=f"{col} Observed")
        ax.plot(model_cycles, modeled, linewidth=2, alpha=0.9, label=f"{col} Model")

    ax.set_xlabel('Cycle')
    ax.set_ylabel('Fluorescence Value')
    ax.set_title('Global Fitting Plot (Optimized Models + Observed)')
    ax.grid(True)
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
    fig.subplots_adjust(right=0.75)

    if save_plot:
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        input_base_name = os.path.basename(file_path).replace('.csv', '')
        plot_filename = f"Global_fitting_plot--{input_base_name}--{timestamp}.png"
        plot_path = os.path.join(output_dir, plot_filename)
        fig.savefig(plot_path, bbox_inches='tight')
        print(f"Saved (Debug/Evaluation) global fitting plot to {plot_path}")

    if show_plot:
        plt.show()
        print("Global fitting plot displayed.")
    else:
        plt.close(fig)

if (eval_flag or debug_flag):
    plot_global_fitting_models(
        seed_optimized_dict,
        fitted_params,
        df_fine_tuned,
        output_dir=output_dir,
        file_path=file_path,
        save_plot=True,
        show_plot=debug_display_flag
    )

# Step 6 (Debug): Export optimized seed-modeled data
if debug_flag and seed_optimized_dict:
    # Model starts at cycle 0, extends to max observed cycle
    max_observed_cycle = int(max(df_fine_tuned.index.values))
    model_cycles = np.arange(0, max_observed_cycle + 1)  # 0, 1, 2, ..., max_cycle
    modeled_df = pd.DataFrame({'Cycle': model_cycles})
    
    for col, seed_opt in seed_optimized_dict.items():
        if seed_opt is None:
            continue
        max_val = fitted_params.get(col, {}).get('final_max_val', None)
        KD = fitted_params.get(col, {}).get('final_KD', None)
        if max_val is None or KD is None:
            continue
        
        # Generate model from cycle 0
        modeled = [seed_opt]  # Cycle 0 = seed
        prev = seed_opt
        for i in range(max_observed_cycle):  # Generate cycles 1 through max_cycle
            current = pcr_model(prev, max_val, KD)
            modeled.append(current)
            prev = current
        modeled_df[col] = modeled

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    input_base_name = os.path.basename(file_path).replace('.csv', '')
    modeled_filename = f"seed_modeled_data--{input_base_name}--{timestamp}.csv"
    modeled_path = os.path.join(output_dir, modeled_filename)
    modeled_df.to_csv(modeled_path, index=False)
    print(f"Saved (Debug) seed-modeled data to {modeled_path}")

# Outputs: report_df, seed_optimized_dict

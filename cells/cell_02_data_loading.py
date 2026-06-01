#==================  Cq and GF Combined - Cell-2   ======================
# Cell-2: Load and prepare qPCR data for qPCR analysis
"""
Purpose:
- Load raw qPCR data from the specified file path (raw data preferred, not machine-corrected input).
- Validate the data format, headers, and data entries.
- Identify columns containing fluorescence data, excluding metadata columns.
- Set the DataFrame index to the 'Cycle' column (if present) or adjust to start at 1.
- Provide a preview of the first 5 data points with the index labeled as 'Cycle' and a summary of maximum fluorescence.
- Plot the raw fluorescence data for all specified columns to visualize the dataset.
- Add detailed debugging output to inspect the loaded data and validation results (Debug mode).

Inputs:
- file_path: String, path to the qPCR data file (from Cell-1).
- eval_flag: True or False, indicating Evaluation mode (from Cell-1).
- debug_flag: True or False, indicating Debug mode (from Cell-1).
- debug_display_flag: True or False, indicating Debug plot display (from Cell-1).
- output_dir: String, path to the outputs directory (from Cell-1).

Outputs:
- df: DataFrame containing the qPCR fluorescence data, with index set to cycles starting at 1.
- columns_to_fit: List of column names with fluorescence data to process in subsequent cells.
- (Optional, Debug): Save the raw data plot as a PNG file.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import datetime

# Verify inputs from Cell-1
try:
    file_path, debug_flag, debug_display_flag, eval_flag, output_dir
except NameError:
    raise NameError("Required variables (file_path, debug_flag, debug_display_flag, eval_flag, output_dir) not defined. Please run Cell-1 first.")

# Load the data
try:
    df = pd.read_csv(file_path)
except FileNotFoundError:
    raise FileNotFoundError(f"Data file not found at {file_path}. Please check the file path and ensure it exists.")
except Exception as e:
    raise Exception(f"Error loading data file: {e}")

# Validate the DataFrame
if df.empty:
    raise ValueError("Loaded DataFrame is empty. Please check the data file and ensure it contains data.")

# Set the index to the 'Cycle' column if present, otherwise adjust the default index
cycle_column = None
for col in df.columns:
    if col.lower() == 'cycle':
        cycle_column = col
        break

if cycle_column is not None:
    try:
        df[cycle_column] = pd.to_numeric(df[cycle_column], errors='coerce')
        if df[cycle_column].isna().any():
            raise ValueError(f"Cycle column '{cycle_column}' contains non-numeric values. Please clean the data.")
        if df[cycle_column].min() != 1:
            print(f"Warning: Cycle column '{cycle_column}' does not start at 1 (min={df[cycle_column].min()}). Adjusting cycles to start at 1.")
            df[cycle_column] = df[cycle_column] - df[cycle_column].min() + 1
        df.set_index(cycle_column, inplace=True)
    except Exception as e:
        raise ValueError(f"Error setting index to Cycle column: {e}")
else:
    df.index = range(1, len(df) + 1)
    print("Warning: No 'Cycle' column found. Using default index starting at 1.")
df.index.name = "Cycle"

# Validate column headers
unnamed_columns = [col for col in df.columns if col.startswith('Unnamed:')]
if unnamed_columns:
    raise ValueError(f"Found columns without headers: {unnamed_columns}. Please provide proper column names.")

# Identify potential metadata columns
metadata_columns = ['cycle', 'index', 'time', 'well', 'sample']
potential_metadata = [col for col in df.columns if col.lower() in metadata_columns]

# Identify fluorescence data columns
columns_to_fit = []
invalid_columns = []
for col in df.columns:
    if col.lower() in metadata_columns:
        continue
    try:
        numeric_series = pd.to_numeric(df[col], errors='coerce')
        if numeric_series.isna().any():
            invalid_rows = df.index[numeric_series.isna()].tolist()
            invalid_columns.append((col, invalid_rows))
        else:
            columns_to_fit.append(col)
    except Exception as e:
        invalid_columns.append((col, f"Error: {e}"))

if invalid_columns:
    error_msg = "Found non-numeric values or errors in the following columns:\n"
    for col, details in invalid_columns:
        error_msg += f" - {col}: {details}\n"
    raise ValueError(error_msg)

if not columns_to_fit:
    raise ValueError("No valid numeric columns found for fluorescence data. Please check the data file for numeric fluorescence values.")

# Check for missing data
if df[columns_to_fit].isna().any().any():
    raise ValueError("Data contains missing values (NaN). Please clean the dataset or replace missing values.")

# Preview: Print the first 5 rows and maximum fluorescence
print("\nFirst 5 rows of loaded data:")
print(df[columns_to_fit].head().to_string(index=True))
print("\nMaximum Fluorescence per Sample (Across All Rows):")
print(df[columns_to_fit].max().round(2).to_string())

# Plot the raw fluorescence data
def plot_raw_data(df, columns, debug_flag=False, debug_display_flag=False, output_dir="outputs", file_path=""):
    """
    Plot the raw fluorescence data for all specified columns.

    Args:
        df (pd.DataFrame): DataFrame with fluorescence data.
        columns (list): List of column names to plot.
        debug_flag (bool): Save debug outputs when True.
        debug_display_flag (bool): Display debug plots when True.
        output_dir (str): Directory to save plots.
        file_path (str): Original file path for naming saved plots.
    """
    if not columns:
        print("No columns to plot.")
        return

    if not debug_display_flag and not debug_flag:
        print("Debug output disabled. Skipping plot generation.")
        return

    fig_size = (12, 8) if len(columns) <= 10 else (15, 10)
    plt.figure(figsize=fig_size, constrained_layout=True)
    cycles = df.index.values  # Use the DataFrame index as cycles
    colors = plt.cm.tab10(np.linspace(0, 1, len(columns)))

    for idx, column in enumerate(columns):
        plt.plot(cycles, df[column], label=column, marker='o', linestyle='-', color=colors[idx])

    plt.axhline(y=0, color='gray', linestyle='--', linewidth=0.5)
    plt.xlabel('Cycle')
    plt.ylabel('Fluorescence Value')
    plt.title('Raw Fluorescence Data: All Cycles')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', ncol=2 if len(columns) > 10 else 1)
    plt.grid(True)

    if debug_flag:
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        input_base_name = os.path.basename(file_path).replace('.csv', '')
        plot_path = os.path.join(output_dir, f'raw_data_plot--{input_base_name}--{timestamp}.png')
        plt.savefig(plot_path, bbox_inches='tight')
        print(f"Saved (Debug) raw data plot to {plot_path}")

    if debug_display_flag:
        plt.show()
        print("Plot displayed successfully.")
    else:
        plt.close()

# Plot the raw data
if columns_to_fit:
    plot_raw_data(df, columns_to_fit, debug_flag, debug_display_flag, output_dir, file_path)
else:
    print("No columns to plot.")

# Outputs: df, columns_to_fit

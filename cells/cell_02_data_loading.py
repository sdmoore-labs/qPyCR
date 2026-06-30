#==================  Cq and GF Combined - Cell-2   ======================
# Cell-2: Load and prepare qPCR data for qPCR analysis
"""
Purpose:
- Load raw qPCR data from the specified file path (raw data preferred, not machine-corrected input).
- Support multiple input formats: CSV (default), RDES (.tsv), and RDML (.rdml).
- Validate the data format, headers, and data entries.
- Identify columns containing fluorescence data, excluding metadata columns.
- Set the DataFrame index to the 'Cycle' column (if present) or adjust to start at 1.
- Provide a preview of the first 5 data points with the index labeled as 'Cycle' and a summary of maximum fluorescence.
- Plot the raw fluorescence data for all specified columns to visualize the dataset.
- Add detailed debugging output to inspect the loaded data and validation results (Debug mode).

Inputs:
- file_path: String, path to the qPCR data file (from Cell-1).
- file_format: String, detected format - 'csv', 'rdes', or 'rdml' (from Cell-1).
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
    file_path, file_format, debug_flag, debug_display_flag, eval_flag, output_dir
except NameError:
    raise NameError("Required variables (file_path, file_format, debug_flag, debug_display_flag, eval_flag, output_dir) not defined. Please run Cell-1 first.")

# ============ Format-specific loading functions ============

def load_csv_data(file_path):
    """Load standard CSV format (Cycle column + sample columns)."""
    return pd.read_csv(file_path)

def load_rdes_data(file_path):
    """
    Load RDES (Real-time PCR Data Essential Spreadsheet) format.
    RDES has wells/samples as rows and cycles as columns.
    This function transposes to qPyCR format (cycles as rows, samples as columns).
    """
    print("Loading RDES format...")
    rdes_df = pd.read_csv(file_path, sep='\t')
    
    # Identify metadata columns vs cycle columns
    metadata_cols = []
    cycle_cols = []
    
    for col in rdes_df.columns:
        try:
            cycle_num = int(col)
            cycle_cols.append((col, cycle_num))
        except (ValueError, TypeError):
            metadata_cols.append(col)
    
    if not cycle_cols:
        raise ValueError("No cycle columns found in RDES file. Expected numeric column headers for cycles.")
    
    # Sort cycle columns by cycle number
    cycle_cols.sort(key=lambda x: x[1])
    cycle_col_names = [c[0] for c in cycle_cols]
    cycle_numbers = [c[1] for c in cycle_cols]
    
    # Determine sample names (use Well column if Sample names are duplicated)
    if 'Sample' in rdes_df.columns:
        if rdes_df['Sample'].duplicated().any() and 'Well' in rdes_df.columns:
            sample_names = [f"{s}_{w}" for s, w in zip(rdes_df['Sample'], rdes_df['Well'])]
        else:
            sample_names = rdes_df['Sample'].values
    elif 'Well' in rdes_df.columns:
        sample_names = rdes_df['Well'].values
    else:
        sample_names = [f"Sample_{i}" for i in range(len(rdes_df))]
    
    # Extract fluorescence data and transpose (samples as rows → samples as columns)
    fluor_data = rdes_df[cycle_col_names].values.T
    
    # Create output DataFrame with cycles as rows
    output_df = pd.DataFrame(fluor_data, index=cycle_numbers, columns=sample_names)
    output_df.index.name = "Cycle"
    
    print(f"RDES import: {len(output_df)} cycles, {len(output_df.columns)} samples")
    print(f"Cycle range: {min(cycle_numbers)} to {max(cycle_numbers)}")
    
    # Reset index to make Cycle a column (for consistency with CSV loading)
    output_df = output_df.reset_index()
    
    return output_df

def load_rdml_data(file_path):
    """
    Load RDML (Real-time PCR Data Markup Language) format.
    Parses the RDML ZIP/XML structure directly to extract raw fluorescence data.
    
    Column names use the well/react ID (unique per reaction) rather than 
    sample names (which may have replicates with the same name).
    """
    print("Loading RDML format...")
    
    import zipfile
    import xml.etree.ElementTree as ET
    
    # RDML namespace
    ns = {'rdml': 'http://www.rdml.org'}
    
    wells_data = {}
    max_cycles = 0
    min_cycles = float('inf')
    
    try:
        with zipfile.ZipFile(file_path, 'r') as zf:
            xml_content = zf.read('rdml_data.xml')
            root = ET.fromstring(xml_content)
            
            # Find experiments
            experiments = root.findall('rdml:experiment', ns)
            if not experiments:
                raise ValueError("No experiments found in RDML file.")
            
            print(f"Found {len(experiments)} experiment(s)")
            
            for exp in experiments:
                exp_id = exp.get('id', 'unknown')
                
                # Find runs in experiment
                runs = exp.findall('rdml:run', ns)
                print(f"Found {len(runs)} run(s) in experiment '{exp_id}'")
                
                for run in runs:
                    # Find reacts (reactions/wells)
                    reacts = run.findall('rdml:react', ns)
                    print(f"Found {len(reacts)} reactions in run")
                    
                    for react in reacts:
                        react_id = react.get('id', 'unknown')
                        
                        # Find data elements (may have multiple if multiplexed)
                        datas = react.findall('rdml:data', ns)
                        
                        for data in datas:
                            # Get target reference (for multiplexed wells)
                            tar_elem = data.find('rdml:tar', ns)
                            target = tar_elem.get('id') if tar_elem is not None else ''
                            
                            # Use react_id as column name (unique well identifier)
                            # Append target only if multiplexed (multiple targets per well)
                            if len(datas) > 1 and target:
                                well_name = f"{react_id}_{target}"
                            else:
                                well_name = react_id
                            
                            # Get amplification data points
                            adps = data.findall('rdml:adp', ns)
                            
                            cycles = []
                            fluor = []
                            for adp in adps:
                                cyc_elem = adp.find('rdml:cyc', ns)
                                fluor_elem = adp.find('rdml:fluor', ns)
                                
                                if cyc_elem is not None and fluor_elem is not None:
                                    try:
                                        cyc_val = int(float(cyc_elem.text))
                                        fluor_val = float(fluor_elem.text)
                                        cycles.append(cyc_val)
                                        fluor.append(fluor_val)
                                    except (ValueError, TypeError):
                                        continue
                            
                            if cycles:
                                wells_data[well_name] = {'cycles': cycles, 'fluor': fluor}
                                max_cycles = max(max_cycles, max(cycles))
                                min_cycles = min(min_cycles, min(cycles))
        
    except zipfile.BadZipFile:
        raise ValueError("Invalid RDML file. Expected ZIP format containing rdml_data.xml")
    
    if not wells_data:
        raise ValueError("No fluorescence data found in RDML file.")
    
    # Build DataFrame - use actual cycle range from data
    all_cycles = list(range(int(min_cycles), int(max_cycles) + 1))
    output_df = pd.DataFrame(index=all_cycles)
    output_df.index.name = "Cycle"
    
    for well_name, data in wells_data.items():
        series = pd.Series(index=data['cycles'], data=data['fluor'])
        output_df[well_name] = series.reindex(all_cycles)
    
    print(f"RDML import: {len(output_df)} cycles (range {min_cycles}-{max_cycles}), {len(output_df.columns)} wells")
    print(f"Wells: {', '.join(list(output_df.columns)[:5])}{'...' if len(output_df.columns) > 5 else ''}")
    
    # Reset index to make Cycle a column
    output_df = output_df.reset_index()
    
    return output_df

# ============ Load data based on detected format ============

try:
    if file_format == 'csv':
        df = load_csv_data(file_path)
    elif file_format == 'rdes':
        df = load_rdes_data(file_path)
    elif file_format == 'rdml':
        df = load_rdml_data(file_path)
    else:
        print(f"Unknown format '{file_format}', attempting CSV load...")
        df = load_csv_data(file_path)
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

# Debug: show what columns we found
print(f"DEBUG: Columns in loaded data: {list(df.columns)[:10]}{'...' if len(df.columns) > 10 else ''}")
print(f"DEBUG: Cycle column found: {cycle_column}")

if cycle_column is not None:
    try:
        print(f"DEBUG: Cycle column values (first 5): {df[cycle_column].head().tolist()}")
        df[cycle_column] = pd.to_numeric(df[cycle_column], errors='coerce')
        if df[cycle_column].isna().any():
            raise ValueError(f"Cycle column '{cycle_column}' contains non-numeric values. Please clean the data.")
        
        # Preserve original cycle numbers from the data file
        start_cycle = int(df[cycle_column].min())
        end_cycle = int(df[cycle_column].max())
        print(f"Note: Data uses cycles {start_cycle} to {end_cycle} (preserving original cycle numbers).")
        
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

# Preview: Print cycle range and first 5 rows
print(f"\n--- Data Summary ---")
print(f"Cycle range: {df.index.min()} to {df.index.max()} ({len(df)} total cycles)")
print(f"Samples: {len(columns_to_fit)}")
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
        input_base_name = os.path.splitext(os.path.basename(file_path))[0]  # Handle any extension
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

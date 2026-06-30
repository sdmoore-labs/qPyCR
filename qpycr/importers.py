"""
qPyCR Data Importers - Support for RDML and RDES formats

This module provides functions to import qPCR data from standard formats
and convert them to the DataFrame format used by qPyCR.

Supported formats:
- RDES (Real-time PCR Data Essential Spreadsheet Format) - TSV files
- RDML (Real-time PCR Data Markup Language) - XML/ZIP files

Usage:
    from qpycr.importers import import_rdes, import_rdml
    
    # Import RDES file
    df = import_rdes("data.tsv")
    
    # Import RDML file (uses built-in XML parsing, no external dependencies)
    df = import_rdml("data.rdml")
"""

import os
import zipfile
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np
from typing import Optional, List, Dict, Union
from pathlib import Path


def import_rdes(file_path: Union[str, Path], 
                sample_column: str = "Sample",
                well_column: str = "Well",
                use_well_as_sample: bool = False,
                target_filter: Optional[str] = None,
                sample_type_filter: Optional[str] = None) -> pd.DataFrame:
    """
    Import qPCR data from an RDES (TSV) file.
    
    RDES format has wells/samples as rows and cycles as columns.
    This function transposes the data to qPyCR format (cycles as rows).
    
    Args:
        file_path: Path to the RDES .tsv file
        sample_column: Column name containing sample identifiers (default: "Sample")
        well_column: Column name containing well positions (default: "Well")
        use_well_as_sample: If True, use well positions instead of sample names
        target_filter: Only include rows with this target (e.g., "Exon 1")
        sample_type_filter: Only include rows with this sample type (e.g., "unkn")
    
    Returns:
        DataFrame with Cycle as index and samples as columns
        
    Example:
        df = import_rdes("RDES_data.tsv")
        # Then use with qPyCR:
        from qpycr import analyze
        # Save to CSV first, or use the analyzer directly
    """
    file_path = Path(file_path)
    
    if not file_path.exists():
        raise FileNotFoundError(f"RDES file not found: {file_path}")
    
    # Read the TSV file
    rdes_df = pd.read_csv(file_path, sep='\t')
    
    # Identify metadata columns (non-numeric, before the cycle data)
    metadata_cols = []
    cycle_cols = []
    
    for col in rdes_df.columns:
        # Try to interpret column as cycle number
        try:
            cycle_num = int(col)
            cycle_cols.append((col, cycle_num))
        except (ValueError, TypeError):
            metadata_cols.append(col)
    
    if not cycle_cols:
        raise ValueError("No cycle columns found in RDES file. "
                        "Expected numeric column headers for cycles.")
    
    # Sort cycle columns by cycle number
    cycle_cols.sort(key=lambda x: x[1])
    cycle_col_names = [c[0] for c in cycle_cols]
    cycle_numbers = [c[1] for c in cycle_cols]
    
    # Apply filters if specified
    filtered_df = rdes_df.copy()
    
    if target_filter and 'Target' in filtered_df.columns:
        filtered_df = filtered_df[filtered_df['Target'] == target_filter]
    
    if sample_type_filter and 'Sample Type' in filtered_df.columns:
        filtered_df = filtered_df[filtered_df['Sample Type'] == sample_type_filter]
    
    if filtered_df.empty:
        raise ValueError("No data remaining after applying filters.")
    
    # Determine sample names
    if use_well_as_sample and well_column in filtered_df.columns:
        sample_names = filtered_df[well_column].values
    elif sample_column in filtered_df.columns:
        # Combine sample and well if samples are not unique
        if filtered_df[sample_column].duplicated().any():
            if well_column in filtered_df.columns:
                sample_names = [f"{s}_{w}" for s, w in 
                               zip(filtered_df[sample_column], filtered_df[well_column])]
            else:
                sample_names = [f"{s}_{i}" for i, s in 
                               enumerate(filtered_df[sample_column])]
        else:
            sample_names = filtered_df[sample_column].values
    else:
        # Fall back to well positions or index
        if well_column in filtered_df.columns:
            sample_names = filtered_df[well_column].values
        else:
            sample_names = [f"Sample_{i}" for i in range(len(filtered_df))]
    
    # Extract fluorescence data and transpose
    fluor_data = filtered_df[cycle_col_names].values.T
    
    # Create output DataFrame with cycles as rows (preserve original cycle numbers)
    output_df = pd.DataFrame(
        fluor_data,
        index=cycle_numbers,
        columns=sample_names
    )
    output_df.index.name = "Cycle"
    
    min_cycle = min(cycle_numbers)
    max_cycle = max(cycle_numbers)
    
    print(f"Imported RDES file: {len(output_df)} cycles (range {min_cycle}-{max_cycle}), {len(output_df.columns)} samples")
    if min_cycle != 1:
        print(f"Note: Data starts at cycle {min_cycle} (preserving original cycle numbers)")
    print(f"Samples: {', '.join(str(s) for s in output_df.columns[:5])}"
          f"{'...' if len(output_df.columns) > 5 else ''}")
    
    return output_df


def import_rdml(file_path: Union[str, Path],
                use_well_id: bool = True) -> pd.DataFrame:
    """
    Import qPCR data from an RDML file.
    
    Uses built-in XML parsing - no external dependencies required.
    RDML files are ZIP archives containing rdml_data.xml.
    
    Args:
        file_path: Path to the RDML file (.rdml)
        use_well_id: If True, use well/react ID as column names (unique per reaction).
                     If False, use sample_target names (may have duplicates for replicates).
    
    Returns:
        DataFrame with Cycle as index and samples/wells as columns
        
    Example:
        df = import_rdml("data.rdml")
    """
    file_path = Path(file_path)
    
    if not file_path.exists():
        raise FileNotFoundError(f"RDML file not found: {file_path}")
    
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
                        
                        # Get sample reference (for non-well-id naming)
                        sample_elem = react.find('rdml:sample', ns)
                        sample_id = sample_elem.get('id') if sample_elem is not None else react_id
                        
                        # Find data elements (may have multiple if multiplexed)
                        datas = react.findall('rdml:data', ns)
                        
                        for data in datas:
                            # Get target reference (for multiplexed wells)
                            tar_elem = data.find('rdml:tar', ns)
                            target = tar_elem.get('id') if tar_elem is not None else ''
                            
                            # Determine column name
                            if use_well_id:
                                # Use react_id (unique well identifier)
                                # Append target only if multiplexed (multiple targets per well)
                                if len(datas) > 1 and target:
                                    well_name = f"{react_id}_{target}"
                                else:
                                    well_name = react_id
                            else:
                                # Use sample_target (may have replicates)
                                well_name = f"{sample_id}_{target}" if target else sample_id
                            
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
    
    # Build DataFrame - use actual cycle range from data (preserve original cycle numbers)
    all_cycles = list(range(int(min_cycles), int(max_cycles) + 1))
    output_df = pd.DataFrame(index=all_cycles)
    output_df.index.name = "Cycle"
    
    for well_name, data in wells_data.items():
        series = pd.Series(index=data['cycles'], data=data['fluor'])
        output_df[well_name] = series.reindex(all_cycles)
    
    print(f"Imported RDML file: {len(output_df)} cycles (range {int(min_cycles)}-{int(max_cycles)}), {len(output_df.columns)} wells")
    if min_cycles != 1:
        print(f"Note: Data starts at cycle {int(min_cycles)} (preserving original cycle numbers)")
    print(f"Wells: {', '.join(str(s) for s in output_df.columns[:5])}"
          f"{'...' if len(output_df.columns) > 5 else ''}")
    
    return output_df


def save_as_qpycr_csv(df: pd.DataFrame, output_path: Union[str, Path]) -> str:
    """
    Save a DataFrame in qPyCR-compatible CSV format.
    
    Args:
        df: DataFrame with Cycle as index and samples as columns
        output_path: Output file path
        
    Returns:
        Path to saved file
    """
    output_path = Path(output_path)
    
    # Reset index to make Cycle a column
    output_df = df.reset_index()
    
    # Ensure Cycle column name
    if output_df.columns[0] != "Cycle":
        output_df = output_df.rename(columns={output_df.columns[0]: "Cycle"})
    
    output_df.to_csv(output_path, index=False)
    print(f"Saved qPyCR-compatible CSV: {output_path}")
    
    return str(output_path)

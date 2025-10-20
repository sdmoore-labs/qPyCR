# Inputs Folder

Place your qPCR data files (CSV format) in this folder for analysis.

## Supported File Formats
- CSV files with Cycle column and sample columns
- Raw fluorescence data (not machine-corrected)
- Standard qPCR machine export formats

## Example File Structure
```
Cycle,A1,A2,A3,B1,B2,B3
1,1000,1200,1100,950,1050,1000
2,1020,1220,1120,970,1070,1020
3,1050,1250,1150,1000,1100,1050
...
```

## Column Naming Guidelines

### ✅ Recommended Names
- **Simple**: `A1`, `B2`, `Sample1`, `Well_1`
- **Descriptive**: `Sample_A1`, `Control_1`, `Treatment_1`
- **With underscores**: `Sample_A1`, `Well_1`, `Control_1`

### ⚠️ Use with Caution
- **Spaces**: `"Sample A"`, `"Well 1"` (works but not recommended)
- **Special characters**: `"Sample-A"`, `"Well.1"` (may cause issues in plots)
- **Very long names**: `"Very_Long_Sample_Name_That_Might_Cause_Display_Issues"`

### ❌ Avoid
- **Reserved words**: `Cycle`, `Index`, `Time`, `Well`, `Sample` (case-insensitive)
- **Empty names**: Unnamed columns will cause errors
- **Non-ASCII characters**: May cause encoding issues

## Notes
- Ensure your CSV has a 'Cycle' column (case-insensitive)
- Sample columns should contain numeric fluorescence values
- Missing values will be flagged during analysis
- Column names are preserved in output files

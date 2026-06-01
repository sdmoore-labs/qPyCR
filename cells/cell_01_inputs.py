#==================  Cq and GF Combined - Cell-1  ======================
# Cell-1: Collect user inputs for qPCR Analysis (Cq and GF Combined)
"""
Purpose:
- Collect user inputs for the qPCR analysis pipeline to produce both Cq and Global-fitting (GF) outputs.
- Gather the path to the qPCR data file (supports root directory filenames in Colab, e.g., 'qpcr_data.csv').
- Assume input data are "raw" (preferred: unadjusted, not platform-corrected).
- Assign optional output flags: Evaluation (-e) and Debug (-d).
- Validate inputs with explanations in prompts and provide retry loops for error recovery (up to 3 attempts per input).
- Create outputs directory for generated files.
- Store validated inputs for use in subsequent cells.

Default Behavior (press Enter for each flag prompt):
- Evaluation (-e) = False: No extra outputs beyond the final report
- Debug (-d) = False: No intermediate files/plots, only final CSV from Cell-11

Optional Flags:
- Enter -e: Evaluation mode (limited set of helpful outputs for inspection)
- Enter -d: Debug mode (extensive set of intermediate data and graphical outputs)

Outputs:
- file_path: String, path to the qPCR data file.
- eval_flag: True or False, indicating whether to generate evaluation outputs (default: False).
- debug_flag: True or False, indicating whether to generate debug outputs (default: False).
- debug_display_flag: True or False, indicating whether to display debug plots (default: False).

- output_dir: String, path to the outputs directory.
"""

import sys
import os

def get_input_with_retry(prompt, validator=None, max_retries=3):
    """
    Helper function to get input with retries on invalid responses.
    - prompt: String to display to user.
    - validator: Optional function that takes input and returns True if valid, False otherwise.
    - max_retries: Max attempts before raising error.
    """
    for attempt in range(max_retries):
        user_input = input(prompt).strip()
        if validator is None or validator(user_input):
            return user_input
        print(f"Invalid input. Please try again. (Attempt {attempt + 1}/{max_retries})")
    raise ValueError(f"Max retries exceeded for input: {prompt}")

# Create outputs directory
output_dir = "outputs"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    print(f"Created outputs directory: {output_dir}")
else:
    print(f"Using existing outputs directory: {output_dir}")

# Collect and validate file path with retry
file_prompt = "Enter the path to your raw qPCR data file (e.g., 'qpcr_data.csv' for root directory in Colab, or full path like '/content/qpcr_data.csv'):\n"
file_path = get_input_with_retry(file_prompt, lambda fp: os.path.isfile(fp))

# Collect evaluation flag (-e)
eval_prompt = "Enter -e for Evaluation output (limited key files + plots), or press Enter for none [default: none]:\n"
eval_input = get_input_with_retry(eval_prompt, lambda ei: ei.lower() in ['-e', ''] or ei == '').lower()
eval_flag = True if eval_input == '-e' else False

# Collect debug flag (-d)
debug_prompt = "Enter -d for Debug output (full intermediate files + plots), or press Enter for minimal output [default: minimal]:\n"
debug_input = get_input_with_retry(debug_prompt, lambda di: di.lower() in ['-d', ''] or di == '').lower()
debug_flag = True if debug_input == '-d' else False

# Debug display flag (intermediate plots shown only in Debug)
debug_display_flag = debug_flag

# Print collected inputs for confirmation
print("\nCollected Inputs:")
print(f"File path: {file_path}")
print(f"Evaluation flag: {eval_flag}")
print(f"Debug flag: {debug_flag}")
print(f"Output directory: {output_dir}")
print("Note: Input data are assumed to be raw (unadjusted, not machine-corrected).")

# Outputs: file_path, debug_flag, debug_display_flag, eval_flag, output_dir

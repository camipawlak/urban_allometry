#!/bin/bash
# run_script.sh - Easily run Python scripts in the allometry conda environment
# Usage: ./run_script.sh <script_name> [script_args...]
# Example: ./run_script.sh 0_clean_names.py
# Example: ./run_script.sh 0_clean_names.py --data-dir data/_raw

# Get the directory of this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Check if user provided a script name
if [ -z "$1" ]; then
    echo "Usage: ./run_script.sh <script_name> [script_args...]"
    echo "Example: ./run_script.sh 0_clean_names.py"
    echo ""
    echo "Available scripts in scripts/py/:"
    ls -1 "$SCRIPT_DIR/scripts/py/"*.py 2>/dev/null || echo "  (no scripts found)"
    exit 1
fi

# Activate the allometry conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate allometry 2>/dev/null

if [ $? -ne 0 ]; then
    echo ""
    echo "ERROR: Could not activate 'allometry' conda environment."
    echo "Please create it first by running:"
    echo "  conda env create -f environment.yml"
    exit 1
fi

# Run the script from the repo root (shift so the script name is not passed as an argument)
SCRIPT="$1"
shift
cd "$SCRIPT_DIR"
python "scripts/py/$SCRIPT" "$@"

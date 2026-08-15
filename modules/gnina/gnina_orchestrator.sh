#!/bin/bash
# Enable strict error handling for the orchestrator setup
set -e

# =====================================================
# 1. Directory Anchoring
# =====================================================
# Forces the script to run from the directory it physically lives in, 
# preventing pathing errors if launched from ~ or /mnt/ceph/
cd "$(dirname "$0")"

# =====================================================
# 2. Input Validation
# =====================================================
if [ "$#" -ne 1 ]; then
    echo "=========================================================="
    echo " ERROR: Missing run list."
    echo " USAGE: ./gnina_orchestrator.sh <run_list.txt>"
    echo " EXAMPLE: ./gnina_orchestrator.sh yef3_grid_benchmark.txt"
    echo "=========================================================="
    exit 1
fi

RUN_LIST="$1"

# Verify the text file actually exists in the folder
if [ ! -f "$RUN_LIST" ]; then
    echo "ERROR: The file '$RUN_LIST' does not exist in $(pwd)."
    exit 1
fi

# Automatically ensure the worker script has execution privileges
chmod +x GNINA_worker.sh

echo "=========================================================="
echo " Starting orchestrator with run list: $RUN_LIST"
echo "=========================================================="

# =====================================================
# 3. The Execution Loop
# =====================================================
line_num=0

# Read the file line by line
while IFS= read -r line || [[ -n "$line" ]]; do
    

    # Skip completely empty lines (allows for visual spacing in the text file)
    if [[ -z "$line" ]]; then
        continue
    fi

    # Skip any line starting with # (ignoring leading whitespace)
    if [[ "$line" =~ ^[[:space:]]*# ]]; then
        continue
    fi
    line_num=$((line_num + 1))
    
    echo " "
    echo ">>> [Job $line_num] ./GNINA_worker.sh $line"
    
    # Temporarily disable strict error handling for the execution step.
    # If a specific GNINA run crashes (e.g., due to a corrupted molecule), 
    # this ensures the orchestrator survives and moves to the next run.
    set +e
    
    # We leave $line unquoted here so Bash automatically splits the text 
    # string into $1, $2, $3, etc., for the worker script.
    ./GNINA_worker.sh $line
    
    EXIT_CODE=$?
    
    # Re-enable strict error handling
    set -e

    if [ $EXIT_CODE -ne 0 ]; then
        echo "WARNING: [Job $line_num] GNINA_worker.sh crashed with exit code $EXIT_CODE."
        echo "Failed parameters: $line"
        sleep 2
    else
        echo ">>> [Job $line_num] Finished successfully."
    fi

done < "$RUN_LIST"

echo "=========================================================="
echo " Orchestrator finished processing jobs in $RUN_LIST."
echo "=========================================================="
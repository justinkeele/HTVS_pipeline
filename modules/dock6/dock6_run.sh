#!/bin/bash

# Define directories
INPUT_DIR="Dock6_ready_molecules"
OUTPUT_DIR="Docking_Results_HEAT"
TEMPLATE="flex_docking.in"

# Create the output folder so your workspace stays clean
mkdir -p "$OUTPUT_DIR"

echo "Starting DOCK6 Flexible Search..."

# Loop through every prepped molecule
for mol2_file in "$INPUT_DIR"/*.mol2; do
    
    # Skip if folder is empty
    [ -e "$mol2_file" ] || continue
    
    # Extract the base name (e.g., "tophit_high_Z1711599189")
    base_name=$(basename "$mol2_file" .mol2)
    
    echo "Docking $base_name"
    
    # Define the unique output prefix for this specific molecule
    unique_output="$OUTPUT_DIR/${base_name}_docked"
    
    # Make a temporary copy of the template, injecting the specific file paths
    # The '|' symbol is used as a delimiter so file paths with '/' don't break the sed command
    sed -e "s|LIGAND_FILE_PATH|$mol2_file|g" \
        -e "s|OUTPUT_PREFIX|$unique_output|g" \
        "$TEMPLATE" > current_run.in
        
    # Run DOCK6 using the dynamically generated parameter file
    dock6 -i current_run.in -o "${unique_output}_terminal_log.out"
    
done

# Clean up the temporary file
rm current_run.in

echo "All validation molecules docked successfully!"
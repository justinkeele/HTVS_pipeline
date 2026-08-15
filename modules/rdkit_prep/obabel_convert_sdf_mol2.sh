#!/bin/bash

# Define your folders
INPUT_DIR="Step_1_rdkit_72pH_molecules"
OUTPUT_DIR="Step_2_72pH_molecules"

# 1. Guarantee the output folder exists before starting
mkdir -p "$OUTPUT_DIR"

# 2. Loop through every .sdf file in the 2D folder
for file in "$INPUT_DIR"/*.sdf; do
    
    # Safety check: If the folder is empty, break the loop
    [ -e "$file" ] || continue

    # 3. Strip the path and extension to get the raw ID (e.g., Z1267337443)
    filename=$(basename "$file")
    id="${filename%.*}"

    echo "Converting $id to mol2"

    # 4. The OpenBabel Execution (Replace flags as needed)
    obabel -isdf "$file" -omol2 -O "$OUTPUT_DIR/${id}.mol2"
done

echo ""
echo "If it worked, 3D molecules are waiting in $OUTPUT_DIR/"

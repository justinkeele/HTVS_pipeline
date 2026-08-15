#!/bin/bash

# Define your folders
INPUT_DIR="2D_molecules"
OUTPUT_DIR="meeko_ready_pipe"

# 1. Guarantee the output folder exists before starting
mkdir -p "$OUTPUT_DIR"

# 2. Loop through every .sdf file in the 2D folder
for file in "$INPUT_DIR"/*.sdf; do
    
    # Safety check: If the folder is empty, break the loop
    [ -e "$file" ] || continue

    # 3. Strip the path and extension to get the raw ID (e.g., Z1267337443)
    filename=$(basename "$file")
    id="${filename%.*}"

    echo "Converting $id to 3D..."

    # 4. The OpenBabel Execution (Replace flags as needed)
    obabel -isdf "$file" -osdf  --p 7.2 | obabel -isdf - -osdf -h | obabel -isdf - -osdf -O "$OUTPUT_DIR/${id}.sdf" --gen3d --best
done

echo ""
echo "If it worked, 3D molecules are waiting in $OUTPUT_DIR/"

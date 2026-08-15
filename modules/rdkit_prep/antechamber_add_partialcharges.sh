#!/bin/bash

#Export amber so it works
export AMBERHOME="/home/justin/dock6/bin/antechamber"
export PATH="$AMBERHOME/bin:$PATH"

# Define your folders
INPUT_DIR="Step_2_72pH_molecules"
OUTPUT_DIR="Dock6_ready_molecules"

# 1. Guarantee the output folder exists before starting
mkdir -p "$OUTPUT_DIR"

# 2. Loop through every .sdf file in the 2D folder
for file in "$INPUT_DIR"/*.mol2; do
    
    # Safety check: If the folder is empty, break the loop
    [ -e "$file" ] || continue

    # 3. Strip the path and extension to get the raw ID (e.g., Z1267337443)
    filename=$(basename "$file")
    id="${filename%.*}"

    echo "Converting $id to mol2"

    # 4. The Antechamber Execution (Replace flags as needed)
    antechamber -i "$file" \
        -fi mol2 \
        -o "$OUTPUT_DIR/$id.mol2" \
        -fo mol2 \
        -c bcc \
        -s 2

done > antechamber.log

echo ""
echo "If it worked, 3D molecules are waiting in $OUTPUT_DIR/"

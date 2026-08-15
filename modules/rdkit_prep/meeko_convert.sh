#!/bin/bash

# Define the pipeline folders
INPUT_DIR="rdkit_prepared_47pH_no_inhibition_molecules"
OUTPUT_DIR="47pH_vina_ready_no_inhibition_molecules"

# 1. Guarantee the output folder exists
mkdir -p "$OUTPUT_DIR"

echo "Converting to .pdbqt..."
echo ""

# 2. Loop through every 3D .sdf file
for file in "$INPUT_DIR"/*.sdf; do
    
    # Safety check: bypass if the folder is empty
    [ -e "$file" ] || continue

    # 3. Extract the clean molecule ID (e.g., tophit_high_Z123)
    filename=$(basename "$file")
    id="${filename%.*}"
    
    echo "Converting $file"
    # 4. The Meeko Execution
    # By default, this assigns Gasteiger charges, AutoDock atom types, and torsion trees.
    mk_prepare_ligand.py -i "$file" -o "$OUTPUT_DIR/${id}.pdbqt"
    
done

echo ""
echo "Maybe it worked!"

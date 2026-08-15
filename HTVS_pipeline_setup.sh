#!/bin/bash


# This block generates the directory architecture

# Use the first command-line argument as the path, OR default to $PWD/HTVS_pipeline
WORKSPACE_DIR="${1:-$PWD/HTVS_pipeline}"
echo "Generating your directories in: $WORKSPACE_DIR"

# Core Executable Script Modules
mkdir -p "${WORKSPACE_DIR}/modules/gnina"
mkdir -p "${WORKSPACE_DIR}/modules/rxdock"
mkdir -p "${WORKSPACE_DIR}/modules/rdkit_prep"
mkdir -p "${WORKSPACE_DIR}/modules/ZINC_download"

# Ligand Libraries
mkdir -p "${WORKSPACE_DIR}/libraries/zinc_subsets"
mkdir -p "${WORKSPACE_DIR}/libraries/enamine_subsets"

# Target Protein Directories
TARGET_DIR="${WORKSPACE_DIR}/targets/yef3"
mkdir -p "${TARGET_DIR}/receptors"
mkdir -p "${TARGET_DIR}/configs/gnina"
mkdir -p "${TARGET_DIR}/ligand_batches/base_test_batch/smiles_pH_H_gen3D_output"

#Results Directory
mkdir -p "${WORKSPACE_DIR}/results/yef3/gnina"

echo "Directory creation is done."
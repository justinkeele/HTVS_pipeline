#!/usr/bin/env python3

import sys
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import RDLogger

INPUT_FILENAME = "all_mols.smi"
OUTPUT_FILENAME = "all_deduplicated_mols.smi"

# Disable RDKit's annoying warning text wall for weird molecules
RDLogger.DisableLog('rdApp.*')

def deduplicate_smiles(input_file, output_file):
    # Initialize the tool that strips away protonation state differences
    uncharger = rdMolStandardize.Uncharger()
    
    seen_molecules = set()
    unique_records = []
    
    print(f"Reading {input_file}...")
    
    with open(input_file, 'r') as f:
        lines = f.readlines()
        
    for line in lines:
        parts = line.strip().split('\t')
        if len(parts) < 2:
            continue
            
        original_smiles = parts[0]
        zinc_id = parts[1]
        
        # 1. Convert text to a chemical graph
        mol = Chem.MolFromSmiles(original_smiles)
        
        # Skip invalid SMILES strings that might have corrupted during download
        if mol is None:
            continue 
            
        # 2. Neutralize the molecule to collapse different protonation states
        try:
            neutral_mol = uncharger.uncharge(mol)
        except Exception:
            # If the uncharger fails on a bizarre metal complex, just use the original
            neutral_mol = mol 
            
        # 3. Generate a strict, canonical SMILES keeping stereochemistry
        # Because we didn't touch tautomers, they remain distinct
        canonical_isomeric_smiles = Chem.MolToSmiles(neutral_mol, isomericSmiles=True)
        
        # 4. Check for duplicates
        if canonical_isomeric_smiles not in seen_molecules:
            seen_molecules.add(canonical_isomeric_smiles)
            # Save the original line so you don't lose the exact ZINC ID
            unique_records.append(f"{original_smiles}\t{zinc_id}\n")

    # Write the cleaned list to a new file
    with open(output_file, 'w') as out_f:
        out_f.writelines(unique_records)
        
    print(f"Started with: {len(lines)} molecules.")
    print(f"Ended with:   {len(unique_records)} unique molecules.")
    print(f"Removed:      {len(lines) - len(unique_records)} duplicates/protonation variants.")
    print(f"Saved to:     {output_file}")
    
deduplicate_smiles(INPUT_FILENAME, OUTPUT_FILENAME)
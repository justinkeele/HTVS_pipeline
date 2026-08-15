#!/usr/bin/env python3

import glob
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from dimorphite_dl import protonate_smiles

# Define your folders
input_dir = "2D_molecules"
output_dir = "Step_1_rdkit_72pH_molecules"

# Make the output folder if it doesn't exist (like 'mkdir -p')
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Tell glob to look specifically inside the input folder
search_path = os.path.join(input_dir, "*.sdf")
# Grab the files
file_list = glob.glob(search_path)


for current_file in file_list:
    # Rename to _3D.sdf
    base_name = os.path.basename(current_file)
    new_name = base_name.replace(".sdf", ".sdf")
    output_path = os.path.join(output_dir, new_name)

    supplier = Chem.SDMolSupplier(current_file)
    original_molecule = supplier[0] 
    
    if original_molecule is not None:
        # Split the molecule into all of its disconnected fragments
        fragments = Chem.GetMolFrags(original_molecule, asMols=True)
        
        # Pick the fragment with the highest atom count
        largest_fragment = max(fragments, key=lambda frag: frag.GetNumAtoms())

        # Feed the isolated drug fragment into the SMILES generator
        smiles_string = Chem.MolToSmiles(largest_fragment)

        # Adjust the protonation state for pH 7.2 FIRST
        # Dimorphite returns a list of possibilities. We grab the most thermodynamicaplly stable one [0].
        protonated_smiles = protonate_smiles(smiles_string, ph_min=7.2, ph_max=7.2, max_variants=1)
        


        if protonated_smiles:
            protonated_molecule = Chem.MolFromSmiles(protonated_smiles[0])

            # Then add Hydrogens
            molecule_with_h = Chem.AddHs(protonated_molecule)
        
            # Generate the 3D coordinates
            AllChem.EmbedMolecule(molecule_with_h, AllChem.ETKDGv3())
        
            # Write to the new file
            writer = Chem.SDWriter(output_path)
            writer.write(molecule_with_h)
            writer.close()

print("3D Batch complete!")

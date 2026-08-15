#!/usr/bin/env python3

import os
from rdkit import Chem
from rdkit.Chem import AllChem
from dimorphite_dl import protonate_smiles
import time
import datetime

# Define your folders
input_file = "/home/justin/Jagdish_lab/YeF3/base_test_batch/new_experimental_2D_smiles.smi"
output_dir = "/home/justin/Jagdish_lab/YeF3/base_test_batch/new_experimental_smiles_pH_H_gen3D_output"
target_ph = 7.2
chunk_size = 10000   
#comment out before HPC run, or set ridiculously high ie: 1,000,000
#max_molecules_to_test = 50

# Make the output folder if it doesn't exist (like 'mkdir -p')
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 3D gen parameters
params = AllChem.ETKDGv3()
params.randomSeed = 42       # Keeps results reproducible (Dr. Patel will like this)



writer = None
batch_index = 1
# Molecule trackers
mol_counter = 0

# Error log helper
def log_error(smiles, zinc_id, error_type, batch_idx):
    batch_folder = os.path.join(output_dir, f"input_batch_{batch_index - 1}")
    os.makedirs(batch_folder, exist_ok=True)
    with open(os.path.join(batch_folder, "failed_mols.txt"), "a") as f:
        f.write(f"{smiles}\t{zinc_id}\t{error_type}\n")

read_line_failure = 0
salt_strip_fail = 0
dimorphite_crash = 0
failed_to_protonate = 0
hydrogen_fail = 0
fail_3D_gen = 0

print(f"Starting 3D gen on {input_file}")
total_start_time = time.perf_counter()
chunk_start_time = time.perf_counter()

# Using 'with open' is RAM-efficient because it only loads one line into memory at a time
with open(input_file, 'r') as f:
    for line in f: # enumerate gives us the line number for debugging
        line = line.strip()
        if not line: 
            read_line_failure += 1
            continue

        # Split the SMILES and the ZINC ID based on the tab character
        parts = line.split()
        if len(parts) != 2: 
            read_line_failure += 1
            continue
        
        smiles_string = parts[0]
        zinc_id = parts[1]

        # Strip Salts
        try:
            # Briefly convert to 2D to measure fragments
            temp_mol = Chem.MolFromSmiles(smiles_string)
            if temp_mol:

                # Split molecule into disconnected fragments (breaks at the '.' in SMILES)
                fragments = Chem.GetMolFrags(temp_mol, asMols=True)
                # Keep only the fragment with the most atoms (the main drug)
                largest_fragment = max(fragments, key=lambda frag: frag.GetNumAtoms())
                # Convert the cleaned drug back into a SMILES string for Dimorphite
                smiles_string = Chem.MolToSmiles(largest_fragment)
            else:
                log_error(smiles_string, zinc_id, "Salt_stripping_failed", batch_index)
                salt_strip_fail += 1
                continue
        except Exception as e:
            log_error(smiles_string, zinc_id, "Salt_stripping_failed", batch_index)
       
       
        #Protonate the de-salted smiles
        try:
            # Dimorphite-DL calculates the pKa at target pH
            protonated_smiles_list = protonate_smiles(smiles_string, ph_min=target_ph, ph_max=target_ph, max_variants=1)
            if not protonated_smiles_list:
                log_error(smiles_string, zinc_id, "Protonation_failed", batch_index)
                dimorphite_crash += 1
                continue
            prot_smiles = protonated_smiles_list[0]
        except Exception as e:
            log_error(smiles_string, zinc_id, "Protonation_failed", batch_index)
            failed_to_protonate += 1
            continue


        #Add hydrogens and make 2D
        try:
            mol_h = Chem.AddHs(Chem.MolFromSmiles(prot_smiles))
        except Exception as e:
            log_error(smiles_string, zinc_id, "Adding_Hydrogens_failed", batch_index)
            hydrogen_fail += 1
            continue

        # Name the molecule by ZINC ID for meeko later
        mol_h.SetProp("_Name", zinc_id)

        # 3D generation
        try:
            if AllChem.EmbedMolecule(mol_h, params) == -1:
                log_error(smiles_string, zinc_id, "3D_generation_failed", batch_index)
                fail_3D_gen += 1
                continue
        except Exception as e:
            log_error(smiles_string, zinc_id, "3D_generation_failed", batch_index)
            fail_3D_gen += 1
            continue

        # Write chunk files, if we hit a multiple of chunk size, ie 10,000
        # then close the current file and open a new one
        if mol_counter % chunk_size == 0:
            if writer:
                writer.close()
            chunk_time = str(datetime.timedelta(seconds=int(time.perf_counter() - chunk_start_time)))
            print(f"Batch took {chunk_time}")
            batch_folder = os.path.join(output_dir, f"input_batch_{batch_index}")
            os.makedirs(batch_folder, exist_ok=True)
            writer = Chem.SDWriter(os.path.join(batch_folder, f"chunk_{batch_index}.sdf"))
            batch_index += 1
            chunk_start_time = time.perf_counter()

        # Write each 3D mol to .sdf
        mol_h.SetProp("_Name", f"{zinc_id}")
        writer.write(mol_h)
        mol_counter += 1            

        #  # If we have reached our max test batches, stop reading the file
        # if max_molecules_to_test is not None and mol_counter >= max_molecules_to_test:
        #     print(f"\n--- Reached test limit of {max_molecules_to_test} molecules. Ending run. ---")
        #     break

#close the last file when loop finishes
if writer: writer.close()
    
total_time = str(datetime.timedelta(seconds=int(time.perf_counter() - total_start_time)))
print(f"===========================================")
print(f"{batch_index - 1} batches created")
print(f"-------------------------------------------")
print(f"{read_line_failure} failed to be read")
print(f"{salt_strip_fail} failed salt removal")
print(f"{failed_to_protonate} failed protonation")
print(f"{dimorphite_crash} crashed dimorphite.dl")
print(f"{hydrogen_fail} failed to add hydrogens")
print(f"{fail_3D_gen} failed 3D generation")
print(f"-------------------------------------------")
print(f"{mol_counter} molecules were converted")
print(f"Total time was {total_time}")
print(f"===========================================")
import pandas as pd
import glob
import os
import csv

print("Compiling detailed scored blocks...")

# 1. Find all TSV files and SORT them alphabetically first!
# This ensures our blocks are grouped together logically before we add blank rows.
tsv_files = sorted(glob.glob("*/all_ligands_master_summary.tsv"))

all_runs_data = []

for file in tsv_files:
    run_name = os.path.dirname(file)
    
    # Read the raw data
    df = pd.read_csv(file, sep='\t', names=['Ligand', 'Vina_Affinity', 'CNN_Pose_Score', 'CNN_Affinity'])
    df = df.dropna()
    
    # Calculate Consensus Ranks
    df['Pose_Rank'] = df['CNN_Pose_Score'].rank(ascending=False, method='min')
    df['Affinity_Rank'] = df['CNN_Affinity'].rank(ascending=False, method='min')
    df['Consensus_Rank'] = (df['Pose_Rank'] + df['Affinity_Rank']) / 2.0
    
    # 2. Sort the data 4 different ways and STRIP THE INDEX
    # reset_index(drop=True) is crucial so they all cleanly align from row 0 to 23
    cons_df = df.sort_values(by=['Consensus_Rank', 'CNN_Pose_Score'], ascending=[True, False]).reset_index(drop=True)
    pose_df = df.sort_values(by='CNN_Pose_Score', ascending=False).reset_index(drop=True)
    aff_df  = df.sort_values(by='CNN_Affinity', ascending=False).reset_index(drop=True)
    vina_df = df.sort_values(by='Vina_Affinity', ascending=True).reset_index(drop=True)
    
    num_mols = len(df)
    
    # 3. Build the highly specific 12-column block for this run
    run_block = pd.DataFrame()
    
    run_block['Run_Name'] = [run_name] * num_mols
    run_block['Ligand_by_Consensus'] = cons_df['Ligand']
    run_block['Consensus_Rank'] = cons_df['Consensus_Rank']
    
    run_block['Spacer_1'] = ""  # Column 4 (Empty)
    
    run_block['Ligand_by_Pose'] = pose_df['Ligand']
    run_block['CNN_Pose_Score'] = pose_df['CNN_Pose_Score']
    
    run_block['Spacer_2'] = ""  # Column 7 (Empty)
    
    run_block['Ligand_by_Affinity'] = aff_df['Ligand']
    run_block['CNN_Affinity'] = aff_df['CNN_Affinity']
    
    run_block['Spacer_3'] = ""  # Column 10 (Empty)
    
    run_block['Ligand_by_Vina'] = vina_df['Ligand']
    run_block['Vina_Affinity'] = vina_df['Vina_Affinity']
    
    # 4. Create the blank row and glue it to the bottom of the block
    blank_row = pd.DataFrame([[""] * len(run_block.columns)], columns=run_block.columns)
    run_block = pd.concat([run_block, blank_row], ignore_index=True)
    
    # Add this completed block to our master list
    all_runs_data.append(run_block)

# 5. Glue all the blocks together
# Because we sorted tsv_files at the top, we DO NOT need to sort the master_df here. 
# If we sorted here, it would push all our carefully placed blank rows to the bottom!
master_df = pd.concat(all_runs_data, ignore_index=True)

# 6. Export to CSV using QUOTE_ALL to prevent Excel Power Query from mangling text
output_filename = "detailed_docking_scores_blocked.csv"
master_df.to_csv(output_filename, index=False, quoting=csv.QUOTE_ALL)

print(f"Success! Exported detailed blocks to {output_filename}")
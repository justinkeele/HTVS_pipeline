import pandas as pd
import glob
import os
import csv

print("Compiling docking data...")

# 1. Find all the TSV summary files across all your run folders
tsv_files = glob.glob("*/all_ligands_master_summary.tsv")

all_runs_data = []

for file in tsv_files:
    # Extract the name of the folder (e.g., "7B7D_HEAT_ex32_rescore_GPU")
    run_name = os.path.dirname(file)
    
    # Read the TSV file
    df = pd.read_csv(file, sep='\t', names=['Ligand', 'Vina_Affinity', 'CNN_Pose_Score', 'CNN_Affinity'])
    df = df.dropna()
    
    # Calculate the Consensus Ranks
    df['Pose_Rank'] = df['CNN_Pose_Score'].rank(ascending=False, method='min')
    df['Affinity_Rank'] = df['CNN_Affinity'].rank(ascending=False, method='min')
    df['Consensus_Rank'] = (df['Pose_Rank'] + df['Affinity_Rank']) / 2.0

    # 2. CREATE THE SORTED LEADERBOARDS
    
    # Sort by Consensus (Ascending: 1 is best). Tie-breaker is Pose Score.
    consensus_sorted = df.sort_values(by=['Consensus_Rank', 'CNN_Pose_Score'], ascending=[True, False])['Ligand'].tolist()
    
    # Sort by CNN Pose Score (Descending: 1.0 is best)
    pose_sorted = df.sort_values(by='CNN_Pose_Score', ascending=False)['Ligand'].tolist()
    
    # Sort by CNN Affinity (Descending: highest pK is best)
    aff_sorted = df.sort_values(by='CNN_Affinity', ascending=False)['Ligand'].tolist()
    
    # Sort by Vina Affinity (Ascending: most negative kcal/mol is best)
    vina_sorted = df.sort_values(by='Vina_Affinity', ascending=True)['Ligand'].tolist()
    
    num_mols = len(consensus_sorted)

    # 3. BUILD THE RUN'S DATAFRAME
    # Instead of listing scores, we list the names of the ligands in their sorted order
    run_df = pd.DataFrame({
        'Run_Name': [run_name] * num_mols,
        'Rank': list(range(1, num_mols + 1)),
        'Sorted_by_Consensus': consensus_sorted,
        'Sorted_by_Pose_Score': pose_sorted,
        'Sorted_by_CNN_Affinity': aff_sorted,
        'Sorted_by_Vina_Affinity': vina_sorted
    })
    
    all_runs_data.append(run_df)

# 4. Combine all runs into the master spreadsheet
master_df = pd.concat(all_runs_data, ignore_index=True)

# Sort the spreadsheet so runs are grouped alphabetically, reading Rank 1 to 24
master_df = master_df.sort_values(by=['Run_Name', 'Rank'])

# 5. Export to CSV using QUOTE_ALL to prevent Excel from mangling the text
output_filename = "hyperparameter_sorting.csv"
master_df.to_csv(output_filename, index=False)

print(f"Success! Exported {len(tsv_files)} runs to {output_filename}")
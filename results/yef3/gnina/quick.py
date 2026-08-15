import pandas as pd

print("Processing telemetry and aligning sort order (keeping duplicates)...")

# 1. Load the raw telemetry logs and combine them
tel1 = pd.read_csv('grid_search_telemetry.csv')
tel2 = pd.read_csv('grid_search_telemetry2.csv')
tel_all = pd.concat([tel1, tel2], ignore_index=True)

# 2. Clean the column names to remove any accidental trailing spaces
tel_all.columns = tel_all.columns.str.strip()

# 3. Drop the broken/empty column
if 'Init_Time_to_start' in tel_all.columns:
    tel_all = tel_all.drop(columns=['Init_Time_to_start'])

# 4. Filter out 0s in time, mols, and CPU columns
tel_filtered = tel_all[
    (tel_all['Mols'] > 0) & 
    (tel_all['Avg_Time_Per_Mol_s'] > 0) & 
    (tel_all['Peak_CPU_Cores'] > 0) & 
    (tel_all['Avg_CPU_Cores'] > 0)
].copy()

# WE ARE NO LONGER DROPPING DUPLICATES. They will be listed sequentially.

# 5. Load the human-readable spreadsheet to get the desired sort order
human_scores = pd.read_csv('human_readable_scores.csv')
ordered_runs = human_scores['Run_Name'].dropna().unique().tolist()

# 6. Rescue the missing runs (like the CPU runs)
# Find any runs in our telemetry that are NOT in the human_scores order
extra_runs = sorted(list(set(tel_filtered['Run_Name']) - set(ordered_runs)))

# Stitch the two lists together: the human_scores order first, followed by the extra runs at the bottom
full_ordered_list = ordered_runs + extra_runs

# 7. Apply the custom sort order to the telemetry data
tel_filtered['Run_Name'] = pd.Categorical(tel_filtered['Run_Name'], categories=full_ordered_list, ordered=True)

# Sort the dataframe. Duplicates will naturally sit right next to each other.
tel_filtered = tel_filtered.sort_values('Run_Name').reset_index(drop=True)

# 8. Output the final aligned telemetry file
output_file = 'master_telemetry.csv'
tel_filtered.to_csv(output_file, index=False)

print(f"Successfully processed {len(tel_filtered)} rows.")
print(f"Saved to {output_file}. Please review manually for duplicates.")
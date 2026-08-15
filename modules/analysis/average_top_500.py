import pandas as pd

# ==========================================
# STEP 0: Define File Paths
# ==========================================
# Define where the input files live relative to ~/Jagdish_lab/YeF3/results/compare_tophits/
SMILES_FILE = "all_deduplicated_mols.smi"
FILE_2IW3   = "../HTVS_runs/array_2IW3_all_mols/mega_sorted_docking_scores.txt"
FILE_FAR    = "../HTVS_runs/array_7B7D_far_all_mols/mega_sorted_docking_scores.txt"
FILE_HEAT   = "../HTVS_runs/array_7B7D_HEAT_all_mols/mega_sorted_docking_scores.txt"

OUTPUT_CSV  = "top500_7B7D_average.csv"

# ==========================================
# STEP 1: Capture Individual Top 500 Sets
# ==========================================
# We read only the first 500 rows (nrows=500) of each sorted docking score file.
# Using set(...) allows us to instantly check membership later using set math.

def get_top500_set(filepath):
    # Load just the first column (ZINC_ID) for the top 500 rows
    df_top = pd.read_csv(filepath, sep=r'\s+', header=None, nrows=500, usecols=[0], names=['ZINC_ID'])
    return set(df_top['ZINC_ID'])

print("Extracting top 500 hit sets from individual runs...")
top500_2IW3 = get_top500_set(FILE_2IW3)
top500_far  = get_top500_set(FILE_FAR)
top500_HEAT = get_top500_set(FILE_HEAT)

# ==========================================
# STEP 2: Load the Complete Datasets
# ==========================================
# sep=r'\s+' tells pandas to split by any amount of whitespace (spaces or tabs).
print("Loading all 310,000 molecules from each file...")

mols_df = pd.read_csv(SMILES_FILE, sep=r'\s+', header=None, names=['SMILES', 'ZINC_ID'])
df_2IW3 = pd.read_csv(FILE_2IW3,   sep=r'\s+', header=None, names=['ZINC_ID', 'Score_2IW3'])
df_far  = pd.read_csv(FILE_FAR,    sep=r'\s+', header=None, names=['ZINC_ID', 'Score_7B7D_far'])
df_HEAT = pd.read_csv(FILE_HEAT,   sep=r'\s+', header=None, names=['ZINC_ID', 'Score_7B7D_HEAT'])

# ==========================================
# STEP 3: Merge Everything on ZINC_ID
# ==========================================
# We chain merges together so all data lines up by ZINC_ID into one master table.
print("Merging datasets into a master table...")

master_df = mols_df.merge(df_far,  on='ZINC_ID', how='inner') \
                   .merge(df_HEAT, on='ZINC_ID', how='inner') \
                   .merge(df_2IW3, on='ZINC_ID', how='inner')

# ==========================================
# STEP 4: Calculate the Averages (Columns 6 & 7)
# ==========================================
# Perform vector math across the entire dataframe at once.
master_df['Avg_7B7D'] = (master_df['Score_7B7D_far'] + master_df['Score_7B7D_HEAT']) / 2.0
master_df['Avg_All3'] = (master_df['Score_2IW3'] + master_df['Score_7B7D_far'] + master_df['Score_7B7D_HEAT']) / 3.0

# round averages to 3 decimal places for a clean presentation
master_df['Avg_7B7D'] = master_df['Avg_7B7D'].round(3)
master_df['Avg_All3'] = master_df['Avg_All3'].round(3)

# ==========================================
# STEP 5: Sort and Slice the Top 500
# ==========================================
# Sort ascending because more negative Vina scores = stronger predicted binding.
# .head(500) grabs exactly the top 500 rows after sorting.
print("Sorting by 7B7D average and selecting the top 500 hits...")

top_500_df = master_df.sort_values(by='Avg_7B7D', ascending=True).head(500).copy()

# ==========================================
# STEP 6: Highlight Multi-Pocket Overlaps (Columns 8 & 9)
# ==========================================
# The '&' operator between sets finds the INTERSECTION (molecules present in BOTH sets).
# .isin() checks if each ZINC_ID is in that intersection set.
# .map({True: 'Yes', False: 'No'}) converts True/False into clean Yes/No labels.

both_7b7d_set = top500_far & top500_HEAT
all_three_set = top500_far & top500_HEAT & top500_2IW3

top_500_df['In_Top500_Both_7B7D'] = top_500_df['ZINC_ID'].isin(both_7b7d_set).map({True: 'Yes', False: '   '})
top_500_df['In_Top500_All_3']     = top_500_df['ZINC_ID'].isin(all_three_set).map({True: 'Yes', False: '   '})

# ==========================================
# STEP 7: Reorder Columns and Export
# ==========================================
# Enforce the exact 9-column order requested for your presentation.
final_columns = [
    'ZINC_ID',              # Col 1
    'SMILES',               # Col 2
    'Score_2IW3',           # Col 3
    'Score_7B7D_far',       # Col 4
    'Score_7B7D_HEAT',      # Col 5
    'Avg_7B7D',             # Col 6 (Sorted by this)
    'Avg_All3',             # Col 7 (Unsorted)
    'In_Top500_Both_7B7D',  # Col 8
    'In_Top500_All_3'       # Col 9
]

top_500_df = top_500_df[final_columns]

# Export to CSV without the pandas numerical index column
top_500_df.to_csv(OUTPUT_CSV, index=False)
print(f"Success! Presentation file saved to: {OUTPUT_CSV}")
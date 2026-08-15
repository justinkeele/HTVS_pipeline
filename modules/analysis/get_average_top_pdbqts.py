import os
import glob
import shutil
import pandas as pd

# ==========================================
# STEP 0: Define Inputs and Target Directories
# ==========================================
CSV_FILE = "top500_7B7D_average.csv"

# Map each friendly pocket name to its source array directory
POCKET_DIRS = {
    "2IW3":      "../HTVS_runs/array_2IW3_all_mols",
    "7B7D_far":  "../HTVS_runs/array_7B7D_far_all_mols",
    "7B7D_HEAT": "../HTVS_runs/array_7B7D_HEAT_all_mols"
}

# ==========================================
# STEP 1: Load the Top 10 and Top 50 IDs
# ==========================================
print(f"Reading top hits from {CSV_FILE}...")
df = pd.read_csv(CSV_FILE)

top_10_ids = df['ZINC_ID'].head(10).tolist()
top_50_ids = df['ZINC_ID'].head(50).tolist()

# Organize into a dictionary so we can loop through both cutoffs easily
CUTOFFS = {
    "top_10": top_10_ids,
    "top_50": top_50_ids
}

# ==========================================
# STEP 2: Find and Copy the PDBQT Files
# ==========================================
def copy_pdbqts_for_group(pocket_name, source_base_dir, cutoff_name, zinc_ids):
    # Create the output directory (e.g., "7B7D_HEAT_top_10")
    dest_dir = f"{pocket_name}_{cutoff_name}"
    os.makedirs(dest_dir, exist_ok=True)
    
    copied_count = 0
    missing_count = 0

    for zinc_id in zinc_ids:
        # Use glob wildcard (*) to check batch_1 through batch_32 instantly
        pattern = f"{source_base_dir}/batch_*/{zinc_id}.pdbqt"
        matches = glob.glob(pattern)
        
        if matches:
            # Grab the first match found and copy it to dest_dir
            source_file = matches[0]
            shutil.copy2(source_file, dest_dir)
            copied_count += 1
        else:
            print(f"  [Warning] Could not find {zinc_id}.pdbqt in {pocket_name}")
            missing_count += 1

    print(f"[{dest_dir}] Copied {copied_count} files ({missing_count} missing).")

# ==========================================
# STEP 3: Run Across All Pockets and Cutoffs
# ==========================================
print("\nStarting PDBQT copy workflow...\n" + "="*40)

for pocket_name, source_dir in POCKET_DIRS.items():
    for cutoff_name, id_list in CUTOFFS.items():
        copy_pdbqts_for_group(pocket_name, source_dir, cutoff_name, id_list)

print("="*40 + "\nAll top 10 and top 50 PDBQT folders are ready!")



# ==========================================
# STEP 1: Load the Ranked ZINC IDs from CSV
# ==========================================
# Use the CSV that contains your sorted top 500 hits
CSV_FILE = "top500_7B7D_average.csv"  # Change to "top500_7B7D_average_presentation.csv" if needed

print(f"Loading rank order from {CSV_FILE}...")
df = pd.read_csv(CSV_FILE)

# Create a dictionary mapping ZINC_ID -> Rank (e.g., {'ZINCou0000009b47': 1, ...})
# enumerate(..., start=1) ensures rank 1 is 1, not 0.
rank_map = {zinc_id: rank for rank, zinc_id in enumerate(df['ZINC_ID'], start=1)}

# ==========================================
# STEP 2: Find All top_10 / top_50 Folders
# ==========================================
# glob("*_top_*") automatically finds 2IW3_top_50, 7B7D_HEAT_top_50, etc.
target_dirs = glob.glob("*_top_*")

for folder in target_dirs:
    if not os.path.isdir(folder):
        continue
    
    print(f"\nOrdering files in: {folder}")
    renamed_count = 0
    
    # Check every .pdbqt file currently in the folder
    for filepath in glob.glob(f"{folder}/*.pdbqt"):
        filename = os.path.basename(filepath)
        
        # Extract the raw ZINC_ID by stripping out the extension
        # (This also works if you re-run the script on already-prefixed files)
        raw_name = filename.replace(".pdbqt", "")
        zinc_id = raw_name.split("_")[-1]  # Grabs ZINC ID even if filename is already "01_ZINC..."
        
        if zinc_id in rank_map:
            rank = rank_map[zinc_id]
            
            # f"{rank:02d}" formats the number with a leading zero (01, 02, ... 09, 10)
            new_filename = f"{rank:02d}_{zinc_id}.pdbqt"
            new_filepath = os.path.join(folder, new_filename)
            
            # Rename the file in-place
            if filepath != new_filepath:
                os.rename(filepath, new_filepath)
                renamed_count += 1
        else:
            print(f"  [Warning] {zinc_id} not found in CSV rank map.")

    print(f"  -> Successfully ordered {renamed_count} files.")

print("\nDone! All PDBQT files now sort by rank automatically.")
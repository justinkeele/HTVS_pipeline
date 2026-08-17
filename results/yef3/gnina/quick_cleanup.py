import os
import glob
import pandas as pd

# Define Output Path
OUTPUT_DIR = "../../../docs/graphs_with_good_paths"

def cleanup_telemetry():
    print("Finding legacy telemetry CSV files...")
    
    # 1. Grab any file that looks like our old telemetry CSVs
    # This securely catches grid_search_telemetry.csv and grid_search_telemetry_2.csv
    csv_files = glob.glob("*telemetry*.csv")
    
    # Exclude 'hardware_telemetry.csv' so we don't accidentally ingest our own output
    csv_files = [f for f in csv_files if "hardware" not in f]
    
    if not csv_files:
        print("Error: Could not find any telemetry CSV files to clean.")
        return
        
    all_data = []
    
    for file in csv_files:
        print(f"Loading {file}...")
        df = pd.read_csv(file)
        # Strip trailing whitespaces from headers just in case
        df.columns = df.columns.str.strip()
        all_data.append(df)
        
    # 2. Combine into a single master dataframe
    master_df = pd.concat(all_data, ignore_index=True)
    
    # 3. Rename the GPU column to match our strict format
    if 'Avg_GPU_Util' in master_df.columns:
        master_df = master_df.rename(columns={'Avg_GPU_Util': 'Avg_GPU_Util_Pct'})
        
    # 4. THE FILTERING GATE (Clean out the junk)
    # We strictly enforce the rules: 24 molecules, no 0 cores, no instant crashes
    clean_df = master_df[
        (master_df['Mols'] == 24) &
        (master_df['Avg_Time_Per_Mol_s'] > 10.0) &
        (master_df['Peak_CPU_Cores'] > 0) &
        (master_df['Avg_CPU_Cores'] > 0) &
        (master_df['Peak_RAM_MB'] > 0) &
        # NEW RULE: If the Run_Name contains 'GPU', the Avg_GPU_Util_Pct MUST be > 0.
        # (~ means NOT). So: It is NOT a GPU run, OR the GPU utilization is > 0.
        (~master_df['Run_Name'].str.contains('GPU') | (master_df['Avg_GPU_Util_Pct'] > 0))
    ].copy()
    
    # 5. Extract strictly the requested columns
    target_cols = [
        'Run_Name', 'Avg_Time_Per_Mol_s', 'Peak_CPU_Cores', 
        'Avg_CPU_Cores', 'Peak_RAM_MB', 'Peak_VRAM_MB', 'Avg_GPU_Util_Pct'
    ]
    
    # Ensure all target columns exist (fill with 0.0 if a column was missing)
    for col in target_cols:
        if col not in clean_df.columns:
            clean_df[col] = 0.0
            
    final_df = clean_df[target_cols].copy()
    
    # 6. DUPLICATE HANDLING (Preserving them so you can manually review)
    # Group by Run_Name and assign a sequence number (0, 1, 2...)
    final_df['Run_Name'] = final_df['Run_Name'] + '_' + final_df.groupby('Run_Name').cumcount().astype(str)
    
    # Clean up the names: Drop the _0 for the first run, change _1 to _run2, _2 to _run3
    final_df['Run_Name'] = final_df['Run_Name'].str.replace(r'_0$', '', regex=True)
    final_df['Run_Name'] = final_df['Run_Name'].str.replace(r'_(\d+)$', lambda m: f"_run{int(m.group(1))+1}", regex=True)
    
    # 7. EXPORT
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    output_file = os.path.join(OUTPUT_DIR, "hardware_telemetry_cleanup.csv")
    
    final_df.to_csv(output_file, index=False)
    
    print(f"\nSuccessfully cleaned out the junk.")
    print(f"Total valid runs preserved: {len(final_df)}")
    print(f"Exported perfectly formatted data to {output_file}")

if __name__ == '__main__':
    cleanup_telemetry()
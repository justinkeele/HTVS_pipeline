import os
import glob
import pandas as pd
import re

# 1. Define your relative paths based on where this script lives
DATA_DIR = "../../../results/yef3/gnina"
OUTPUT_DIR = "../../../docs/graphs_with_good_paths"

def parse_run_name(run_name):
    '''
    Breaks down a folder name like "7B7D_HEAT_ex32_rescore_GPU_repeat_6"
    into its individual metadata components.
    '''
    # --- CHANGE 1: THE FILTERING GATE ---
    # If the folder doesn't have "repeat_" in the name, abort and return None.
    # This prevents your old "test_4" runs from infecting the new data.
    if "repeat_" not in run_name:
        return None

    pocket = "Unknown"
    exhaustiveness = "Unknown"
    mode = "Unknown"
    hardware = "Unknown"
    repeat_number = "Unknown" 

    parts = run_name.split('_')
    
    if "2IW3" in run_name:
        pocket = "2IW3"
    elif "HEAT" in run_name:
        pocket = "7B7D_HEAT"
    elif "far" in run_name:
        pocket = "7B7D_far"
        
    for p in parts:
        if p.startswith("ex") and p[2:].isdigit():
            exhaustiveness = p
            break
            
    if "rescore" in run_name:
        mode = "rescore"
    elif "refinement" in run_name:
        mode = "refinement"
        
    if "GPU" in run_name:
        hardware = "GPU"
    elif "CPU" in run_name:
        hardware = "CPU"
        
    # --- CHANGE 2: DYNAMIC REGEX EXTRACTION ---
    # Automatically extracts the integer following "repeat_"
    match = re.search(r'repeat_(\d+)', run_name)
    if match:
        repeat_number = match.group(1)

    return pocket, exhaustiveness, mode, hardware, repeat_number

def process_file(target_filename, output_csv):
    '''
    Hunts down every instance of target_filename in the DATA_DIR, 
    extracts the scores, applies the metadata, and saves to OUTPUT_DIR.
    '''
    all_data = []
    
    search_pattern = os.path.join(DATA_DIR, "*", target_filename)
    files_found = glob.glob(search_pattern)
    
    for file_path in files_found:
        full_dir_path = os.path.dirname(file_path)
        run_name = os.path.basename(full_dir_path)
        
        parsed = parse_run_name(run_name)
        
        # If parse_run_name returned None (because it's an old test run), skip it!
        if parsed is None:
            continue
            
        pocket, ex, mode, hw, repeat_num = parsed
        
        try:
            df = pd.read_csv(file_path, sep='\t', header=None, names=['Ligand', 'Score'])
            
            df['Run_Name'] = run_name
            df['Pocket'] = pocket
            df['Exhaustiveness'] = ex
            df['Mode'] = mode
            df['Hardware'] = hw
            
            # --- CHANGE 3: NEW COLUMN NAME ---
            df['Repeat_Number'] = repeat_num
            
            cols = ['Ligand', 'Run_Name', 'Pocket', 'Exhaustiveness', 'Mode', 'Hardware', 'Repeat_Number', 'Score']
            df = df[cols]
            
            all_data.append(df)
            
        except Exception as e:
            print(f"Skipped {file_path} due to error: {e}")
            
    if all_data:
        master_df = pd.concat(all_data, ignore_index=True)
        
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        
        final_output_path = os.path.join(OUTPUT_DIR, output_csv)
        master_df.to_csv(final_output_path, index=False)
        print(f"Successfully created {final_output_path} with {len(master_df)} rows.")
    else:
        print(f"No data found for {target_filename} matching the 'repeat_' criteria.")

if __name__ == '__main__':
    print("Starting data extraction for 10-run replicates...")
    # --- CHANGE 4: NEW TARGET OUTPUT NAMES ---
    process_file("sorted_cnn_affinity.txt", "repeat_runs_cnn_affinity.csv")
    process_file("sorted_cnn_pose_score.txt", "repeat_runs_cnn_pose_score.csv")
    process_file("sorted_vina_affinity.txt", "repeat_runs_vina_affinity.csv")
    print("Extraction complete!")
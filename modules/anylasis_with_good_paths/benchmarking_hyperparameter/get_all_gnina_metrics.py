import os
import glob
import pandas as pd

# 1. Define your relative paths based on where this script lives
DATA_DIR = "../../../results/yef3/gnina"
OUTPUT_DIR = "../../../docs/graphs_with_good_paths"

def parse_run_name(run_name):
    '''
    Breaks down a folder name like "7B7D_far_ex128_rescore_GPU_test_3"
    into its individual metadata components based on your naming convention.
    '''
    pocket = "Unknown"
    exhaustiveness = "Unknown"
    mode = "Unknown"
    hardware = "Unknown"
    test_version = "1" 

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
        
    if "test_2" in run_name:
        test_version = "2"
    elif "test_3" in run_name:
        test_version = "3"
    elif "test_4" in run_name:
        test_version = "4"
    else:
        test_version = "1" 

    return pocket, exhaustiveness, mode, hardware, test_version

def process_file(target_filename, output_csv):
    '''
    Hunts down every instance of target_filename in the DATA_DIR, 
    extracts the scores, applies the metadata, and saves to OUTPUT_DIR.
    '''
    all_data = []
    
    # 2. Build the search path to look inside the data directory
    search_pattern = os.path.join(DATA_DIR, "*", target_filename)
    files_found = glob.glob(search_pattern)
    
    for file_path in files_found:
        # 3. Extract JUST the folder name, ignoring the "../../../" pathing
        # os.path.dirname gets the full folder path
        # os.path.basename isolates just the final folder name
        full_dir_path = os.path.dirname(file_path)
        run_name = os.path.basename(full_dir_path)
        
        pocket, ex, mode, hw, test_ver = parse_run_name(run_name)
        
        try:
            df = pd.read_csv(file_path, sep='\t', header=None, names=['Ligand', 'Score'])
            
            df['Run_Name'] = run_name
            df['Pocket'] = pocket
            df['Exhaustiveness'] = ex
            df['Mode'] = mode
            df['Hardware'] = hw
            df['Test_Version'] = test_ver
            
            cols = ['Ligand', 'Run_Name', 'Pocket', 'Exhaustiveness', 'Mode', 'Hardware', 'Test_Version', 'Score']
            df = df[cols]
            
            all_data.append(df)
            
        except Exception as e:
            print(f"Skipped {file_path} due to error: {e}")
            
    if all_data:
        master_df = pd.concat(all_data, ignore_index=True)
        
        # 4. Ensure the output directory actually exists before trying to save
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        
        # 5. Build the final output path and save
        final_output_path = os.path.join(OUTPUT_DIR, output_csv)
        master_df.to_csv(final_output_path, index=False)
        print(f"Successfully created {final_output_path} with {len(master_df)} rows.")
    else:
        print(f"No data found for {target_filename} in {DATA_DIR}")

if __name__ == '__main__':
    print("Starting data extraction...")
    process_file("sorted_cnn_affinity.txt", "master_cnn_affinity.csv")
    process_file("sorted_cnn_pose_score.txt", "master_cnn_pose_score.csv")
    process_file("sorted_vina_affinity.txt", "master_vina_affinity.csv")
    print("Extraction complete!")
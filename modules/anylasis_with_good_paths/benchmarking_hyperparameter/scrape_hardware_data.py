import os
import glob
import re
import pandas as pd

INPUT_DIR = "../../gnina"
OUTPUT_DIR = "../../../docs/graphs_with_good_paths"

def scrape_telemetry():
    print("Scanning log files for hardware telemetry...")
    

    # Find all log files (ONLY .out) in the specified input directory
    # CHANGED: Use os.path.join to target INPUT_DIR and removed the .txt search
    search_pattern = os.path.join(INPUT_DIR, "*.out")
    log_files = glob.glob(search_pattern)
    
    all_runs = []
    
    for file_path in log_files:
        with open(file_path, 'r') as file:
            content = file.read()
            
            # 2. THE REGEX ENGINE
            # We use re.finditer to find EVERY telemetry block in the file.
            # The regex looks for the exact start and end borders of your ASCII tables.
            table_pattern = re.compile(r'\|--.*?\|(.*?)\|--.*?\|', re.DOTALL)
            
            # This splits the file into chunks based on the "Run Name:" header
            # We use a trick to split the text every time a new run block starts
            blocks = content.split('| Run Name:')
            
            for block in blocks[1:]:  # Skip the first chunk before the first table
                try:
                    # Re-attach the text we split on so the regex can find it
                    block = '| Run Name:' + block
                    
                    # Extract the raw numbers using specific regex patterns
                    # (\S+) grabs words without spaces, ([\d.]+) grabs numbers with decimals
                    run_name = re.search(r'Run Name:\s+(\S+)', block).group(1)
                    mols = float(re.search(r'Total Mols:\s+(\d+)', block).group(1))
                    
                    # Some overheads are 0, which is normal.
                    overhead = float(re.search(r'Init. Overhead Time:\s+([\d.]+)', block).group(1))
                    avg_time = float(re.search(r'Avg Time / Molecule:\s+([\d.]+)', block).group(1))
                    peak_cpu = float(re.search(r'Peak CPU Cores Used:\s+([\d.]+)', block).group(1))
                    avg_cpu = float(re.search(r'Avg CPU Cores Used:\s+([\d.]+)', block).group(1))
                    peak_ram = float(re.search(r'Peak System RAM:\s+([\d.]+)', block).group(1))
                    peak_vram = float(re.search(r'Peak GPU VRAM:\s+([\d.]+)', block).group(1))
                    avg_gpu = float(re.search(r'Avg GPU Utilization:\s+([\d.]+)', block).group(1))
                    
                    # 3. THE VALIDATION GATE
                    # Reject failed runs based on your exact criteria
                    if mols != 24:
                        continue
                    if avg_time <= 10.0:
                        continue
                    if peak_cpu == 0 or avg_cpu == 0 or peak_ram == 0:
                        continue
                        
                    # If it survives the gauntlet, save the data
                    run_data = {
                        'Run_Name': run_name,
                        'Avg_Time_Per_Mol_s': avg_time,
                        'Peak_CPU_Cores': peak_cpu,
                        'Avg_CPU_Cores': avg_cpu,
                        'Peak_RAM_MB': peak_ram,
                        'Peak_VRAM_MB': peak_vram,
                        'Avg_GPU_Util_Pct': avg_gpu
                    }
                    all_runs.append(run_data)
                    
                except AttributeError:
                    # If the regex fails to find a number (because the table crashed halfway), skip it
                    pass

    # 4. DUPLICATE HANDLING & DATAFRAME CREATION
    if not all_runs:
        print("No valid telemetry data found.")
        return
        
    df = pd.DataFrame(all_runs)
    
    # Check for duplicate Run_Names. If found, append an index (_run0, _run1)
    # cumcount() counts how many times it has seen that specific Run_Name before
    df['Run_Name'] = df['Run_Name'] + '_' + df.groupby('Run_Name').cumcount().astype(str)

    # CSV EXPORT
    # CHANGED: Ensure the output directory exists and build the full path
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    final_output_path = os.path.join(OUTPUT_DIR, "hardware_telemetry.csv")
    df.to_csv(final_output_path, index=False)
    print(f"Successfully scraped {len(df)} valid runs.")
    print(f"Exported clean data to {final_output_path}")

if __name__ == '__main__':
    scrape_telemetry()
#!/bin/bash
# Enable strict error handling
set -e

#=====================================================
# Global Variables (edit when cloning to new machine) It points to the root folder HTVS_pipeline
#=====================================================
export WORKSPACE_DIR="/home/justin/Jagdish_lab/HTVS_pipeline"

time_begin=$(date +%s)

# Dynamic Conda Auto-Detector
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
else
    echo "Error: Conda not found. Please install Conda to a path like: \n$HOME/miniconda3/etc/profile.d/conda.sh or $HOME/anaconda3/etc/profile.d/conda.sh and try again."
    exit 1
fi
conda activate gnina
# Export the active Conda PATH so GNU Parallel subshells inherit it natively
export PATH

#=====================================================
# Enviroment detector
#=====================================================

# Anchor the script to the directory it was submitted from
if [[ -n "$SLURM_SUBMIT_DIR" ]]; then
    cd "$SLURM_SUBMIT_DIR" || exit
else
    LAUNCH_DIR=$(pwd)
    cd "$(dirname "$0")" || exit
fi

if [[ -n "$SLURM_CPUS_PER_TASK" ]]; then
    cpu_cores=$SLURM_CPUS_PER_TASK
    echo "SLURM detected, runing with $cpu_cores cores."
    # Force GNU Parallel to use Bash
    export PARALLEL_SHELL=/bin/bash
    export bash_shells="$SLURM_NTASKS"
else
    cpu_cores=20
    echo "No SLURM detected, running with $cpu_cores cores."
    export bash_shells=4
    # We will add a safety check here to make sure that the script dosen't ask for more cores than avilable. 
fi

# =====================================================
# Argument Capture & Validation
# =====================================================

# 1. Check if the user provided the correct number of arguments
if [ "$#" -lt 6 ]; then
    echo "=========================================================="
    echo " ERROR: Missing arguments"
    echo " USAGE: ./GNINA_worker.sh <PROTIEN> <TARGET_ID> <EXHAUSTIVENESS> <CNN_SCORING> <GPU or CPU> <LIGAND_BATCH>"
    echo " EXAMPLE: ./GNINA_worker.sh yef3 2IW3 16 rescore GPU test_batch"
    echo "=========================================================="
    exit 1
fi

# Capture the command-line arguments
TARGET_ID="$1"
POCKET="$2"
EXHAUSTIVENESS="$3"
CNN_SCORING="$4"
HARDWARE_TYPE="$5"
LIGAND_BATCH="$6"
#=====================================================
# Validate Exhaustiveness 
#=====================================================

# Check if EXHAUSTIVENESS is a whole positive number
if ! [[ "$EXHAUSTIVENESS" =~ ^[0-9]+$ ]]; then
    echo "ERROR: EXHAUSTIVENESS must be a whole positive number (e.g., 8, 16, 32)."
    exit 1
fi

# Warn if exaustivness is outside recommended range
if [ "$EXHAUSTIVENESS" -lt 8 ] || [ "$EXHAUSTIVENESS" -gt 32 ]; then
    echo "WARNING: EXHAUSTIVENESS is set to $EXHAUSTIVENESS."
    echo "Values below 8 may result in poor poses. Values above 32 are computationally expensive with diminishing returns."
    sleep 4 # Pause briefly so the user sees the warning
fi

# =====================================================
# Dependency Check
# =====================================================
for req_cmd in gnina mk_prepare_ligand.py parallel; do
    if ! command -v "$req_cmd" &> /dev/null; then
        echo "FATAL ERROR: Required command '$req_cmd' is not installed or not in PATH."
        echo "Make sure your Conda environment is active and Apptainer wrappers are configured. The required dependecies _______ are missing."
        exit 1
    fi
done

#=====================================================
# Path Setup
#=====================================================
export batches_dir_path="${WORKSPACE_DIR}/targets/${TARGET_ID}/ligand_batches/${LIGAND_BATCH}/smiles_pH_H_gen3D_output"


# Note: Right now this only accepts 7B7D_HEAT, 7B7D_far, or 2IW3. 
# We need to change this later to dynamically identify pockets and configs.
case "$POCKET" in
    "7B7D_HEAT")
        export receptor_path="${WORKSPACE_DIR}/targets/${TARGET_ID}/receptors/7B7D_72pH.pdbqt"
        export config_gnina_path="${WORKSPACE_DIR}/targets/${TARGET_ID}/configs/gnina/7B7D_HEAT_config.txt"
        ;;
    "7B7D_far")
        export receptor_path="${WORKSPACE_DIR}/targets/${TARGET_ID}/receptors/7B7D_72pH.pdbqt"
        export config_gnina_path="${WORKSPACE_DIR}/targets/${TARGET_ID}/configs/gnina/7B7D_far_config.txt"
        ;;
    "2IW3")
        export receptor_path="${WORKSPACE_DIR}/targets/${TARGET_ID}/receptors/2IW3_47pH.pdbqt"
        export config_gnina_path="${WORKSPACE_DIR}/targets/${TARGET_ID}/configs/gnina/2IW3_config.txt"
        ;;
    *)
        echo "ERROR: Unknown pocket '$POCKET'. Must be 7B7D_HEAT, 7B7D_far, or 2IW3."
        exit 1
        ;;
esac
#=====================================================
# Run names and directory setup
#=====================================================

export run_name="${POCKET}_ex${EXHAUSTIVENESS}_${CNN_SCORING}_${HARDWARE_TYPE}"
export output_dir="${WORKSPACE_DIR}/results/${TARGET_ID}/gnina/${run_name}"
export meeko_temp="${output_dir}/meeko_temp"

mkdir -p "$output_dir"
mkdir -p "$output_dir/results"
mkdir -p "$meeko_temp"

# =====================================================
# Exporting Variables to GNU Parallel for GNINA Docking
# =====================================================
export EXHAUSTIVENESS
export CNN_SCORING
export HARDWARE_TYPE


#=====================================================
#GNU function
#=====================================================

process_ligand() {
    full_path_to_ligand_var=$1   # $1 will be the full path caught by parallel
    clean_name=$(basename "$full_path_to_ligand_var" .pdbqt) #strips the dir off the name
    output_sdf="${output_dir}/${clean_name}_result.sdf"  #output file path
    
    temp_log="${output_dir}/${clean_name}_temp.log"
    final_log="${output_dir}/${clean_name}_results.log"

    # Sanitize garbage ZINC forcefield types (like CG0, CG1) back to standard Carbon (C) so GNINA quits throwing a fit
    sed -i -E 's/ C?G[0-9][[:space:]]*$/ C/g' "$full_path_to_ligand_var"

    echo "starting gnina for: $clean_name"
    
    # Toggles srun based on the environment
    if [[ -n "$SLURM_CPUS_PER_TASK" ]]; then
        gnina_cores="$SLURM_CPUS_PER_TASK"
        EXEC_CMD="srun --exclusive -N1 -n1 --cpus-per-task=$SLURM_CPUS_PER_TASK gnina"
    else
        gnina_cores=5
        EXEC_CMD="gnina"
    fi

    # Hardware toggle (CPU vs GPU)
    # The ^^ syntax forces the user input to uppercase so "cpu", "Cpu", and "CPU" all match perfectly.
    if [[ "${HARDWARE_TYPE^^}" == "CPU" ]]; then
        gpu_flag="--no_gpu"
    else
        gpu_flag=""
    fi

    # Actual gnina stuff
    $EXEC_CMD \
        -r "$receptor_path" \
        -l "$full_path_to_ligand_var" \
        --config "$config_gnina_path" \
        --exhaustiveness "$EXHAUSTIVENESS" \
        --scoring vinardo \
        --cnn_scoring "$CNN_SCORING" \
        $gpu_flag \
        -o "$output_sdf" \
        --cpu "$gnina_cores" \
        > "$temp_log" 2>&1
        
    awk '/mode | affinity/ {flag=1} /ERROR|Error|WARNING|Fatal/ {print} flag && NF > 0 {print} /^$/ && flag {flag=0}' "$temp_log" > "$final_log"

    if [ -s "$output_sdf" ]; then
        # comment this line out for extra debugging
        #rm -f "$temp_log"
        true
    else
        echo -e "${clean_name}\t failed gnina docking" >> "$batch_failed_log"
    fi

}
export -f process_ligand  #export the function to be used by GNU parallel

#=====================================================
# Hardware Monitor Daemon
#=====================================================
start_hardware_monitor() {
    local stat_file=$1
    local process_name=$2
    > "$stat_file" # Clear the file if it exists

    while true; do
        # 1. CPU & RAM (RSS in MB) isolated to your target process
        # Divides CPU% by 100 to get exact number of active cores
        read cpu_cores ram_mb <<< $(ps -eo %cpu,rss,command | grep "[$process_name:0:1]${process_name:1}" | awk '{cpu+=$1; ram+=$2} END {if(cpu=="") cpu=0; if(ram=="") ram=0; print cpu/100, ram/1024}')
        
        # 2. GPU metrics (Only runs if nvidia-smi is available on the node)
        if command -v nvidia-smi &> /dev/null; then
            read vram_mb gpu_util <<< $(nvidia-smi --query-gpu=memory.used,utilization.gpu --format=csv,noheader,nounits | awk '{vram+=$1; util+=$2} END {if(vram=="") vram=0; if(util=="") util=0; print vram, util}')
        else
            vram_mb=0
            gpu_util=0
        fi
        
        # Write to our hidden tracking file
        echo "$cpu_cores $ram_mb $vram_mb $gpu_util" >> "$stat_file"
        sleep 1
    done
}

echo "sacrificing ligands to the machine gods, this may take a while..."

for i in $(seq 1 1); do
    export batch_failed_log="$output_dir/results/failed_mols_batch_${i}.txt"
    export meeko_log="$output_dir/meeko_batch_${i}.log"
    echo "================================="
    echo "|       Starting batch $i       |"
    echo "================================="
    
    input_chunk="$batches_dir_path/input_batch_${i}/chunk_${i}.sdf"
    if [ ! -f "$input_chunk" ]; then
        echo "Error: $input_chunk not found."
        exit 1
    fi

    batch_begin=$(date +%s)
    
    echo "Exploding batch $i with meeko"
    mk_prepare_ligand.py -i "$input_chunk" --multimol_outdir "$meeko_temp" > "$meeko_log" 2>&1  || true
    grep -Ei "processed:|written:|error|duplicate|traceback|exception|fatal|warning" "$meeko_log" || true

    # Benchmarking triggers!
    docking_begin=$(date +%s)
    init_time=$((docking_begin - batch_begin))
    
    stat_file="$output_dir/hardware_stats_batch_${i}.tmp"
    
    # Launch the daemon in the background to watch 'gnina' and save its PID
    start_hardware_monitor "$stat_file" "gnina" &
    MONITOR_PID=$!

    echo "Docking batch $i in Gnina" 
    find "$meeko_temp" -name "*.pdbqt" | parallel --line-buffer --no-notice --joblog "$output_dir/results/batch_${i}_joblog.log" --resume -j "$bash_shells" process_ligand {} || true

    # Kill the background daemon the second GNU Parallel finishes
    kill $MONITOR_PID 2>/dev/null || true
    wait $MONITOR_PID 2>/dev/null || true


    # batch time
    batch_end=$(date +%s)
    batch_time=$((batch_end - batch_begin))
    batch_formatted=$(printf '%02d:%02d:%02d' $((batch_time/3600)) $((batch_time%3600/60)) $((batch_time%60)))
    echo "Batch $i completed in ${batch_formatted}"

    # Determine the environment for the output block
    if [[ -n "$SLURM_JOB_ID" ]]; then
        hw_env="This script is running on SLURM Node: $SLURMD_NODENAME (Allocated Cores: $SLURM_CPUS_PER_TASK)"
    else
        hw_env="This script is running on: $(hostname) (Allocated Cores: $cpu_cores)"
    fi

    # Extract Max and Average hardware stats
    read max_cpu avg_cpu max_ram max_vram avg_gpu <<< $(awk '
        BEGIN {max_c=0; sum_c=0; max_r=0; max_v=0; sum_u=0; count=0}
        {
            if ($1 > max_c) max_c = $1;
            sum_c += $1;                # Accumulate CPU cores
            if ($2 > max_r) max_r = $2;
            if ($3 > max_v) max_v = $3;
            sum_u += $4;
            count++;
        }
        END {
            if (count > 0) print max_c, sum_c/count, max_r, max_v, sum_u/count;
            else print 0, 0, 0, 0, 0;
        }' "$stat_file")

    # Get Average Time per Molecule from the GNU Parallel joblog
    joblog_file="$output_dir/results/batch_${i}_joblog.log"
    total_mols=$(awk 'NR>1 {count++} END {print count+0}' "$joblog_file")
    avg_mol_time=$(awk 'NR>1 {sum+=$4} END {if(NR>1) print sum/(NR-1); else print 0}' "$joblog_file")

    # Print the Pretty ASCII Benchmark Block
    echo " "
    echo "|----------------------------------------------------------------|"
    printf "| %-60s |\n" "Run Name:             $run_name"
    printf "| %-60s |\n" "Hardware Environment: $hw_env"
    echo "|----------------------------------------------------------------|"
    printf "| %-60s |\n" "Total Mols:           $total_mols"
    printf "| %-60s |\n" "$(printf 'Init. Overhead Time:  %s seconds' "$init_time")"
    printf "| %-60s |\n" "$(printf 'Avg Time / Molecule:  %.2f seconds' "$avg_mol_time")"
    echo "|----------------------------------------------------------------|"
    printf "| %-60s |\n" "$(printf 'Peak CPU Cores Used:  %.1f cores' "$max_cpu")"
    printf "| %-60s |\n" "$(printf 'Avg CPU Cores Used:   %.1f cores' "$avg_cpu")"
    printf "| %-60s |\n" "$(printf 'Peak System RAM:      %.1f MB' "$max_ram")"
    printf "| %-60s |\n" "$(printf 'Peak GPU VRAM:        %.0f MB' "$max_vram")"
    printf "| %-60s |\n" "$(printf 'Avg GPU Utilization:  %.1f%%' "$avg_gpu")"
    echo "|----------------------------------------------------------------|"
    echo " "

    # Clean up the tracker file
    rm -rf "$stat_file"
    rm -rf "${meeko_temp:?}/"*
done



summary_file="$output_dir/all_ligands_master_summary.tsv"
> "$summary_file"

# Single pass to extract $2 (Vina), $4 (CNNpose), and $5 (CNNaffinity) for Pose 1
find "$output_dir" -name "*_results.log" | while read -r file; do
    awk '/^ *1 / {print $2, $4, $5; exit}' "$file" | while read -r vina cnn_pose cnn_aff; do
        if [ -n "$vina" ]; then
            clean_name=$(basename "$file" _results.log)
            echo -e "${clean_name}\t${vina}\t${cnn_pose}\t${cnn_aff}" >> "$summary_file"
        fi
    done
done

# 1. Sort by Vina / Vinardo Empirical Affinity ($2 ascending -> most negative first)
sort -k2,2n "$summary_file" | awk '{print $1 "\t" $2}' > "$output_dir/sorted_vina_affinity.txt"

# 2. Sort by CNN Pose Score ($3 descending -> closest to 1.0 first)
sort -k3,3nr "$summary_file" | awk '{print $1 "\t" $3}' > "$output_dir/sorted_cnn_pose_score.txt"

# 3. Sort by CNN Affinity ($4 descending -> highest pK first)
sort -k4,4nr "$summary_file" | awk '{print $1 "\t" $4}' > "$output_dir/sorted_cnn_affinity.txt"


# Final time
time_end=$(date +%s)
run_time=$((time_end - time_begin))
formatted_time=$(printf '%02d:%02d:%02d' $((run_time/3600)) $((run_time%3600/60)) $((run_time%60)))
echo "======================================================================"
echo "Total run time: ${formatted_time}"
# Give the file system a second to flush the final echo text to disk
sleep 2 

#=====================================================
# Extract the telemetry into a CSV report - we might change this to a diffrent report type later.
#=====================================================

csv_file="${WORKSPACE_DIR}/results/${TARGET_ID}/gnina/grid_search_telemetry.csv"

# Add a header row if the file doesn't exist yet
if [ ! -f "$csv_file" ]; then
    echo "Run_Name,Hardware,Mols,Init_Time_s,Avg_Time_Per_Mol_s,Peak_CPU_Cores,Avg_CPU_Cores,Peak_RAM_MB,Peak_VRAM_MB,Avg_GPU_Util" > "$csv_file"
fi

# Cleanly extract just the numeric values and first words to avoid CSV formatting errors
clean_hw_env=$(echo "$hw_env" | awk -F':' '{print $1}')

# Append this run's data cleanly to the next row
echo "${run_name},${clean_hw_env},${total_mols},${init_time},${avg_mol_time},${max_cpu},${avg_cpu},${max_ram},${max_vram},${avg_gpu}" >> "$csv_file"
#!/bin/bash

time_begin=$(date +%s)
source ~/miniconda3/etc/profile.d/conda.sh
conda activate vina
#Export the active Conda PATH so GNU Parallel subshells inherit it natively
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
    cpu_cores=35
    echo "No SLURM detected, running with $cpu_cores cores."
    export bash_shells=7
fi

#=====================================================
# Paths and vars
#=====================================================

timestamp=$(date +"%Y-%m-%d-%H")

run_name="2IW3_72pH_ex64"

## GNU parallel paths for HPC
batches_dir_path="/mnt/ceph/keel4205/Vina/test_batch_72pH"
export receptor_path="/mnt/ceph/keel4205/Vina/receptors/2IW3_72pH.pdbqt"
export output_dir="/mnt/ceph/keel4205/Vina/results/${run_name}"
export config_vina_path="/mnt/ceph/keel4205/Vina/receptors/2IW3_config.txt"

## GNU parallel paths for local
#batches_dir_path="/home/justin/Jagdish_lab/vina_testing/inputs_batches"
#export receptor_path="/home/justin/Jagdish_lab/vina_testing/test_receptor.pdbqt"
#export output_dir="/home/justin/Jagdish_lab/vina_testing/outputs_batches"
#export config_vina_path="/home/justin/Jagdish_lab/bash_scripts/vina_test_config.txt"



mkdir -p "$output_dir"

#=====================================================
#GNU function
#=====================================================

process_ligand() {

    full_path_to_ligand_var=$1   # $1 will be the full path caught by parallel

    clean_name=$(basename "$full_path_to_ligand_var" .pdbqt) #strips the dir off the name
    output_pdbqt="${output_dir}/${clean_name}_result.pdbqt"  #output file path
    
    temp_log="${output_dir}/${clean_name}_temp.log"
    final_log="${output_dir}/${clean_name}_results.log"

    echo "starting vina for: $clean_name"
    
    # Toggles srun based on the environment
    if [[ -n "$SLURM_CPUS_PER_TASK" ]]; then
        vina_cores="$SLURM_CPUS_PER_TASK"
        EXEC_CMD="srun --exclusive -N1 -n1 --cpus-per-task=$SLURM_CPUS_PER_TASK vina"
    else
        vina_cores=5
        EXEC_CMD="vina"
    fi

    #actual vina stuff
    $EXEC_CMD \
        --receptor "$receptor_path" \
        --ligand "$full_path_to_ligand_var" \
        --config "$config_vina_path" \
        --out "$output_pdbqt" \
        --scoring vinardo \
        --exhaustiveness=64 \
        --cpu "$vina_cores" \
        > "$temp_log" 2>&1
        
        #--my other parameters here

        #exaustiveness is search depth, linear to compute time, default=8

        #num_modes default is 9, it says how many binding modes to output
        #--num_modes=9

        #energy_range default is 3, says how far from the best binding score the other modes can be to be included in the output.
        #--energy_range=3

        #verbosity default is 1, it says how much info to print to the log file. 0 is no info, 1 is only the final results, 2 is more info about the docking process
        #--verbosity=1

        awk 'NR >= 20 && NR <= 28 || NR >= 35 || /Error/ || /error/ || /WARNING/ || /Fatal/' "$temp_log" > "$final_log"
        rm "$temp_log" 
    #use meeko to convert to .sdf
    #mk_export.py --input "$output_pdbqt" --output "${output_dir}/${clean_name}_result.sdf"
    #rm "$output_pdbqt" #emoves temp .pdbqt files

}
export -f process_ligand  #export the function to be used by GNU parallel


echo "sacrificing ligands to the machine gods, this may take a while..."

for i in $(seq 1 2); do
    current_batch="$batches_dir_path/batch_$i"

    if [ -d "$current_batch" ]; then #check if the batch directory exists
        echo "================================="
        echo "|      Processing batch $i      |"
        echo "================================="
        batch_begin=$(date +%s)
        # Point find at the current batch, so GNU parallel only sees one batch at a time
        find "$current_batch" -name "*.pdbqt" | parallel --line-buffer --no-notice --joblog "$output_dir/batch_${i}_run.log" --resume -j "$bash_shells" process_ligand {}
        rm -f "$output_dir"/*_temp.log #cleanup any temp log files that might not have been removed due to errors

        batch_end=$(date +%s)
        batch_time=$((batch_end - batch_begin))
        batch_formatted=$(printf '%02d:%02d:%02d' $((batch_time/3600)) $((batch_time%3600/60)) $((batch_time%60)))
        echo "Batch $i completed in ${batch_formatted}"
    else 
        echo "Warning: $current_batch not found. Skipping."
    fi
done

#copy quick_sort.sh into output dir and run it. 

cp /mnt/ceph/keel4205/Vina/quick_sort.sh "$output_dir/"
cd "$output_dir" || exit
chmod +x quick_sort.sh
./quick_sort.sh

# 2. Calculate time BEFORE moving the logs
time_end=$(date +%s)
run_time=$((time_end - time_begin))
formatted_time=$(printf '%02d:%02d:%02d' $((run_time/3600)) $((run_time%3600/60)) $((run_time%60)))
echo "======================================================================"
echo "Total run time: ${formatted_time}"

# Give the file system a second to flush the final echo text to disk
sleep 2 

# 3. Dynamically move the correct log file
if [[ -n "$SLURM_JOB_ID" ]]; then
    # We are in SLURM. Grab the auto-generated slurm log.
    mv "${SLURM_SUBMIT_DIR}/slurm-${SLURM_JOB_ID}.out" "$output_dir/${run_name}_${timestamp}.out"
else
    # We are local. Grab the nohup log.
    mv "${LAUNCH_DIR}/nohup.out" "$output_dir/${run_name}_${timestamp}.out"
fi
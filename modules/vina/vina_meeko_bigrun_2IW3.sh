#!/bin/bash

#SBATCH -J big_vina_run_2IW3
#SBATCH -p eight
#SBATCH --mem-per-cpu=125M
#SBATCH --ntasks=130
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --exclude=n081

# Array directive: Run jobs numbered 1 through 32
#SBATCH --array=1-32
# Create a separate log file for each batch (%A is the master job ID, %a is the array ID)
#SBATCH -o array_%A_batch_%a.out

#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=keel4205@vandals.uidaho.edu



#-----------------------------------------------------
#Safety check to prevent running on head node
if [[ -z "$SLURM_CPUS_PER_TASK" ]]; then
    echo "====================================================================="
    echo "FATAL ERROR: Slurm variables are missing!"
    echo "====================================================================="
    exit 1
fi

time_begin=$(date +%s)
source ~/miniconda3/etc/profile.d/conda.sh
conda activate vina
#Export the active Conda PATH so GNU Parallel subshells inherit it natively
export PATH
#=====================================================
# Paths and vars
#=====================================================
export bash_shells="$SLURM_NTASKS"
export vina_cores="$SLURM_CPUS_PER_TASK"
export total_cores=$(( bash_shells * vina_cores ))

timestamp=$(date +"%Y-%m-%d-%H")

run_name="array_2IW3_all_mols"
i=$SLURM_ARRAY_TASK_ID

## GNU parallel paths for HPC using $SLURM_TMPDIR
export inputs="/mnt/ceph/keel4205/HTVS/inputs"
export results="/mnt/ceph/keel4205/HTVS/results/${run_name}"
export batch_out_dir="$results/batch_${i}"
export meeko_temp="$batch_out_dir/meeko_temp"
export batch_failed_log="$results/failed_mols_batch_${i}.txt"
export meeko_log="$batch_out_dir/meeko_batch_${i}.log"

export receptor_path="/mnt/ceph/keel4205/Vina/receptors/2IW3_72pH.pdbqt"
export config_vina_path="/mnt/ceph/keel4205/Vina/receptors/2IW3_config.txt"

mkdir -p "$meeko_temp"
mkdir -p "$batch_out_dir"

#=====================================================
#GNU worker function
#=====================================================

process_ligand() {
    ligand_path=$1   # $1 will be the full path caught by parallel
    clean_name=$(basename "$ligand_path" .pdbqt) #strips the dir off the name
    
    output_pdbqt="${batch_out_dir}/${clean_name}.pdbqt"  #output file path
    temp_log="${batch_out_dir}/${clean_name}_temp.log"
    final_log="${batch_out_dir}/${clean_name}_results.log"

    #echo "starting vina for: $clean_name"
    
    #actual vina stuff
    srun --exclusive -N1 -n1 --cpus-per-task=$vina_cores vina \
        --receptor "$receptor_path" \
        --ligand "$ligand_path" \
        --config "$config_vina_path" \
        --out "$output_pdbqt" \
        --scoring vinardo \
        --exhaustiveness=16 \
        --cpu "$vina_cores" \
        > "$temp_log" 2>&1

        # Filter junk but keep errors
        awk 'NR >= 20 && NR <= 28 || NR >= 35 || /Error/ || /error/ || /WARNING/ || /Fatal/' "$temp_log" > "$final_log"
        rm -f "$temp_log" 

        if [ ! -s "$output_pdbqt" ]; then
            echo -e "${clean_name}\t failed Vina docking" >> "$batch_failed_log"
        fi

}
export -f process_ligand  #export the function to be used by GNU parallel




echo "Sacrificing ligands to the machine gods, this may take a while..."

input_chunk="$inputs/input_batch_${i}/chunk_${i}.sdf"

if [ ! -f "$input_chunk" ]; then
    echo "Batch $i not found."
    exit 1
fi

batch_begin=$(date +%s)
echo "================================="
echo "|       Starting batch $i       |"
echo "================================="

echo "Exploding chunk $i with meeko"
mk_prepare_ligand.py -i "$input_chunk" --multimol_outdir "$meeko_temp" > "$meeko_log" 2>&1
grep -Ei "processed:|written:|error|duplicate|traceback|exception|fatal|warning" "$meeko_log"

echo "Docking chunk $i in Vina" 
# Point find at the current batch, so GNU parallel only sees one batch at a time
find "$meeko_temp" -name "*.pdbqt" | parallel --line-buffer --no-notice --joblog "$results/batch_${i}_joblog.log" --resume -j "$bash_shells" process_ligand {}

echo "Sorting and tarballing batch $i results"
tarball_name="batch_${i}_results.tar.gz"

# Find files one at a time to prevent memory overflow, then feed to the loop
find "$batch_out_dir" -name "*_results.log" | while read -r file; do
    
    # Search for the line that starts with "   1" and grab the 2nd column (the score)
    # The 'exit' tells awk to stop searching once it finds the top hit, saving CPU time
    score=$(awk '/^ *1 / {print $2; exit}' "$file")

    if [ -n "$score" ]; then
        # Strip off the path and the extension so the final list is just clean ZINC IDs
        clean_name=$(basename "$file" _results.log)
        echo -e "${clean_name}\t${score}"
    fi

done | sort -k2,2n > "$batch_out_dir/sorted_docking_scores.txt"

echo "Top 10 hits:"
head -n 10 "$batch_out_dir/sorted_docking_scores.txt"

tar -czf "$results/$tarball_name" -C "$results" "batch_${i}"

rm -rf "$meeko_temp"/*
rm -rf "$batch_out_dir"

batch_end=$(date +%s)
batch_time=$((batch_end - batch_begin))
batch_formatted=$(printf '%02d:%02d:%02d' $((batch_time/3600)) $((batch_time%3600/60)) $((batch_time%60)))
echo "Batch $i completed in ${batch_formatted}"

# Final time
time_end=$(date +%s)
run_time=$((time_end - time_begin))
formatted_time=$(printf '%02d:%02d:%02d' $((run_time/3600)) $((run_time%3600/60)) $((run_time%60)))
echo "======================================================================"
echo "Total run time: ${formatted_time}"
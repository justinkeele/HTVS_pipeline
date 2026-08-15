#!/bin/bash

#SBATCH -J big_vina_run_test
#SBATCH -p eight-short
#SBATCH --mem=4G
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00

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

run_name="test2_2IW3_ex16_all_molecules"

## GNU parallel paths for HPC using $SLURM_TMPDIR
export ceph_inputs="/mnt/ceph/keel4205/HTVS/inputs"
export ceph_results="/mnt/ceph/keel4205/HTVS/results/${run_name}"
export receptor_path="/mnt/ceph/keel4205/Vina/receptors/2IW3_72pH.pdbqt"
export config_vina_path="/mnt/ceph/keel4205/Vina/receptors/2IW3_config.txt"

mkdir -p "$ceph_results"

export NODE_TMP="${SLURM_TMPDIR:-/tmp/$USER/$SLURM_JOB_ID}"
mkdir -p "$NODE_TMP"

export tmp_meeko_ready="$NODE_TMP/meeko_ready"
export tmp_vina_results="$NODE_TMP/vina_results"

#=====================================================
# Failsafe save
#=====================================================
# When SLURM terminates early it sends SIGTERM warning
failsafe_save() {
    echo "Warning, triggering failsafe"
    tar -czf "$NODE_TMP/failsafe_save.tar.gz" -C "$tmp_vina_results" .
    cp "$NODE_TMP/failsafe_save.tar.gz" "$ceph_results/"
    echo "failsafe save succeded"
    exit 1
}
trap 'failsafe_save' SIGTERM

#=====================================================
#GNU worker function
#=====================================================

process_ligand() {
    ligand_path=$1   # $1 will be the full path caught by parallel
    clean_name=$(basename "$ligand_path" .pdbqt) #strips the dir off the name
    
    output_pdbqt="${tmp_vina_results}/${clean_name}.pdbqt"  #output file path
    temp_log="${tmp_vina_results}/${clean_name}_temp.log"
    final_log="${tmp_vina_results}/${clean_name}_results.log"

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

for i in $(seq 1 3); do

    ceph_chunk="$ceph_inputs/input_batch_${i}/chunk_${i}.sdf"

    if [ ! -f "$ceph_chunk" ]; then
        echo "Batch $i not found."
        continue
    fi
    batch_begin=$(date +%s)
    echo "================================="
    echo "|       Starting batch $i       |"
    echo "================================="

    #setup directories in tmpdir
    mkdir -p "$tmp_meeko_ready"
    mkdir -p "$tmp_vina_results"

    export batch_failed_log="$ceph_results/failed_mols_batch_${i}.txt"

    echo "Exploding chunk $i with meeko"
    mk_prepare_ligand.py -i "$ceph_chunk" --multimol_outdir "$tmp_meeko_ready"

    echo "Docking chunk $i in Vina"
            # Point find at the current batch, so GNU parallel only sees one batch at a time
    find "$tmp_meeko_ready" -name "*.pdbqt" | head -n 5 | parallel --line-buffer --no-notice --joblog "$tmp_vina_results/batch_${i}_joblog.log" --resume -j "$bash_shells" process_ligand {}


    echo "Sorting and tarballing batch $i results"
    tarball_name="batch_${i}_results.tar.gz"

    mv "$tmp_vina_results" "$NODE_TMP/batch_${i}"

    # Find files one at a time to prevent memory overflow, then feed to the loop
    find . -maxdepth 1 -name "*_results.log" | while read -r file; do
        
        # Search for the line that starts with "   1" and grab the 2nd column (the score)
        # The 'exit' tells awk to stop searching once it finds the top hit, saving CPU time
        score=$(awk '/^ *1 / {print $2; exit}' "$file")

        if [ -n "$score" ]; then
            # Strip off the path and the extension so the final list is just clean ZINC IDs
            clean_name=$(basename "$file" _results.log)
            echo -e "${clean_name}\t${score}"
        fi

    done | sort -k2,2n > sorted_docking_scores.txt

    echo "Top 10 hits:"
    head -n 10 sorted_docking_scores.txt

    tar -czf "$NODE_TMP/$tarball_name" -C "$NODE_TMP" "batch_${i}"

    cp "$NODE_TMP/$tarball_name" "$ceph_results/"
    echo "Batch $i saved to ceph"

    rm -rf "$tmp_meeko_ready"/*
    rm -rf "$NODE_TMP/batch_${i}"/*
    rm -f "$NODE_TMP/$tarball_name"
    mkdir -p "$tmp_vina_results"

    batch_end=$(date +%s)
    batch_time=$((batch_end - batch_begin))
    batch_formatted=$(printf '%02d:%02d:%02d' $((batch_time/3600)) $((batch_time%3600/60)) $((batch_time%60)))
    echo "Batch $i completed in ${batch_formatted}"

done

# Final time
time_end=$(date +%s)
run_time=$((time_end - time_begin))
formatted_time=$(printf '%02d:%02d:%02d' $((run_time/3600)) $((run_time%3600/60)) $((run_time%60)))
echo "======================================================================"
echo "Total run time: ${formatted_time}"
#!/bin/bash

#SBATCH -J vina_run_multiple_exaustivness
#SBATCH -p eight-long
#SBATCH --mem=20G
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=10
#SBATCH --time=120:00:00

#SBATCH --test-only

#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=keel4205@vandals.uidaho.edu

#-----------------------------------------------------
# main variables
cpu_cores=$SLURM_CPUS_PER_TASK

#Safety check to prevent running on head node
if [[ -z "$SLURM_CPUS_PER_TASK" ]]; then
    echo "====================================================================="
    echo "FATAL ERROR: Slurm variables are missing!"
    echo "====================================================================="
    exit 1
fi

# Paranoia Check: Scream and die if we accidentally try to over-allocate on the cluster
if [[ "$cpu_cores" -gt "$SLURM_CPUS_PER_TASK" ]]; then
    echo "====================================================================="
    echo "FATAL ERROR: GNU Parallel cores ($cpu_cores) exceeds Slurm allocation ($SLURM_CPUS_PER_TASK)!"
    echo "Screaming and dying to prevent HPC node throttling."
    echo "====================================================================="
    exit 1
fi

echo "Safety check passed. Running with $cpu_cores cores."

bash vina_worker.sh

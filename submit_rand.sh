#!/bin/bash
# Job array
# MUST MAKE -o directory before submitting.
#SBATCH --job-name=cmssm_random
#SBATCH -p free
#SBATCH --array=0-1000
#SBATCH --cpus-per-task=8
#SBATCH -o /dir/for/standard/output/files/  # STDOUT

omp_threads=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$omp_threads

cd absolute/path/to/MSSM_Scan/src

# Run the object file
absolute/path/to/MSSM_Scan/src/dm_rand.o $SLURM_ARRAY_TASK_ID

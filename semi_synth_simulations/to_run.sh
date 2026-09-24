#!/usr/bin/bash -l
#SBATCH --job-name=arrayJob
#SBATCH --time=02:00:00
#SBATCH --ntasks=7
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=200GB
#SBATCH --output=arrayJob_%A_%a.out
#SBATCH --error=arrayJob_%A_%a.err
#SBATCH --array=1-50

module load rstudio/4.5.2

# Print this sub-job's task ID
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
Rscript --vanilla simulation_run.R $SLURM_ARRAY_TASK_ID

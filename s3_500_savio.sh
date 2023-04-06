#!/bin/bash
#SBATCH --job-name=sys3_n500
#SBATCH --account=co_biostat
#SBATCH --partition=savio3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --time=50:00:30
#SBATCH --mail-user=junming_shi@berkeley.edu
## Command(s) to run:

R CMD BATCH --no-save scripts/scripts_v3/2_simulations_sys3_500.R scripts/scripts_v3/2_simulations_sys3_500.Rout
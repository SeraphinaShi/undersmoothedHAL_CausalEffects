#!/bin/bash
#SBATCH --job-name=sys1_n200
#SBATCH --account=co_biostat
#SBATCH --partition=savio3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --time=50:00:30
#SBATCH --mail-user=junming_shi@berkeley.edu
## Command(s) to run:

R CMD BATCH --no-save scripts/scripts_v5/2_simulations_sys1_200.R scripts/scripts_v5/2_simulations_sys1_200.Rout
#!/bin/bash
#
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -M junming_shi@berkeley.edu
#$ -m base

R CMD BATCH --no-save scripts/scripts_v5_final/2_simulations_sys6.R scripts/scripts_v5_final/2_simulations_sys6.Rout
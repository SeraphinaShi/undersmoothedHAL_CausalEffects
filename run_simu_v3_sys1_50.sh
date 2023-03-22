#!/bin/bash
#
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -M junming_shi@berkeley.edu
#$ -m base

R CMD BATCH --no-save scripts/scripts_v3/2_simulations_sys1_50.R scripts/scripts_v3/2_simulations_sys1_50.Rout
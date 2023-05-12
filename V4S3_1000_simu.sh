#!/bin/bash
#
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -M junming_shi@berkeley.edu
#$ -m base

R CMD BATCH --no-save scripts/scripts_v4/2_simulations_sys3_1000.R scripts/scripts_v4/2_simulations_sys3_1000.Rout
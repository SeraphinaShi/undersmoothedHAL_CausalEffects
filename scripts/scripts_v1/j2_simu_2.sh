#!/bin/bash
#
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -M junming_shi@berkeley.edu
#$ -m base
#

R CMD BATCH --no-save 2_simulations_2.R 2_simulations_2.Rout


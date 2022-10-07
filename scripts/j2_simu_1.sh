#!/bin/bash
#
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -M junming_shi@berkeley.edu
#$ -m base
#

R CMD BATCH --no-save 2_simulation_1.R 2_simulation_1.Rout

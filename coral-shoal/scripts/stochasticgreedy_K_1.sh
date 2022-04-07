#!/bin/bash -l
#
#PBS -N CORAL-GREEDY-K-1-MEDIAN
#PBS -l ncpus=10
#PBS -l cputype=6140
#PBS -l mem=150gb
#PBS -l walltime=2:00:00
#PBS -e scripts/stochasticgreedy_K_1_error.txt
#PBS -o scripts/stochasticgreedy_K_1_output.txt

cd $PBS_O_WORKDIR
module load gdal/3.2.1-foss-2020b
module load r/4.0.3-foss-2020b
R -e "j <- 1; source('./stochastic_greedy_exchange.R');"

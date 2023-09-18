#!/bin/sh
#SBATCH --time=00-24:00:00 # days-hh:mm:ss
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=20
#SBATCH --mem=100GB
#SBATCH --partition nesi_prepost
module load Anaconda2/5.2.0-GCC-7.1.0

#python PP_TEST.py
CMIP6_Data_Processing.py

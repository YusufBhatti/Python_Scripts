#!/bin/sh
#SBATCH --time=00-24:00:00 # days-hh:mm:ss
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=20
#SBATCH --mem=100GB
#SBATCH --partition nesi_prepost
#module load Anaconda2/5.2.0-GCC-7.1.0
export PATH=/nesi/project/niwa02757/ybh10/miniconda3/bin:$PATH
source activate master

#python /Scripts/PP_Chem_DMS.ipynb
#python post_processes.py
#python further_postproc.py
#python MODEL_Stash_PP.py
####python CMIP6_Data_Processing.py
#python CMIP6_HadGEM_PP.py 
#python yes.py
#python Untitled.py
#python DMS_SSA.py
python COPY_of_Interpolating_CMIP6_Data_4_slurm.py

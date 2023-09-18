#!/bin/sh
#SBATCH --time=00-24:00:00 # days-hh:mm:ss
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=20
#SBATCH --mem=400GB
#SBATCH --partition nesi_prepost
#module load Anaconda2/5.2.0-GCC-7.1.0
export PATH=/nesi/project/niwa02757/ybh10/miniconda3/bin:$PATH
source activate master

#python /Scripts/PP_Chem_DMS.ipynb
#python CDF_Calculation.py
#python RNG_Stochastic_Model_Generator.py
#python Ancillary_MODIS_CLIM.py
#python Ancillary_MODIS_CLIM-DAILY.py
python MODIS_SORT_DATA_3D_data.py

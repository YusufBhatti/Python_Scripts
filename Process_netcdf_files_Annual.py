# This script processes daily averaged datasets and compiles it into one file.

import os
os.chdir("/home/ybh10/Scripts/")
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas as pd
import xarray as xr
from netCDF4 import Dataset
import glob
from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.stats import stats
import gc
import math
from scipy.stats import sem

#from my_functions import *
import datetime
dms_ppt=(29/62.13)*1e12 # molecular mass of air (g/mol). molecular mass of DMS (g/mol).
h2so4_mmr=1e14
h2so4_aod=1e15
so2_ppb=(29/64.066)*1e11 # ppb
aitken_h2so4=1e11
accoum_h2so4=1e11
coarse_h2so4=1e13
nucleation_h2so4=1e15

Control_file='/home/ybh10/Objective_2/Post_Processed_Data/MEDUSA_LM86/Years/Processed/'
MEDUSA_LM86_file='/home/ybh10/Objective_2/Post_Processed_Data/MEDUSA_LM86/Years/Processed/'

#########################################################################################
MODIS_LM86_file='/home/ybh10/Objective_2/Post_Processed_Data/MODIS_LM86/Years/Processed/'
CHEM_new_file='/home/ybh10/Objective_2/Post_Processed_Data/CHEM3_MOD_New/Years/Processed/'
MODIS_N00_file='/home/ybh10/Objective_2/Post_Processed_Data/MODIS_N00/Years/Processed/'
MODIS_W14_file='/home/ybh10/Objective_2/Post_Processed_Data/MODIS_W14/Years/Processed/'
MODIS_GM16_file='/home/ybh10/Objective_2/Post_Processed_Data/MODIS_GM16/Years/Processed/'
MODIS_B17_file='/home/ybh10/Objective_3/Postprocess_Data/MODIS_B17/Processed/'
MODIS_B17_CLIM_file='/home/ybh10/Objective_3/Postprocess_Data/MODIS_B17_CLIM/Processed/'

#Perturbed_file='/home/ybh10/Objective_2/Post_Processed_Data/MODIS/Years/Processed/'

#########################################################################################
Lana_file='/home/ybh10/Objective_2/Post_Processed_Data/Lana_LM86/Years/Processed/'
Lana_LM86_file='/home/ybh10/Objective_2/Post_Processed_Data/Lana_LM86/Years/Processed/'

Lana_GM16_file='/home/ybh10/Objective_2/Post_Processed_Data/Lana_GM16/Years/Processed/'
Lana_W14_file='/home/ybh10/Objective_2/Post_Processed_Data/Lana_W14/Years/Processed/'
Lana_B17_file='/home/ybh10/Objective_2/Post_Processed_Data/Lana_B17/Years/Processed/'

Lana_chem_file='/home/ybh10/Objective_2/Post_Processed_Data/Lana_CHEM3_NEW/Years/Processed/'
#########################################################################################
Hulswar_LM86_file='/home/ybh10/Objective_2/Post_Processed_Data/Hulswar_LM86/Years/Processed/'

def standard(data):
    serror=sem(data,axis=0)
    positive=np.nanmean(data,axis=0)
    negative=np.nanmean(data,axis=0)
    pos_err=positive+serror
    neg_err=negative-serror
    return pos_err, neg_err

def find_diffs(run10,run1,axis1):
    t = stats.ttest_ind(run10,run1,axis=axis1); tt=t[1]
    diff = np.where(tt>0.05, 1, 0) #1 indicates it is not statistically significant; 0 indicates it is
    return diff

def sea_ice_mask (sea_ice_data,data_2d,lat,lon,num):
    var = np.where(sea_ice_data>num,100000,data_2d)#mask sea ice
    data=np.empty((np.shape(data_2d)[0],np.shape(data_2d)[1])); data[:]=np.nan
    lat_data=np.empty((np.shape(data_2d)[0],np.shape(data_2d)[1])); data[:]=np.nan
    lon_data=np.empty((np.shape(data_2d)[0],np.shape(data_2d)[1])); data[:]=np.nan
    
    for xla in range(0,lat.shape[0]):
        for ylon in range(lon.shape[1]):
            if var[xla,ylon] > 90000:
                data[xla,ylon]=var[xla,ylon]
            else:
                data[xla,ylon]=np.nan
    return data

def sea_ice_bias(data,pre_num,post_num,sea_ice_data):
    sic_new=np.nanmean(sea_ice_data[-post_num:],axis=(0))
#     pre=np.nanmedian(data[:,:pre_num],axis=(1))
    post=np.nanmedian(data[-post_num:],axis=(0))
#     var_old = np.where(sic==0,np.nan,pre)#mask sea ice
    var_new = np.where(sic_new==0,np.nan,post)#mask sea ice
#     model_diff=var_new-var_old
    return post
def for_sea_ice(data,pre_num,post_num,sea_ice_data):
    sic=np.nanmean(data[:pre_num],axis=(0))
    sic_new=np.nanmean(data[-post_num:],axis=(0))
    var_new = np.where(sic_new==0,np.nan,sic_new)#mask sea ice
    var_old = np.where(sic==0,np.nan,sic)#mask sea ice
    model_diff=var_new-var_old
    
def to_monthly(dats):
    day = dats.time.dt.day
    month = dats.time.dt.month

    # assign new coords
    ds = dats.assign_coords(day=("time", day.data), month=("time", month.data))
    # reshape the array to (..., "month", "year")
    return ds.set_index(time=("day", "month")).unstack("time")  

def to_yearly(dats):
    day = dats.time.dt.day

    year = dats.time.dt.year
    month = dats.time.dt.month

    # assign new coords
    ds = dats.assign_coords(year=("time", year.data), month=("time", month.data))
    # reshape the array to (..., "month", "year")
    return ds.set_index(time=("year", "month")).unstack("time")  
def custom_round(x, base=1):
        return int(base * round(float(x)/base))
    
def to_day_mon_yr(dats):
    day = dats.time.dt.day

    year = dats.time.dt.year
    month = dats.time.dt.month

    # assign new coords
    ds = chl.assign_coords(year=("time", year.data), month=("time", month.data),day=("time", day.data))
    # reshape the array to (..., "month", "year")
    return ds.set_index(time=("year", "month","day")).unstack("time")  
months = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
# file_types=['Lana_LM86','MEDUSA_B17','MODIS_B17','MODIS_B17_CHEM3','MODIS_B17_GEOS_CHEM','MODIS_B17_GFDL'
#            ,'MODIS_B17_MIROC','MODIS_B17_NORESM','MODIS_B17_SOCOL','MODIS_LM86','REV3_B17']
#file_types=['Hulswar_LM86','Lana_B17','Lana_W14','MEDUSA_LM86','MODIS_LM86','MODIS_W14','MODIS_N00']

# file_types=['MEDUSA_B17','MODIS_B17_CHEM3','MODIS_B17_GEOS_CHEM','MODIS_B17_GFDL','MODIS_B17_MIROC','MODIS_B17_NORESM','MODIS_B17_SOCOL','MODIS_LM86',
#            'REV3_B17']
AOD_file=[]
file_types=['MEDUSA_LM86','MODIS_LM86']#,'Lana_LM86','Hulswar_LM86',
        #'MODIS_W14','Lana_W14','MODIS_B17','Lana_B17']
file_types=['MODIS_B17']

for i in range(0,len(file_types)):
    AOD_file.append(globals().get(file_types[i]+'_file'))

#[,
#            ,,'','MODIS_B17_SOCOL','MODIS_LM86',]
for file_type in (AOD_file):
    run_length='Years'
    file_input='{}'.format(file_type)

    year=np.arange(2009,2019,1) 
    direct=file_type

    months = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'] 
    dates=[]
    for y in (year):
        for m in (months):
             dates.append('{}{}'.format(y,m))    
#    dates=dates[12:36]

    Variable_File=['DMS_Mass_Mixing_Ratio','DMS_SURF_EMISSIONS','H2SO4_Mass_Mixing_Ratio',
                   'AITKEN_MODE_SOL_H2SO4_MMR','ACCUMULATION_MODE_SOL_H2SO4_MMR',
                   'COARSE_MODE_SOL_H2SO4_MMR','ACCUMULATION_MODE_SOL_SSA_MMR',
                   'COARSE_MODE_SOL_SSA_MMR','SO2_MASS_MIXING_RATIO','SO2_SURF_EMISSIONS','Sea_Ice_Fraction','NUCLEATION_MODE_SOL_H2SO4_MMR','DRY_PARTICLE_DIAMETER_nucleation_sol', 'DRY_PARTICLE_DIAMETER_aitken_sol', 'DRY_PARTICLE_DIAMETER_accumulation_sol',
              'DRY_PARTICLE_DIAMETER_coarse_sol', 'DRY_PARTICLE_DIAMETER_aitken_insol','number_mixing_ratio_nucleation_sol_mode','number_mixing_ratio_aitken_sol_mode','number_mixing_ratio_accumulation_sol_mode',
             'number_mixing_ratio_coarse_sol_mode','number_mixing_ratio_aitken_insol_mode','Number_nucleation_sol_mode','Number_aitken_sol_mode','Number_accumulation_sol_mode',
           'Number_coarse_sol_mode','Number_aitken_insol_mode',
                   'atmosphere_optical_thickness_due_to_insoluble_aitken_mode_sulphate_aerosol',
                   'atmosphere_optical_thickness_due_to_soluble_accumulation_mode_sulphate_aerosol','atmosphere_optical_thickness_due_to_soluble_aitken_mode_sulphate_aerosol',
                   'atmosphere_optical_thickness_due_to_soluble_coarse_mode_sulphate_aerosol','atmosphere_optical_thickness_due_to_dust_ambient_aerosol','CLEAR_SKY_OUTGOING_SW_FLUX_TOA','ALL_SKY_OUTGOING_SW_FLUX_TOA',
                   'DROPLET_NUMBER_CONCENTRATION','CDNC_Cloud_Top','WEIGHT_FOR_CDNC','ALL_SKY_OUTGOING_LW_FLUX_TOA','CLEAR_SKY_OUTGOING_LW_FLUX_TOA','10M_WIND_U-COMP',
                   'ccn_concentration_25r','H2SO4_nucleation_sol','H2SO4_aitken_sol','H2SO4_accum_sol','H2SO4_coarse_sol','O3_Mass_Mixing_Ratio',
            'CS2+OH=>SO2+COS_2','CS2+OH=>SO2+COS_P','DMS+BrO=>DMSO+Br_2','DMS+BrO=>DMSO+Br_P',
           'DMS+Cl=>0.5SO2+0.5DMSO+0.5HCl+0.5ClO_2','DMS+Cl=>0.5SO2+0.5DMSO+0.5HCl+0.5ClO_P','DMS+NO3=>SO2+HONO2+MeOO+HCHO_2','DMS+NO3=>SO2+HONO2+MeOO+HCHO_P','DMS+O3=>NULL3_2','DMS+O3=>NULL3_P','DMS+O3=>SO2_2',
           'DMS+O3=>SO2_P','DMS+OH=>SO2+DMSO+MeOO_2','DMS+OH=>SO2+DMSO+MeOO+P','DMS+OH=>SO2+MeOO+HCHO_2','DMS+OH=>SO2+MeOO+HCHO_P','DMSO+OH=>MSIA+SO2_2','DMSO+OH=>MSIA+SO2_P','MSIA+O3=>MSA_2','MSIA+O3=>MSA_P',
           'MSIA+O3=>NULL4_2','MSIA+O3=>NULL4_P','MSIA+O3=>NULL5_2','MSIA+O3=>NULL5_P','MSIA+OH=>SO2+MSA_2','MSIA+OH=>SO2+MSA_P','SO2+H2O2=>NULL0_2','SO2+H2O2=>NULL0_P','SO2+HOBr=>NULL6_2','SO2+HOBr=>NULL6_P',
           'SO2+HOBr=>NULL7_2','SO2+HOBr=>NULL7_P','SO2+O3=>NULL1_2','SO2+O3=>NULL1_P','SO2+O3=>NULL2_2','SO2+O3=>NULL2_P','SO2+OH=>HO2+H2SO4_2','SO2+OH=>HO2+H2SO4_P','SO2_natural_emissions','Total_BC_load',
           'Total_H2SO4_load','Total_OM_load','Total_Sea_Salt_load','DMSO_Mass_Mixing_Ratio']
    #Variable_File=['Sea_Ice_Fraction']
    #'

    Variable=['mass_fraction_of_dimethyl_sulfide_in_air','m01s50i214','mass_fraction_of_sulfuric_acid_in_air','mass_fraction_of_sulfuric_acid_in_soluble_aitken_mode_dry_aerosol_in_air', 'mass_fraction_of_sulfuric_acid_in_soluble_accumulation_mode_dry_aerosol_in_air','mass_fraction_of_sulfuric_acid_in_soluble_coarse_mode_dry_aerosol_in_air','mass_fraction_of_seasalt_in_soluble_accumulation_mode_dry_aerosol_in_air','mass_fraction_of_seasalt_in_soluble_coarse_mode_dry_aerosol_in_air',
              'mass_fraction_of_sulfur_dioxide_in_air','m01s50i215','sea_ice_area_fraction','mass_fraction_of_sulfuric_acid_in_soluble_nucleation_mode_dry_aerosol_in_air','m01s38i401','m01s38i402','m01s38i403','m01s38i404','m01s38i405','number_of_particles_per_air_molecule_of_soluble_nucleation_mode_aerosol_in_air','number_of_particles_per_air_molecule_of_soluble_aitken_mode_aerosol_in_air',
                'number_of_particles_per_air_molecule_of_soluble_accumulation_mode_aerosol_in_air','number_of_particles_per_air_molecule_of_soluble_coarse_mode_aerosol_in_air',
                'number_of_particles_per_air_molecule_of_insoluble_aitken_mode_aerosol_in_air','m01s38i504','m01s38i505','m01s38i506','m01s38i507','m01s38i508',
              'atmosphere_optical_thickness_due_to_insoluble_aitken_mode_sulphate_aerosol','atmosphere_optical_thickness_due_to_soluble_accumulation_mode_sulphate_aerosol',
              'atmosphere_optical_thickness_due_to_soluble_aitken_mode_sulphate_aerosol','atmosphere_optical_thickness_due_to_soluble_coarse_mode_sulphate_aerosol','atmosphere_optical_thickness_due_to_dust_ambient_aerosol',
              'toa_outgoing_shortwave_flux_assuming_clear_sky','toa_outgoing_shortwave_flux','product_of_number_concentration_of_stratiform_cloud_liquid_water_particles_and_stratiform_cloud_liquid_water_area_fraction_and_sunlit_binary_mask',
              'm01s01i298','m01s01i299','toa_outgoing_longwave_flux','toa_outgoing_longwave_flux_assuming_clear_sky','x_wind','m01s38i439','m01s38i485','m01s38i486','m01s38i487',
              'm01s38i488','mass_fraction_of_ozone_in_air',
             'm01s50i144','m01s52i144','m01s50i138','m01s52i138','m01s50i137','m01s52i137','m01s50i338','m01s52i338','m01s50i331','m01s52i331','m01s50i339','m01s52i339',
             'm01s50i337','m01s52i337','m01s50i336','m01s52i336','m01s50i136','m01s52i136','m01s50i139','m01s52i139','m01s50i332','m01s52i332','m01s50i333','m01s52i333','m01s50i135','m01s52i135',
             'm01s50i151','m01s52i151','m01s50i334','m01s52i334','m01s50i335','m01s52i335','m01s50i152','m01s52i152','m01s50i153','m01s52i153','m01s50i150','m01s52i150','m01s50i217','m01s38i525',
             'm01s38i520','m01s38i531','m01s38i539','mass_fraction_of_dimethyl_sulfoxide']                   

    mfdata=['DMS_Mass_Mixing_Ratio','DMS_SURF_EMISSIONS','H2SO4_Mass_Mixing_Ratio','AITKEN_MODE_SOL_H2SO4_MMR','ACCUMULATION_MODE_SOL_H2SO4_MMR','COARSE_MODE_SOL_H2SO4_MMR','ACCUMULATION_MODE_SOL_SSA_MMR',
                   'COARSE_MODE_SOL_SSA_MMR','SO2_MASS_MIXING_RATIO',
     'NUCLEATION_MODE_SOL_H2SO4_MMR','atmosphere_optical_thickness_due_to_insoluble_aitken_mode_sulphate_aerosol','atmosphere_optical_thickness_due_to_soluble_coarse_mode_sulphate_aerosol',
     'atmosphere_optical_thickness_due_to_soluble_aitken_mode_sulphate_aerosol','atmosphere_optical_thickness_due_to_soluble_accumulation_mode_sulphate_aerosol','atmosphere_optical_thickness_due_to_dust_ambient_aerosol','CLEAR_SKY_OUTGOING_SW_FLUX_TOA',
    'ALL_SKY_OUTGOING_SW_FLUX_TOA','DROPLET_NUMBER_CONCENTRATION','CDNC_Cloud_Top','WEIGHT_FOR_CDNC','ALL_SKY_OUTGOING_LW_FLUX_TOA','CLEAR_SKY_OUTGOING_LW_FLUX_TOA',
    '10M_WIND_U-COMP','ccn_concentration_25r','H2SO4_nucleation_sol','H2SO4_aitken_sol','H2SO4_accum_sol','H2SO4_coarse_sol','O3_Mass_Mixing_Ratio','CS2+OH=>SO2+COS_2','CS2+OH=>SO2+COS_P','DMS+BrO=>DMSO+Br_2','DMS+BrO=>DMSO+Br_P',
           'DMS+Cl=>0.5SO2+0.5DMSO+0.5HCl+0.5ClO_2','DMS+Cl=>0.5SO2+0.5DMSO+0.5HCl+0.5ClO_P','DMS+NO3=>SO2+HONO2+MeOO+HCHO_2','DMS+NO3=>SO2+HONO2+MeOO+HCHO_P','DMS+O3=>NULL3_2','DMS+O3=>NULL3_P','DMS+O3=>SO2_2',
           'DMS+O3=>SO2_P','DMS+OH=>SO2+DMSO+MeOO_2','DMS+OH=>SO2+DMSO+MeOO+P','DMS+OH=>SO2+MeOO+HCHO_2','DMS+OH=>SO2+MeOO+HCHO_P','DMSO+OH=>MSIA+SO2_2','DMSO+OH=>MSIA+SO2_P','MSIA+O3=>MSA_2','MSIA+O3=>MSA_P',
           'MSIA+O3=>NULL4_2','MSIA+O3=>NULL4_P','MSIA+O3=>NULL5_2','MSIA+O3=>NULL5_P','MSIA+OH=>SO2+MSA_2','MSIA+OH=>SO2+MSA_P','SO2+H2O2=>NULL0_2','SO2+H2O2=>NULL0_P','SO2+HOBr=>NULL6_2','SO2+HOBr=>NULL6_P',
           'SO2+HOBr=>NULL7_2','SO2+HOBr=>NULL7_P','SO2+O3=>NULL1_2','SO2+O3=>NULL1_P','SO2+O3=>NULL2_2','SO2+O3=>NULL2_P','SO2+OH=>HO2+H2SO4_2','SO2+OH=>HO2+H2SO4_P','SO2_natural_emissions','Total_BC_load',
           'Total_H2SO4_load','Total_OM_load','Total_Sea_Salt_load','DMSO_Mass_Mixing_Ratio','DRY_PARTICLE_DIAMETER_nucleation_sol', 'DRY_PARTICLE_DIAMETER_aitken_sol', 'DRY_PARTICLE_DIAMETER_accumulation_sol',
              'DRY_PARTICLE_DIAMETER_coarse_sol', 'DRY_PARTICLE_DIAMETER_aitken_insol']
    if file_type == 'Control':
        mfdata.append('Sea_Ice_Fraction')
    else:
        pass
    #aod_files=sorted(glob.glob('/home/ybh10/Objective_2/Post_Processed_Data/{}/*atmosphere**{}*'.format(model_type,daty)))
    ## DMS_MMR, H2SO4_Mass_Mixing_Ratio, AITKEN_MODE_SOL_H2SO4_MMR, ACCUMULATION_MODE_SOL_H2SO4_MMR, SO2_MASS_MIXING_RATIO, Sea_Ice
    #Variable_File=Variable_File[40:40]
    for i in range(0,len(Variable_File)):
        files=[]
        #if Variable_File[i] == 'DMS_Mass_Mixing_Ratio':
        #     try:
        for d in (dates):
            file=('{}../Raw/{}_{}.nc'.format(direct,Variable_File[i],d))
            files.append(file)


        try:
            try:
                xr.open_mfdataset('{}/../Annual_Processed/{}.nc'.format(direct,Variable_File[i]))
                print('{} = Exists files'.format(Variable_File[i]))

                pass
            except:
                try:
                    data=xr.open_mfdataset(files[:])
                except:
                    try:
                        data=xr.open_mfdataset(files[:],combine='nested',concat_dim='time')
                    except:
                        try:
                            data=xr.open_mfdataset('{}../Raw/{}*'.format(direct,Variable_File[i]))
                        except:
                            data=xr.open_mfdataset('{}../Raw/{}*'.format(direct,Variable_File[i]),combine='nested',concat_dim='time')

                try:
                    data.level_height
                    data = data[Variable[i]].sel(model_level_number=slice(0, 16))


                    print(Variable[i])
                except:
                    pass   
                try: 
                    data.pseudo_level
                    data=data[Variable[i]][2]
                    print(Variable[i])
                except:
                    pass
                print(Variable_File[i])
            #    end
                if Variable_File[i] == 'Sea_Ice_Fraction':
                    pass
                else:
                    data=data.load()
                data.to_netcdf('{}/../Annual_Processd_Obj_2/{}.nc'.format(direct,Variable_File[i]))
        except:
            print('{} = no files'.format(Variable_File[i]))
#     pass

    #  except:

        gc.collect()

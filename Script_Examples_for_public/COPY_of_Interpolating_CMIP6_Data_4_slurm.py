## CMIP 6 
import os
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from netCDF4 import Dataset
import glob
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import cartopy.feature as cfeature
AOD=[]
import xarray as xr
#os.chdir("/home/ybh10/Scripts/")
#from my_functions import *
import gc
First=0
first1=0
Last=0
def pick_a_season(season):
    month_data=seas(season) # information for the script below.
    if season == 'DJF' !=-1: # The first = first month in slice ## IF NOT DJF
        First=month_data[0]
        first1=month_data[1] ### ONLY USE IF DJF - It should be 0
        Last=month_data[2] # The last month in the slice. Should be the last month.
        return First,first1,Last
    else:
        First=month_data[0]
        Last=month_data[1]
        return First,Last

def seas (season): # For when defining First to Last in months for the seasons
    if season == 'DJF' !=-1:
        First=11 # The first = first month in slice ## IF NOT DJF
        first1=0 ## ONLY USE IF DJF - It should be 0
        Last=1 # The last month in the slice. Should be the last month.
        return First, first1, Last
    if season == 'MAM' !=-1:
        First=2
        Last=5
    if season == 'JJA' !=-1:
        First=5
        Last=8
    if season == 'SON' !=-1:
        First=8
        Last=11
    return First, Last

def CMIP6(globe,dms_PPT,global_array):
    """
        This function will give me all seasons when specified the first and last!
    """
    season=np.arange(0,len(globe),12)
    if (First)==11 !=-1: # Only activated when DJF seasonality is being used.
        for seas in (season):
            p=((globe[First+seas])*dms_PPT)
            j=((globe[first1+seas])*dms_PPT)
            i=((globe[Last+seas])*dms_PPT)
            globally=np.nanmean((p,j,i),axis=0)
            global_array.append(globally)
    else: # Activated when MAM, JA, SON is being created.
        for seas in (season):
                globally=(np.nanmean(globe[First+seas:Last+seas],axis=0)*dms_PPT)
                global_array.append(globally)
    return global_array

def process_cmip6(no_of_variants,files,variable_name,season_type): 
    """
        Processes RAW CMIP6 taking in the number of 
        files within one variant_name (one), concentating 
        it within one dataset (full model run in one array)
        and sorts each year into seasons and compiles each model
        ensemble member as a seperate dataset. This function also
        incorporates ONE SEASON (WHICH_SEASON), for a mean
        Resulting in two datasets saved: 
                        each year with each ONE season mean
                        each year with monthly mean.
    """
    hi=0
    no_season_globe=[]
    season_array=[]
    globe=[]
    Seasoned_ensemble=[]
    #ODMS_Global_array_year=[]
#     for var in (variable_name):
    for f in (files):
        if f.find (variable_name) !=-1:  # Oceanic DMS
            print(f)
            vn_rawdata=xr.open_dataset(f)
            vn_annual_data=vn_rawdata[variable_name]
            if vn_annual_data.ndim==4:
                vn_annual_data=vn_annual_data[:,0]
            else:
                pass
            #no_season_globe.append(vn_annual_data)
            seasonal_data=CMIP6(vn_annual_data.data,1,globe)
            season_array.append(globe)
            globe=[]
            hi=hi+1
            #checking number of iterations within ensemble member. 
            print(hi)
            print('------------')
            if hi == no_of_variants !=-1: #if number of iterations is = number in the file ==:
                if hi == 1 !=-1:     
                    total=((season_array[0]))
                if hi == 2 !=-1:
                    total=np.concatenate((season_array[0],season_array[1]))
                if hi == 3 !=-1:
                    total=np.concatenate((season_array[0],season_array[1],season_array[2]))
                if hi == 4 !=-1:
                    total=np.concatenate((season_array[0],season_array[1],season_array[2],season_array[3]))
                if hi == 5 !=-1:
                    total=np.concatenate((season_array[0],season_array[1],season_array[2],season_array[3],season_array[4]))
                if hi == 8 !=-1:
                    total=np.concatenate((season_array[0],season_array[1],season_array[2],season_array[3],season_array[4],season_array[5],season_array[6],season_array[7]))
                #total2=np.concatenate((no_season_globe[0],no_season_globe[1]))
                Seasoned_ensemble.append(total)
                #ODMS_Global_array_year.append(total2)
                test=[]
                #no_season_globe=[]
                hi=0 
    return Seasoned_ensemble

# end_sentence=''
# file_output='/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Ozone_Region/{}/Ensamble/'.format(seasons)  
# global_file_output='/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/{}/Ensamble/'.format(seasons)     

seasons=['SON','DJF'] ## Choose the season you want, then the code and function will collate the correct 

    
#from my_functions import *


data_type='Sea_Ice'
#name=['uas','ua','uas']

#name=['uas','ua','uas']
Route='Indirect'

if Route == 'Direct':  #### DIRECT ####
    if data_type == 'Oceanic_DMS' or 'DMS_flux '!=-1:
    ##### OCEANIC_DMS & DMS_FLUX #####
        models=['UKESM1','NorESM2-LM','MIROC-ES2L']
        itertions=[2,8,1]
    if data_type == 'DMS'!=-1:
    ##### Atmospheric DMS #####
            models=['UKESM1','NorESM2-LM','MIROC-ES2L']
            itertions=[3,8,1]
    if data_type == 'uas'!=-1:
        print('uas_direct')
        models=['UKESM1','NorESM2-LM','MIROC-ES2L']
        itertions=[2,8,1]
        name=['uas','ua','uas']
    if data_type == 'Ozone_Column' !=-1:
##### Column Ozone #####
        models=['UKESM1','CNRM-CM6-1','GISS-E2-1-G','MRI-ESM2','GFDL','IPSL-CM6A-LR-INCA']
        itertions=[2,1,4,1,2,1] 

if Route == 'Indirect':  #### INDIRECT ####
    if data_type == 'AOD' or data_type == 'SSA_MMR' !=-1:
        print('aod/ss')
        models=['UKESM1','BCC','CESM2-WACCM','GFDL-CM4','GISS-E2-1-G','HadGEM3','MIROC-ES2L','MPI','MRI-ESM2-0','NorESM2-LM','GFDL-ESM4','CESM2_']
        itertions=[2,1,1,2,3,2,1,5,1,8,2,1]
    if data_type == 'DMS' !=-1:    
        models=['UKESM1','BCC','CESM2-WACCM','GFDL-CM4','GISS-E2-1-G','HadGEM3','MIROC-ES2L','MPI','MRI-ESM2-0','NorESM2-LM','GFDL-ESM4','CESM2_']
        itertions=[2,1,1,2,3,2,1,6,8,8,2,4] 
    if data_type == 'uas' !=-1:
        models=['UKESM1','BCC','CESM2-WACCM','GFDL-CM4','GISS-E2-1-G','HadGEM3','MIROC-ES2L','MPI','MRI-ESM2-0','NorESM2-LM','GFDL-ESM4','CESM2_']
        name=['uas','uas','ua','uas','uas','uas','uas','uas','uas','ua','uas','ua']
        itertions=[2,1,1,2,3,2,1,5,1,8,2,4] 
    if data_type == 'Sea_Ice' !=-1:
    ##### Column Ozone #####
        models=['UKESM1','BCC','CESM2-WACCM','GFDL-CM4','GISS-E2-1-G','HadGEM3','MIROC-ES2L','MPI','MRI-ESM2-0','NorESM2-LM','GFDL-ESM4','CESM2_']
        itertions=[2,1,1,2,3,2,1,5,1,8,2,1]
        name=['siconca','siconc','siconc','siconc','siconca','siconc','siconc','siconca','siconc','siconc','siconc','siconc']


i=0
SON_data=[]
DJF_data=[]
for season_mode in (seasons):
    print(season_mode)
    for model,number,count in zip(models,itertions,range(0,len(models))):
        if data_type == 'uas' or data_type == 'Sea_Ice' !=-1:
            i = count
        file = sorted(glob.glob('/nesi/project/niwa02757/ybh10/CMIP6/{}/{}/Interpolated_Grid/*{}*.nc'.format(Route,data_type,model)))
        if season_mode == 'SON'!=-1:

            First,Last=pick_a_season(season_mode)
            print(name[i])

            if data_type == 'Ozone_Column' and count == 5 :
                data_array=np.array(process_cmip6(number,file,name[i],season_mode))/100000
                SON_data.append(data_array)
            else:
                SON_data.append(process_cmip6(number,file,name[i],season_mode))                                    
        else:
            First,Last,first1=pick_a_season(season_mode)
            if data_type == 'Ozone_Column' and count == 5 :
            #    data_array=np.array(process_cmip6(number,file,name[i],season_mode))/100000
                DJF_data.append(data_array)
            else:
                DJF_data.append(process_cmip6(number,file,name[i],season_mode)) 


for i in range(0,len(SON_data)):
    print(np.shape(SON_data[i]))
    print(models[i])
    print('-------')
    
file_output='/home/ybh10/CMIP6/{}/{}/Processed/'.format(Route,data_type)  
# #######################################################################################
# # LATITUDE MEAN FOR THIS DATASET
## np.save('{}{}_Mean_Latitude.npy'.format(file_output,data_type),mean_lat)
## np.save('{}{}_Mean_Longitude.npy'.format(file_output,data_type),mean_lon)
# ## Latitude & LONGITUDE MEAN
## np.save('{}{}_each_models_Latitude.npy'.format(file_output,data_type),latitudes)
## np.save('{}{}_each_models_Longitude.npy'.format(file_output,data_type),longitudes)
## np.save('{}{}_Mean_Coords.npy'.format(file_output,data_type),mean_lon)
# #######################################################################################
#Saving the data from above within ONE file and 4 seperate files.
np.save('{}SON/RAW/{}_SON_DATA.npy'.format(file_output,data_type),SON_data)
np.save('{}DJF/RAW/{}_DJF_DATA.npy'.format(file_output,data_type),DJF_data)
for model,i in zip(models,range(0,len(models))):
    np.save('{}SON/RAW/{}_{}_SON_DATA.npy'.format(file_output,data_type,model),SON_data[i])
    np.save('{}DJF/RAW/{}_{}_DJF_DATA.npy'.format(file_output,data_type,model),DJF_data[i])
gc.collect()
#np.save('{}{}CHEG_Global{}.npy'.format(global_file_output,Ensamble_output,end_sentence),Ensamble_Global_CHEG)

## INTERPOLATION OF MODEL DATA ## 

SON_data=np.load('/home/ybh10/CMIP6/{}/{}/Processed/SON/RAW/{}_SON_DATA.npy'.format(Route,data_type,data_type),allow_pickle=True)
DJF_data=np.load('/home/ybh10/CMIP6/{}/{}/Processed/DJF/RAW/{}_DJF_DATA.npy'.format(Route,data_type,data_type),allow_pickle=True)
lat=xr.open_dataset('/nesi/project/niwa02757/ybh10/CMIP6/DMS_flux/fgdms_Omon_UKESM1-0-LL_historical_r12i1p1f2_gn_185001-194912.nc').lat
lont=xr.open_dataset('/nesi/project/niwa02757/ybh10/CMIP6/DMS_flux/fgdms_Omon_UKESM1-0-LL_historical_r12i1p1f2_gn_185001-194912.nc').lon

# UKESM=np.array(SON_data[0])
# CNRM=np.array(SON_data[1])
# NOR=np.array(SON_data[2])
# MIROC=np.array(SON_data[3])
m1=0;m2=0;m3=0;m4=0;m5=0;m6=0
m7=0;m8=0;m9=0;m10=0;m11=0;m12=0

models=[m1,m2,m3]#,m4,m5,m6,m7,m8,m9,m10,m11,m12]
models_mean=[m1,m2,m3]#,m4,m5,m6,m7,m8,m9,m10,m11,m12]
season_data=[SON_data,DJF_data]
years=np.arange(1850,2015,1)
for season_type,seaso in zip(season_data,seasons):
    for model,i in zip(models,range(0,len(models))):
        models[i]=np.array(season_type[i])
        models_mean[i]=np.nanmean(season_type[i],axis=0)
        print('mid of file {}'.format(i))
        print('mid file {}'.format(seaso))
        models[i] = models[i][:,-75::]
        models_mean[i] = models_mean[i][-75::]

    print('end of file {}'.format(i))
    print('end of file {}'.format(seaso))

    np.save('{}{}/{}_{}_model_mean.npy'.format(file_output,seaso,data_type,seaso),models_mean)
    
# UKESM=models[0];NOR=models[1];MIROC=models[2]
# UKESM_mean=models[0];NOR_mean=models_mean[1];MIROC_mean=models_mean[2]

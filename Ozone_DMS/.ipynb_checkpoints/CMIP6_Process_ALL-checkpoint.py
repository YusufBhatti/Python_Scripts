# CMIP 6 
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
os.chdir("/home/ybh10/Scripts/")
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
            vn_annual_data=vn_rawdata[variable_name][:]
            #no_season_globe.append(vn_annual_data)
            seasonal_data=CMIP6(vn_annual_data,1,globe)
            season_array.append(globe)
            globe=[]
            hi=hi+1
            #checking number of iterations within ensemble member. 
            print(hi)
            if hi == no_of_variants !=-1: #if number of iterations is = number in the file ==:
                if hi == 2 !=-1:
                    total=np.concatenate((season_array[0],season_array[1]))
                if hi == 1 !=-1:     
                    total=((season_array[0]))
                if hi == 8 !=-1:
                    total=np.concatenate((season_array[0],season_array[1],season_array[2],season_array[3],season_array[4],season_array[5],season_array[6],season_array[7]))
                #total2=np.concatenate((no_season_globe[0],no_season_globe[1]))
                Seasoned_ensemble.append(total)
                #ODMS_Global_array_year.append(total2)
                test=[]
                #no_season_globe=[]
                hi=0 
    return Seasoned_ensemble

end_sentence=''
file_output='/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Ozone_Region/{}/Ensamble/'.format(seasons)  
global_file_output='/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/{}/Ensamble/'.format(seasons)     

seasons=['SON','DJF'] ## Choose the season you want, then the code and function will collate the correct 

    
#from my_functions import *
os.chdir("/nesi/project/niwa02757/ybh10/CMIP6/UKESM1/Historic/")
Historic='Historic'
models=['UKESM1','CNRM-ESM2-1','NorESM2-LM','MIROC-ES2L']
itertions=[2,1,8,1]
#DMS_flux_file = sorted(glob.glob('/nesi/project/niwa02757/ybh10/CMIP6/DMS_flux/*{}*.nc'.format(model)))


files=[ODMS_file,DMS_flux_file]
name=['fgdms']
data_type='DMS_flux'
SON_data=[]
DJF_data=[]
for season_mode in (seasons):
    print(season_mode)
    for model,number in zip(models,itertions):
        for i in range(0,len(name)):
            file = sorted(glob.glob('/nesi/project/niwa02757/ybh10/CMIP6/{}/*{}*.nc'.format(data_type,model)))
            if season_mode == 'SON':
                First,Last=pick_a_season(season_mode)
                SON_data.append(process_cmip6(number,file,name[i],season_mode)) 
            else:
                First,Last,first1=pick_a_season(season_mode)
                DJF_data.append(process_cmip6(number,file,name[i],season_mode)) 

end_sentence=''
file_output='/home/ybh10/CMIP6/{}/Processed/'.format(data_type)  
#global_file_output='/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/{}/'.format(seasons)     
#Ensamble_output='Ensamble/'

# #######################################################################################
# CHEG
#np.save('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/Ensamble/CHEG_Global.npy',CHEG_Global_array)      
np.save('{}{}_SON_DATA.npy'.format(file_output,data_type),SON_data)
np.save('{}{}_DJF_DATA.npy'.format(file_output,data_type),DJF_data)

#np.save('{}{}CHEG_Global{}.npy'.format(global_file_output,Ensamble_output,end_sentence),Ensamble_Global_CHEG)
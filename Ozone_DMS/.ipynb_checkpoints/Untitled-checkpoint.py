# CMIP 6 
import os
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
from netCDF4 import Dataset
import glob
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.stats import stats
import gc
from scipy.stats import sem

AOD=[]
import xarray as xr
os.chdir("/home/ybh10/Scripts/")

#from my_functions import *
os.chdir("/nesi/project/niwa02757/ybh10/CMIP6/UKESM1/Historic/")

vn1=xr.open_dataset('/nesi/project/niwa02757/ybh10/CMIP6/UKESM1/Historic/Raw_Data/DMS/dms_AERmon_UKESM1-0-LL_historical_r10i1p1f2_gn_185001-189912.nc')
lats=vn1.lat
lons=vn1.lon
time=vn1.time

DMS_ppt=((29/62.13)*1e12)
ODMS_nm=1e6
ODMS_nm_nor=1e9
#ODMS_trill=1e12
#ODMS_ppt=(((18/62.13)*1e12))
O3_ppm=((29/48)*1e6)
years=np.arange(1850,2015,1)
years_ozone=np.arange(1940,2015,1)
direct_season_SON='SON'
Seasons=['SON','DJF']
Direct_variables_of_interest=['Oceanic_DMS','DMS','Ozone_Column','uas']
Indirect_variables_of_interest=['AOD','DMS','Sea_Ice','SSA_MMR','uas','DMS_SSA']
models=['UKESM1','BCC','CESM2-WACCM','GFDL-CM4','GISS-E2-1-G','HadGEM3','MIROC-ES2L','MPI','MRI-ESM2-0','NorESM2-LM','GFDL-ESM4','CESM2_']
direct_models=['UKESM1','NorESM2-LM','MIROC-ES2L']

#Indirect_variables_of_interest=['UKESM1','BCC','CESM2-WACCM','GFDL-CM4','GISS-E2-1-G','HadGEM3','MIROC-ES2L','MPI','MRI-ESM2-0','NorESM2-LM','GFDL-ESM4','CESM2_']
ppt=[ODMS_nm,DMS_ppt,1,1,1]
ppt_indirect=[1,DMS_ppt,1,1,1,1]

ODMS_direct_SON=[];DMS_direct_SON=[]; Ozone_Column_direct_SON=[];uas_direct_SON=[]
direct_data_SON=[ODMS_direct_SON,DMS_direct_SON,Ozone_Column_direct_SON,uas_direct_SON]

ODMS_direct_DJF=[];DMS_direct_DJF=[]; Ozone_Column_direct_DJF=[]; uas_direct_DJF=[]
direct_data_DJF=[ODMS_direct_DJF,DMS_direct_DJF,Ozone_Column_direct_DJF,uas_direct_DJF]

SON_direct_data_ensemble=[]
DJF_direct_data_ensemble=[]
# SON_Direct_Data=[]
# DJF_Direct_Data=[]
# Seasons_Direct_Data=[SON_Direct_Data,DJF_Direct_Data]
######################################################################################
################################## DIRECT VARIABLES ##################################
######################################################################################
i=0
m=0
for Season in (Seasons):
     for data_type,data_input_son,data_input_djf in zip(Direct_variables_of_interest,direct_data_SON,direct_data_DJF):
        directory_path='/home/ybh10/CMIP6/Direct/{}/Processed/{}'.format(data_type,Season)
        file_load=np.load('{}/{}_{}_model_mean.npy'.format(directory_path,data_type,Season))
        #multi_file_load=np.load('{}/RAW/{}_{}_DATA.npy'.format(directory_path,data_type,direct_season_SON))
        data=np.squeeze(file_load)
        data_mean=np.nanmean(data,axis=0)
        print('loading {} onto the datarray.... np.shape({})'.format(data_type,np.squeeze(np.shape(data_mean))))
   #     data_input.append(data_mean)
        if Season == 'SON':
            SON_direct_data_ensemble.append(file_load*ppt[i])
            i=i+1

        if Season == 'DJF':
            DJF_direct_data_ensemble.append(file_load*ppt[m])
            m=m+1
SON_direct_data_ensemble[0][1]=(SON_direct_data_ensemble[0][1]*1e3)
DJF_direct_data_ensemble[0][1]=(DJF_direct_data_ensemble[0][1]*1e3)

SON_direct_data_ensemble=np.array(SON_direct_data_ensemble)
DJF_direct_data_ensemble=np.array(DJF_direct_data_ensemble)
    #seasons_data.append(data_input)
#direct_data=np.squeeze(direct_data)

######################################################################################
################################## INDIRECT VARIABLES ##################################
######################################################################################

SON_indirect_data_ensemble=[]
DJF_indirect_data_ensemble=[]
# SON_Direct_Data=[]
# DJF_Direct_Data=[]
# Seasons_Direct_Data=[SON_Direct_Data,DJF_Direct_Data]
i=0
m=0
for Season in (Seasons):
    for data_type in (Indirect_variables_of_interest):
        directory_path='/home/ybh10/CMIP6/Indirect/{}/Processed/{}'.format(data_type,Season)
        file_load=np.load('{}/{}_{}_model_mean.npy'.format(directory_path,data_type,Season))
        #multi_file_load=np.load('{}/RAW/{}_{}_DATA.npy'.format(directory_path,data_type,direct_season_SON))
        data=np.squeeze(file_load)
        data_mean=np.nanmean(data,axis=0)
        print('loading {} onto the datarray.... np.shape({})'.format(data_type,np.squeeze(np.shape(data_mean))))
        if Season == 'SON':
            SON_indirect_data_ensemble.append(file_load*ppt_indirect[i])
            i=i+1
            #data_input_son.append(data_mean)
        if Season == 'DJF':
            DJF_indirect_data_ensemble.append(file_load*ppt_indirect[m])
            print('m={}'.format(m))
            m=m+1
#data_input_djf.append(data_mean)

    #seasons_data.append(data_input)
#direct_data=np.squeeze(direct_data)

SON_indirect_data_ensemble=np.array(SON_indirect_data_ensemble)
DJF_indirect_data_ensemble=np.array(DJF_indirect_data_ensemble)

################################################################
######################## DMS / SSA Code ########################
################################################################

# dms_ssa_djf=np.empty((12, 75, 144,192)); dms_ssa_djf[:]=np.nan
# dms_ssa_son=np.empty((12, 75, 144,192)); dms_ssa_son[:]=np.nan

# for n in range(0,12):
#     for i in range(0,75):
#         for lo in range(0,144):
#             for la in range(0,192):
#                 dms_ssa_djf[n,i,lo,la]=DJF_indirect_data_ensemble[1][n,i,lo,la]/DJF_indirect_data_ensemble[3][n,i,lo,la]
#                 dms_ssa_son[n,i,lo,la]=SON_indirect_data_ensemble[1][n,i,lo,la]/SON_indirect_data_ensemble[3][n,i,lo,la]

#     print(np.nanmean(dms_ssa_djf))
# daty=[dms_ssa_son,dms_ssa_djf]
# # SON_indirect_data_ensemble.append(dms_ssa_son)
# # DJF_indirect_data_ensemble.append(dms_ssa_djf)
# for season,dats in zip(Seasons,daty):
#     np.save('/home/ybh10/CMIP6/Indirect/DMS_SSA/Processed/{}/DMS_SSA_{}_model_mean.npy'.format(season,season),dats)

################################################################
################################################################
################################################################



gc.collect()

def standard(data):
    serror=sem(data)
    pos_err=data+serror
    neg_err=data-serror
    if np.ndim(data)==2 !=-1:
        positive=np.nanmean(pos_err,axis=0)
        negative=np.nanmean(neg_err,axis=0)
    return positive, negative

def CMIP6(regional,globe,dms_PPT,array,global_array):  # This function will give me all seasons when specified the first and last!
    season=np.arange(0,len(regional),12)
    if (First)==11 !=-1: # Only activated when DJF seasonality is being used.
        for seas in (season):
            p=((globe[First+seas]))
            j=((globe[first1+seas]))
            i=((globe[Last+seas]))
            globally=np.nanmean((p,j,i))
            global_array.append(globally)
    else: # Activated when MAM, JA, SON is being created.
        for seas in (season):
                regionally=(np.nanmean(regional[First+seas:Last+seas],axis=0)*dms_PPT)
                globally=(np.nanmean(globe[First+seas:Last+seas],axis=0)*dms_PPT)
                array.append(regionally)
                global_array.append(globally)
    return array,global_array

# ##########################################################################################################
# ######################################## RELATIVE PERCENTAGE #############################################
# ##########################################################################################################
# ######################################## PERCENTILES DATA  #############################################

import copy

directs=['Direct']
seasos=['SON']
for dirs, seas in zip(directs,seasos):
    Direct_or_Indirect= dirs ######## CHANGABLE ATTRIBUTION ########
    Season_type=seas           ######## CHANGABLE ATTRIBUTION ########
    dataset=[]


    if Direct_or_Indirect == 'Direct':
        variables_for_calculation=copy.deepcopy(Direct_variables_of_interest) ######## CHANGABLE ATTRIBUTION ########
        if Season_type =='SON':
            data_in= copy.deepcopy(SON_direct_data_ensemble)                     ######## CHANGABLE ATTRIBUTION ########
        else:
            data_in= copy.deepcopy(DJF_direct_data_ensemble)                     ######## CHANGABLE ATTRIBUTION ########
    else:
        variables_for_calculation=copy.deepcopy(Indirect_variables_of_interest) ######## CHANGABLE ATTRIBUTION ########
        if Season_type =='SON':
            data_for_calculations= copy.deepcopy(SON_indirect_data_ensemble)                     ######## CHANGABLE ATTRIBUTION ########
        else:
            data_for_calculations= copy.deepcopy(DJF_indirect_data_ensemble)                     ######## CHANGABLE ATTRIBUTION ########
        nums=[0,1,3,5,7,8,9,10]
        data_in=np.zeros((6,len(nums),75,144,192)) ;  data_in[:]=np.nan
        model_var=[]
        for t in range(0,len(data_in)):
            for i,l in zip(nums,range(0,len(nums))):
                data_in[t,l]=data_for_calculations[t][i]
                model_var.append(models[i])
                if t == 0:
                    print('models = {}'.format(models[i]))

    PERCENTILE=84

    for datas,variable_label in zip(data_in,variables_for_calculation):
        print(variable_label) # iterating through each variables in the loop.
        print(np.shape(datas))
        old_variable=np.nanpercentile(datas[:,:21], PERCENTILE,axis=(0,1))
        #old_variable=np.nanmedian(percentile_data[:21,:,:],axis=(0)) # Taking the mean for the 1940 - 60 climatology .
        vari=np.empty((65,144,192)); vari[:]=np.nan # creating empty array for data through from 1950 --> 2014(5)
        new_data=np.nanpercentile(datas, PERCENTILE,axis=(0))

       # new_data=np.nanmedian(datas,axis=0)
        #data=np.nanmean(datas[num],axis=(0))
      #  for t in range(0,datas[num].shape[0]):
        for i in range(0,vari.shape[0]):
            for la in range(0,new_data.shape[1]):
                for lo in range(0,new_data.shape[2]):
                      #  print('hi')
                    if variable_label == 'uas' !=-1:
                        vari[i,la,lo]=(new_data[i+10,la,lo]-old_variable[la,lo])
                    else:
                        if variable_label == 'dms_ssa' !=-1:
                            vari[i,la,lo]=(new_data[i+10,la,lo]-old_variable[la,lo])
                        else:
                            vari[i,la,lo]=((new_data[i+10,la,lo]/old_variable[la,lo])-1)*100    
        dataset.append(vari)

    np.save('/home/ybh10/CMIP6/{}/Rel_Diff/Relative_Difference_{}_Dataset_{}_MEDIAN_{}%'.format(Direct_or_Indirect,Direct_or_Indirect,Season_type,PERCENTILE),dataset)    
    # ##########################################################################################################
    # ######################################## Absolute Difference #############################################
    # ##########################################################################################################
    dataset_abs=[]

    for datas,variable_label in zip(data_in,variables_for_calculation):
        print(variable_label) # iterating through each variables in the loop.
        print(np.shape(datas))
        old_variable=np.nanpercentile(datas[:,:21,:,:], PERCENTILE,axis=(0,1))
        #old_variable=np.nanmedian(percentile_data[:21,:,:],axis=(0)) # Taking the mean for the 1940 - 60 climatology .
        vari=np.empty((65,144,192)); vari[:]=np.nan # creating empty array for data through from 1950 --> 2014(5)
        new_data=np.nanpercentile(datas, PERCENTILE,axis=(0))

        #data=np.nanmean(datas[num],axis=(0))
      #  for t in range(0,datas[num].shape[0]):
        for i in range(0,vari.shape[0]):
            for la in range(0,new_data.shape[1]):
                for lo in range(0,new_data.shape[2]):
                      #  print('hi')
                    vari[i,la,lo]=((new_data[i+10,la,lo]-old_variable[la,lo]))    
        dataset_abs.append(vari)

    np.save('/home/ybh10/CMIP6/{}/Abs_Diff/Absolute_Difference_{}_Dataset_{}_MEDIAN_{}%'.format(Direct_or_Indirect,Direct_or_Indirect,Season_type,PERCENTILE),dataset_abs)    

#     def find_diffs(run10,run1,axis1):
#         t = stats.ttest_ind(run10,run1,axis=axis1); tt=t[1]
#         diff = np.where(tt>0.05, 1, 0) #1 indicates it is not statistically significant; 0 indicates it is
#         return diff

#     spatial_dots=np.empty((len(data_in), 144,192)); spatial_dots[:]=np.nan
#     for i in range(0,len(data_in)):
#         print("Loading and mapping {} for Statistical Significance".format(variables_for_calculation[i]))
#         percentile_data=np.percentile(data_in[i], PERCENTILE,axis=0)
#         Post_Ozone=percentile_data[-21:]
#         Pre_Ozone=percentile_data[:21]
#         data_diff=find_diffs(Post_Ozone,Pre_Ozone,0)
#         spatial_dots[i,:,:]=(data_diff)
#     #data_dots=np.array(data_dots)

#     Zonal_dots=np.empty((len(data_in), 75, 144)); Zonal_dots[:]=np.nan

#     for var in range(0,len(data_in)):
#         print(variables_for_calculation[var])
#         for i in range(0,75):
#             Post_Ozone=np.nanpercentile(data_in[var], PERCENTILE,axis=0)
#     #        Post_Ozone=(np.nanmedian(data_in[var][:,:],axis=0))
#             percentile_data=np.nanpercentile(data_in[var][:,:21], PERCENTILE,axis=0)
#             Pre_Ozone=np.nanmedian(percentile_data,axis=(0))
#             data_diff=find_diffs(Post_Ozone[i],Pre_Ozone,1)
#             Zonal_dots[var,i]=(data_diff)
#     Zonal_dots=np.array(Zonal_dots)

#     np.save('/home/ybh10/CMIP6/{}/Dots/{}_Zonal_Dots_MEDIAN_{}%.npy'.format(Direct_or_Indirect,Season_type,PERCENTILE),Zonal_dots)
#     np.save('/home/ybh10/CMIP6/{}/Dots/{}_Spatial_Dots_MEDIAN_{}%.npy'.format(Direct_or_Indirect,Season_type,PERCENTILE),spatial_dots)

#     gc.collect()


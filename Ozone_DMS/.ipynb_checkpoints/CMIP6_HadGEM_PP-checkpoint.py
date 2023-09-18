#!/usr/bin/env python
# coding: utf-8

# In[21]:


# CMIP 6 
import os
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xrvn1
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
from my_functions import *
os.chdir("/nesi/project/niwa02757/ybh10/CMIP6/UKESM1/Historic/")
Oz='Ozone_loss'
SO2='SO2'
DMS='DMS'
Historic='Historic'
ssp245='SSP245'
model='HadGEM'
uk='UKESM1'
first_part='/nesi/project/niwa02757/ybh10/CMIP6'
final_part='Historic/Raw_Data'
dms_file = sorted(glob.glob('{}/{}/{}/DMS/*.nc'.format(first_part,model,final_part)))
SO2_file = sorted(glob.glob('{}/{}/{}/SO2/*.nc'.format(first_part,model,final_part)))
ODMS_file = sorted(glob.glob('{}/{}/{}/Ocean_DMS/*.nc'.format(first_part,uk,final_part)))
O3_file = sorted(glob.glob('{}/{}/{}/Ozone_loss//*.nc'.format(first_part,uk,final_part)))
AOD_file = sorted(glob.glob('{}/{}/{}/AOD/*.nc'.format(first_part,model,final_part)))
WIND_file = sorted(glob.glob('{}/{}/{}/Wind/*.nc'.format(first_part,model,final_part)))
RFSW_CS_file = sorted(glob.glob('{}/{}/{}/Radiative_Forcing/Clear_Sky/TOA/SW/*.nc'.format(first_part,model,final_part)))
Surface_RFSW_CS_file = sorted(glob.glob('{}/{}/{}/Radiative_Forcing/Clear_Sky/Surface/SW/*.nc'.format(first_part,model,final_part)))

# ssp_dms_file = sorted(glob.glob('/nesi/project/niwa02757/ybh10/CMIP6/UKESM1/{}/DMS/*.nc'.format(ssp245)))
# ssp_SO2_file = sorted(glob.glob('/nesi/project/niwa02757/ybh10/CMIP6/UKESM1/{}/SO2/*.nc'.format(ssp245)))
# ssp_ODMS_file = sorted(glob.glob('/nesi/project/niwa02757/ybh10/CMIP6/UKESM1/{}/Ocean_DMS/*.nc'.format(ssp245)))
# ssp_O3_file = sorted(glob.glob('/nesi/project/niwa02757/ybh10/CMIP6/UKESM1/{}/Ozone_loss//*.nc'.format(ssp245)))
# ssp_AOD_file = sorted(glob.glob('/nesi/project/niwa02757/ybh10/CMIP6/UKESM1/{}/AOD/*.nc'.format(ssp245)))
# ssp_WIND_file = sorted(glob.glob('/nesi/project/niwa02757/ybh10/CMIP6/UKESM1/{}/Wind/*.nc'.format(ssp245)))

MODIS_Clim = np.load('/home/ybh10/DMS_Emissions/Chemistry_Scheme/CHEM/MODIS2007_2014.npy')
aod_variables=['atmosphere_optical_thickness_due_to_dust_ambient_aerosol','atmosphere_optical_thickness_due_to_insoluble_aitken_mode_sulphate_aerosol',
'atmosphere_optical_thickness_due_to_soluble_accumulation_mode_sulphate_aerosol','atmosphere_optical_thickness_due_to_soluble_aitken_mode_sulphate_aerosol',
'atmosphere_optical_thickness_due_to_soluble_coarse_mode_sulphate_aerosol']

os.chdir("/nesi/nobackup/niwa02757/ybh10/DMS")
inpath_lana=sorted(glob.glob('Lana/2003/*.nc'))
inpath_medusa=sorted(glob.glob('MEDUSA/2003/*.nc'))
months = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']

vn2=xr.open_dataset(ODMS_file[3])
olat=vn2.latitude
olon=vn2.longitude

vn1=xr.open_dataset('/nesi/project/niwa02757/ybh10/CMIP6/UKESM1/Historic/Raw_Data/DMS/dms_AERmon_UKESM1-0-LL_historical_r10i1p1f2_gn_185001-189912.nc')
lat=vn1.lat
lon=vn1.lon
time=vn1.time
# Dates=['January', 'February','March','April','May','June','July','August','September','October','November','December']
lat_bnds = [-90, -60]
inpath_modis=sorted(glob.glob('/nesi/nobackup/niwa02757/ybh10/Observational_Data/MODIS/AOD/*.nc'))
latty_o=xr.open_dataset(inpath_modis[0])
lat_o=latty_o['lat']
lon_o=latty_o['lon']
# vn1=xr.open_dataset('/nesi/project/niwa02757/ybh10/CMIP6/UKESM1/Historic/Raw_Data/Ocean_DMS/dmsos_Omon_UKESM1-0-LL_historical_r10i1p1f2_gn_185001-194912.nc')
# Olat=vn1.latitude
# Dates=['January', 'February','March','April','May','June','July','August','September','October','November','December']

# In[2]:


names=['DMS','SO2','Ozone_Column','AOD','Wind']
def file_load(files):
    ozone_file=[]
    if np.size(files)==4:
        for i in range(0,4):
            a=np.load(files[i])
            ozone_file.append(a)
        return ozone_file
Historic='Historic'
ssp245='SSP245'
season=['DJF','MAM','JJA','SON']
#ozone_area_file=[]
ozone_filey=[]
global_filey=[]

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

def CMIP6(regional,globe,dms_PPT,array,global_array):  # This function will give me all seasons when specified the first and last!
    season=np.arange(0,len(regional),12)
    if (First)==11 !=-1: # Only activated when DJF seasonality is being used.
        for seas in (season):
            l=((regional[First+seas])*dms_PPT)
            h=((regional[first1+seas])*dms_PPT)
            m=((regional[Last+seas])*dms_PPT)
            regionally=np.nanmean((l,h,m),axis=0)
            array.append(regionally)
            p=((globe[First+seas])*dms_PPT)
            j=((globe[first1+seas])*dms_PPT)
            i=((globe[Last+seas])*dms_PPT)
            globally=np.nanmean((p,j,i),axis=0)
            global_array.append(globally)
    else: # Activated when MAM, JA, SON is being created.
        for seas in (season):
                regionally=(np.nanmean(regional[First+seas:Last+seas],axis=0)*dms_PPT)
                globally=(np.nanmean(globe[First+seas:Last+seas],axis=0)*dms_PPT)
                array.append(regionally)
                global_array.append(globally)
    return array,global_array

# for var in (names):
#     ozone_area_file=[]
#     global_area_file=[]
#     for s in (season):
#         ozone_file_name=sorted(glob.glob("/nesi/project/niwa02757/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Ozone_Region/{}/{}_60_90S.npy".format(s,var)))[0]
#         global_file_name=sorted(glob.glob("/nesi/project/niwa02757/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/{}/{}_Global.npy".format(s,var)))[0]
#         ozone_area_file.append(ozone_file_name)
#         global_area_file.append(global_file_name)
#     ozone_file=file_load(ozone_area_file)
#     global_file=file_load(global_area_file)
#     ozone_filey.append(ozone_file)
#     global_filey.append(global_file)


  #  return
## To COMPILE all ensamble data for a VARIABLE into ONE npy ARR
def ensamble_mon(files,file,variable,regional_variable): # files = f (all files in that variable); file = x (1 file)

    if len(files) ==26 !=-1: # This function assumes TWO ensamble files in ONE run
      #  print(file)
        if file.find ('1850')!=-1:
            monthy_1850.append(variable)
            monthy_reg_1850.append(regional_variable)
        if file.find ('2014')!=-1:
            monthy_2014.append(variable) # e.g. variable=vn_o3.dmsos[:]
            monthy_reg_2014.append(regional_variable) # e.g. regional_variable=vn_o3.dmsos[:,0:94,:]
        Total=[monthy_1850,monthy_2014] # neatly puts appended data into one GLOBAL VARIABLE.
        Total_reg=[monthy_reg_1850,monthy_reg_2014] # for REGIONAL VARIABLE.
    if len(files) >27 !=-1:# This function assumes FOUR ensamble files in ONE run
   #     print('hi')
        if file.find ('1850')!=-1: 
            monthy_1850.append(variable)
            monthy_reg_1850.append(regional_variable)
        if file.find ('1900')!=-1:
            monthy_1900.append(variable)
            monthy_reg_1900.append(regional_variable)
        if file.find ('1950')!=-1:
            monthy_1950.append(variable)
            monthy_reg_1950.append(regional_variable)
        if file.find ('2014')!=-1:
            monthy_2014.append(variable)
            monthy_reg_2014.append(regional_variable)
        Total=[monthy_1850,monthy_1900,monthy_1950,monthy_2014]
        Total_reg=[monthy_reg_1850,monthy_reg_1900,monthy_reg_1950,monthy_reg_2014]
    return Total,Total_reg
#    return global_var,regional_var

def ensamble_mean(Total):
    monthy_1850=[]; monthy_1900=[]; monthy_1950=[]; monthy_2014=[]
    monthy_reg_1850=[]; monthy_reg_1900=[]; monthy_reg_1950=[]; monthy_reg_2014=[]
    if len(Total[0]) ==2 !=-1:
        global_1850=np.nanmean(Total[0][0],axis=0); reg_1850=np.nanmean(Total[1][0],axis=0)
        global_2014=np.nanmean(Total[0][1],axis=0); reg_2014=np.nanmean(Total[1][1],axis=0)
        global_var=np.concatenate((global_1850,global_2014)); regional_var=np.concatenate((reg_1850,reg_2014))
    if len(Total[0]) ==4 !=-1:
        global_1850=np.nanmean(Total[0][0],axis=0); reg_1850=np.nanmean(Total[1][0],axis=0)
        global_1900=np.nanmean(Total[0][1],axis=0); reg_1900=np.nanmean(Total[1][1],axis=0)
        global_1950=np.nanmean(Total[0][2],axis=0); reg_1950=np.nanmean(Total[1][2],axis=0)
        global_2014=np.nanmean(Total[0][3],axis=0); reg_2014=np.nanmean(Total[1][3],axis=0)
        global_var=np.concatenate((global_1850,global_1900,global_1950,global_2014)); 
        print('global_var')
        regional_var=np.concatenate((reg_1850,reg_1900,reg_1950,reg_2014))
 
    return global_var,regional_var
# # Oceanic_DMS_Ozone_files=file_load(sorted(glob.glob("/nesi/project/niwa02757/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Ozone_Region/*/*{}*".format(Oceanic_name))))
# # Oceanic_DMS_Global_files=file_load(sorted(glob.glob("/nesi/project/niwa02757/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/*/*{}*".format(Oceanic_name))))

names=['DMS','SO2','Ozone_Column','AOD','Wind']

def ensamble_total(f,runs,var,var_reg,regional_array,global_array):
    if len(f) ==26 !=-1:
        print('--------------------')
        model_run=np.concatenate((var[0],var[1]))
        model_run_reg=np.concatenate((var_reg[0],var_reg[1]))
        ensamble.append(model_run)
        runs.append(model_run_reg)
        hi=0
        if x.find (f[len(f)-1]) !=-1:
            print('-----  LAST ------------------')
            annual_mean=np.nanmean(ensamble,axis=0)
            regional_mean=np.nanmean(ensamble_reg,axis=0)
            dmsx=CMIP6(regional_mean,annual_mean,1,regional_array,global_array)
            return dmsx
    if len(f) >27 !=-1:
        print('--------------------')
        model_run=np.concatenate((var[0],var[1],var[2],var[3]))
        model_run_reg=np.concatenate((var_reg[0],var_reg[1],var_reg[2],var_reg[3]))

        ensamble.append(model_run)
        runs.append(model_run_reg)

        if x.find (f[len(f)-1]) !=-1:
            print('-----  LAST ------------------')
#             var=[]
#             var_reg=[]
            annual_mean=np.nanmean(ensamble,axis=0)
            regional_mean=np.nanmean(ensamble_reg,axis=0)
            dmsx=CMIP6(regional_mean,annual_mean,1,regional_array,global_array)
            
            return dmsx
#count=0
def ensambley(count,test,no_season_SH,no_season_globe,test1,xr_variable,xr_variable_regional,Ensamble_global_array,Ensamble_regional_array,Ensamble_no_season_global_array,Ensamble_no_season_regional_array):
    globe=[]
    reg=[]
    no_season_globe.append(xr_variable)
    no_season_SH.append(xr_variable_regional)
    dmsx=CMIP6(xr_variable_regional,xr_variable,1,reg,globe)
    test.append(globe)
    test1.append(reg)
 #   count=count+1
    print(count)
    if count == 4 !=-1:
        total=np.concatenate((test[0],test[1],test[2],test[3]))
        total1=np.concatenate((test1[0],test1[1],test1[2],test1[3]))
        Ensamble_global_array.append(total)
        Ensamble_regional_array.append(total1)
        total2=np.concatenate((no_season_globe[0],no_season_globe[1],no_season_globe[2],no_season_globe[3]))
        total3=np.concatenate((no_season_SH[0],no_season_SH[1],no_season_SH[2],no_season_SH[3]))
        Ensamble_no_season_global_array.append(total2)
        Ensamble_no_season_regional_array.append(total3)
        test1=[]
        test=[]
        no_season_SH=[]
        no_season_globe=[]
#        count=0
        if x.find (f[len(f)-1]) !=-1:
            print('-----  LAST ------------------')
            return dmsx
    


# In[ ]:





# In[31]:

# PROCESSING THE RAW DATA IF NEW DATA IS AVAILABLE #######
########## 
## This code 

seasons='SON' ## Choose the season you want, then the code and function will collate the correct 
month_data=seas(seasons) # information for the script below.
if seasons == 'DJF' !=-1: # The first = first month in slice ## IF NOT DJF
    First=month_data[0]
    first1=month_data[1] ### ONLY USE IF DJF - It should be 0
    Last=month_data[2] # The last month in the slice. Should be the last month.
else:
    First=month_data[0]
    Last=month_data[1]

end_sentence=''
file_output='/home/ybh10/CMIP6/{}/Historic/Numpy_Array/Ozone_Region/{}/Ensamble/'.format(model,seasons)  
global_file_output='/home/ybh10/CMIP6/{}/Historic/Numpy_Array/Global/{}/Ensamble/'.format(model,seasons)     

monthy_1850=[]; monthy_1900=[]; monthy_1950=[]; monthy_2014=[]
monthy_reg_1850=[]; monthy_reg_1900=[]; monthy_reg_1950=[]; monthy_reg_2014=[]
name=['dms_','so2_','dmsos_Omon','o3_','toz','od550','dmsos_60','uas','rsutcs','rlutcs','rsdscs']
    #DMS, SO2, Sea-DMS, Ozone, Column Ozone, AOD, N/A,         Wind,  TOA SWCS, TOA LWCS, Surface SW,

files=[dms_file,SO2_file,ODMS_file,O3_file,AOD_file,WIND_file,RFSW_CS_file,Surface_RFSW_CS_file]
var=[]
var_reg=[]
ensamble=[]
ensamble_reg=[]

DMS_SH=[]; DMS_Global_array=[]; 
SO2_SH=[]; SO2_Global_array=[];Ozone_SH=[];  Ozone_Global_array=[]
ODMS_SH=[]; ODMS_Global_array=[]
OZ_SH=[];  OZ_Global_array=[]
AOD_SH=[];  AOD_Global_array=[]
Wind_SH=[];  Wind_Global_array=[]
RFSW_SH=[];  RFSW_Global_array=[]
Surface_RFSW_SH=[]; Surface_RFSW_Global_array=[]

reg=[]; globe=[]; test=[]; test1=[]
Ensamble_Global_DMS=[]; Ensamble_Regional_DMS=[]
Ensamble_Global_SO2=[]; Ensamble_Regional_SO2=[]
Ensamble_Global_Ozone=[]; Ensamble_Regional_Ozone=[]
Ensamble_Global_AOD=[]; Ensamble_Regional_AOD=[]
Ensamble_Global_Wind=[]; Ensamble_Regional_Wind=[]
Ensamble_Global_Ocean_DMS=[]; Ensamble_Regional_Ocean_DMS=[]
Ensamble_Global_RFSW=[]; Ensamble_Regional_RFSW=[]
Ensamble_Global_Surface_RFSW=[]; Ensamble_Regional_Surface_RFSW=[]

test=[]; test1=[]; no_season_SH=[]; no_season_globe=[]
hi=0
dms_runs=[]; so2_runs=[];ozone_runs=[];odms_runs=[];aod_runs=[];wind_runs=[];
for f in (files):
    for x in (f):
        lat_inds = np.where((lat > lat_bnds[0]) & (lat < lat_bnds[1]))
        lat_inds = np.squeeze(lat_inds)
#         if x.find (name[0]) !=-1:  # DMS
#             print(x)
#             vn_dms=xr.open_dataset(x)
#             dms_globe=vn_dms.dms[:,0,:,:]
#             dms=vn_dms.dms[:,0,lat_inds[0]:lat_inds[(len(lat_inds)-1)]+1,:]
#             no_season_SH.append(dms) ## These two lines just get ALL the months of the year
#             no_season_globe.append(dms_globe) # They only need to be saved ONCE
#             dmsx=CMIP6(dms,dms_globe,1,reg,globe)
#             test.append(globe)
#             test1.append(reg)
#             globe=[]
#             reg=[]
#             hi=hi+1
#             print(hi)
#             if hi == 4 !=-1:
#                 total=np.concatenate((test[0],test[1],test[2],test[3]))
#                 total1=np.concatenate((test1[0],test1[1],test1[2],test1[3]))
#                 total2=np.concatenate((no_season_globe[0],no_season_globe[1],no_season_globe[2],no_season_globe[3]))
#                 total3=np.concatenate((no_season_SH[0],no_season_SH[1],no_season_SH[2],no_season_SH[3]))
#                 Ensamble_Global_DMS.append(total)
#                 Ensamble_Regional_DMS.append(total1)
#                 DMS_Global_array.append(total2)
#                 DMS_SH.append(total3)
#                 test1=[]
#                 test=[]
#                 no_season_SH=[]
#                 no_season_globe=[]
#                 hi=0
#         if x.find (name[1]) !=-1: # SO2
#             print(x)
#             vn_so2=xr.open_dataset(x)
#             so2_globe=vn_so2.so2[:,0,:,:]
#             so2=vn_so2.so2[:,0,lat_inds[0]:lat_inds[(len(lat_inds)-1)]+1,:]
#             no_season_SH.append(so2)
#             no_season_globe.append(so2_globe)
#             dmsx=CMIP6(so2,so2_globe,1,reg,globe)
#             test.append(globe)
#             test1.append(reg)
#             globe=[]
#             reg=[]
#             hi=hi+1
#             print(hi)
#             if hi == 4 !=-1:
#                 total=np.concatenate((test[0],test[1],test[2],test[3]))
#                 total1=np.concatenate((test1[0],test1[1],test1[2],test1[3]))
#                 Ensamble_Global_SO2.append(total)
#                 Ensamble_Regional_SO2.append(total1)
#                 total2=np.concatenate((no_season_globe[0],no_season_globe[1],no_season_globe[2],no_season_globe[3]))
#                 total3=np.concatenate((no_season_SH[0],no_season_SH[1],no_season_SH[2],no_season_SH[3]))
#                 SO2_Global_array.append(total2)
#                 SO2_SH.append(total3)
#                 test1=[]
#                 test=[]
#                 no_season_SH=[]
#                 no_season_globe=[]
#                 hi=0
#         if x.find (name[4]) !=-1: # Ozone Column
#             print(x)
#             vn_o3=xr.open_dataset(x)
#             o3_globe=vn_o3.toz[:,:,:]
#             o3=vn_o3.toz[:,lat_inds[0]:lat_inds[(len(lat_inds)-1)]+1,:]
#             months=np.arange(0,len(o3_globe),30)
#             aod_mon=[]
#             aod_gl_mon=[]
#             for m in (months):
#                 aod_globe=np.nanmean(o3_globe[m:m+30,:,:],axis=0)
#                 aod=np.nanmean(o3[m:m+30,:,:],axis=0)
#                 aod_gl_mon.append(aod_globe)
#                 aod_mon.append(aod)
#             no_season_SH.append(aod_mon)
#             no_season_globe.append(aod_gl_mon)
#             dmsx=CMIP6(aod_mon,aod_gl_mon,1e5,reg,globe)
#             test.append(globe)
#             test1.append(reg)
#             globe=[]
#             reg=[]
#             hi=hi+1
#             print(hi)
#             if hi == 2 !=-1:
#                 total=np.concatenate((test[0],test[1]))
#                 total1=np.concatenate((test1[0],test1[1]))
#                 total2=np.concatenate((no_season_globe[0],no_season_globe[1]))
#                 total3=np.concatenate((no_season_SH[0],no_season_SH[1]))
#                 Ensamble_Global_Ozone.append(total)
#                 Ensamble_Regional_Ozone.append(total1)
#                 Ozone_Global_array.append(total2)
#                 Ozone_SH.append(total3)
#                 test1=[]
#                 test=[]
#                 no_season_SH=[]
#                 no_season_globe=[]
#                 hi=0             
        if x.find (name[5]) !=-1:   #AOD
            print(x)
  #          print(f)
            vn_o3=xr.open_dataset(x)
            o3_globe=vn_o3.od550aer[:,:,:]
#             o3=vn_o3.od550aer[:,lat_inds[0]:lat_inds[(len(lat_inds)-1)]+1,:]
            months=np.arange(0,len(o3_globe),30)
            aod_mon=[]
            aod_gl_mon=[]
            for m in (months):
                aod_globe=np.nanmean(o3_globe[m:m+30,:,:],axis=0)
           #     aod=np.nanmean(o3[m:m+30,:,:],axis=0)
                aod_gl_mon.append(aod_globe)
     #           aod_mon.append(aod)
        #    no_season_SH.append(aod_mon)
            no_season_globe.append(aod_gl_mon)
            dmsx=CMIP6(aod_gl_mon,aod_gl_mon,1,reg,globe)
            test.append(globe)
      #      test1.append(reg)
            globe=[]
            reg=[]
            hi=hi+1
            print(hi)
            if hi == 2 !=-1:
                total=np.concatenate((test[0],test[1]))
          #      total1=np.concatenate((test1[0],test1[1]))
                total2=np.concatenate((no_season_globe[0],no_season_globe[1]))
         #       total3=np.concatenate((no_season_SH[0],no_season_SH[1]))
                Ensamble_Global_AOD.append(total)
         #       Ensamble_Regional_AOD.append(total1)
                AOD_Global_array.append(total2)
         #       AOD_SH.append(total3)
                test1=[]
                test=[]
                no_season_SH=[]
                no_season_globe=[]
                hi=0  
#         if x.find (name[7]) !=-1:   #Wind
#             print(x)
#             vn_o3=xr.open_dataset(x)
#             o3_globe=vn_o3.uas[:]
#             o3=vn_o3.uas[:,lat_inds[0]:lat_inds[(len(lat_inds)-1)]+1,:]
#             months=np.arange(0,len(o3_globe),30)
#             aod_mon=[]
#             aod_gl_mon=[]
#             for m in (months):
#                 aod_globe=np.nanmean(o3_globe[m:m+30,:,:],axis=0)
#                 aod=np.nanmean(o3[m:m+30,:,:],axis=0)
#                 aod_gl_mon.append(aod_globe)
#                 aod_mon.append(aod)
#             no_season_SH.append(aod_mon)
#             no_season_globe.append(aod_gl_mon)
#             dmsx=CMIP6(aod_mon,aod_gl_mon,1,reg,globe)
#             test.append(globe)
#             test1.append(reg)
#             globe=[]
#             reg=[]
#             hi=hi+1
#             print(hi)
#             if hi == 2 !=-1:
#                 total=np.concatenate((test[0],test[1]))
#                 total1=np.concatenate((test1[0],test1[1]))
#                 total2=np.concatenate((no_season_globe[0],no_season_globe[1]))
#                 total3=np.concatenate((no_season_SH[0],no_season_SH[1]))
#                 Ensamble_Global_Wind.append(total)
#                 Ensamble_Regional_Wind.append(total1)
#                 Wind_Global_array.append(total2)
#                 Wind_SH.append(total3)
#                 test1=[]
#                 test=[]
#                 no_season_SH=[]
#                 no_season_globe=[]
#                 hi=0 
#         if x.find (name[2]) !=-1: # FOR OCEANIC DMS --> NEMO MODEL 
#             vn_o3=xr.open_dataset(x)
#             abc=vn_o3.dmsos[:]
#             abc_reg=vn_o3.dmsos[:,0:94,:] # Requires changing IF changing latitude constraint
#             no_season_SH.append(abc_reg)
#             no_season_globe.append(abc)
#             dmsx=CMIP6(abc_reg,abc,1,reg,globe)
#             test.append(globe)
#             test1.append(reg)
#             globe=[]
#             reg=[]
#             hi=hi+1
#             print(hi)
#             if hi == 2 !=-1:
#                 total=np.concatenate((test[0],test[1]))
#                 total1=np.concatenate((test1[0],test1[1]))
#                 total2=np.concatenate((no_season_globe[0],no_season_globe[1]))
#                 total3=np.concatenate((no_season_SH[0],no_season_SH[1]))
#                 Ensamble_Global_Ocean_DMS.append(total)
#                 Ensamble_Regional_Ocean_DMS.append(total1)
#                 ODMS_Global_array.append(total2)
#                 ODMS_SH.append(total3)
#                 test1=[]
#                 test=[]
#                 no_season_SH=[]
#                 no_season_globe=[]
#                 hi=0       
#         if x.find (name[8]) !=-1:   #RFSW
#             print(x)
#             vn_o3=xr.open_dataset(x)
#             o3_globe=vn_o3.rsutcs[:]
#             o3=vn_o3.rsutcs[:,lat_inds[0]:lat_inds[(len(lat_inds)-1)]+1,:]
#             months=np.arange(0,len(o3_globe),30)
#             aod_mon=[]
#             aod_gl_mon=[]
#             for m in (months):
#                 aod_globe=np.nanmean(o3_globe[m:m+30,:,:],axis=0)
#                 aod=np.nanmean(o3[m:m+30,:,:],axis=0)
#                 aod_gl_mon.append(aod_globe)
#                 aod_mon.append(aod)
#             no_season_SH.append(aod_mon)
#             no_season_globe.append(aod_gl_mon)
#             dmsx=CMIP6(aod_mon,aod_gl_mon,1,reg,globe)
#             test.append(globe)
#             test1.append(reg)
#             globe=[]
#             reg=[]
#             hi=hi+1
#             print(hi)
#             if hi == 2 !=-1:
#                 total=np.concatenate((test[0],test[1]))
#                 total1=np.concatenate((test1[0],test1[1]))
#                 total2=np.concatenate((no_season_globe[0],no_season_globe[1]))
#                 total3=np.concatenate((no_season_SH[0],no_season_SH[1]))
#                 Ensamble_Global_RFSW.append(total)
#                 Ensamble_Regional_RFSW.append(total1)
#                 RFSW_Global_array.append(total2)
#         #        Wind_SH.append(total3)
#                 test1=[]
#                 test=[]
#                 no_season_SH=[]
#                 no_season_globe=[]
#                 hi=0 
#         if x.find (name[10]) !=-1:   #Surface SW
#             print(x)
#             vn_o3=xr.open_dataset(x)
#             o3_globe=vn_o3.rsdscs[:]
#             o3=vn_o3.rsdscs[:,lat_inds[0]:lat_inds[(len(lat_inds)-1)]+1,:]
#             months=np.arange(0,len(o3_globe),30)
#             aod_mon=[]
#             aod_gl_mon=[]
#             for m in (months):
#                 aod_globe=np.nanmean(o3_globe[m:m+30,:,:],axis=0)
#                 aod=np.nanmean(o3[m:m+30,:,:],axis=0)
#                 aod_gl_mon.append(aod_globe)
#                 aod_mon.append(aod)
#             no_season_SH.append(aod_mon)
#             no_season_globe.append(aod_gl_mon)
#             dmsx=CMIP6(aod_mon,aod_gl_mon,1,reg,globe)
#             test.append(globe)
#             test1.append(reg)
#             globe=[]
#             reg=[]
#             hi=hi+1
#             print(hi)
#             if hi == 2 !=-1:
#                 total=np.concatenate((test[0],test[1]))
#                 total1=np.concatenate((test1[0],test1[1]))
#                 total2=np.concatenate((no_season_globe[0],no_season_globe[1]))
#                 total3=np.concatenate((no_season_SH[0],no_season_SH[1]))
#                 Ensamble_Global_Surface_RFSW.append(total)
#                 Ensamble_Regional_Surface_RFSW.append(total1)
#                 Surface_RFSW_Global_array.append(total2)
#                 Surface_RFSW_SH.append(total3)
#                 test1=[]
#                 test=[]
#                 no_season_SH=[]
#                 no_season_globe=[]
#                 hi=0 


sw=sorted(glob.glob('{}/{}/{}/Radiative_Forcing/TOA/SW/*'.format(first_part,model,final_part)))
lw=sorted(glob.glob('{}/{}/{}/Radiative_Forcing/TOA/LW/*'.format(first_part,model,final_part)))
hi=0
rf_1850=[]
globe=[]
reg=[]
test=[]; test1=[]; no_season_SH=[]; no_season_globe=[]
Ensamble_Global_RF=[]; Ensamble_Regional_RF=[];RF_Global_array=[]
for s,l in zip(sw,lw):
    hi=hi+1
    print(hi)
    if s.find ('rsutcs_C') !=-1:
        vn_rf=xr.open_dataset(s)
        swrf=vn_rf.rsutcs
        print(s)
    if l.find ('rlutcs_C') !=-1:
        vn_rf=xr.open_dataset(l)
        lwrf=vn_rf.rlutcs
        rf=swrf+lwrf
        print(l)

        months=np.arange(0,len(rf),30)
        aod_mon=[]
        aod_gl_mon=[]
        for m in (months):
            aod_globe=np.nanmean(rf[m:m+30,:,:],axis=0)
       #     aod=np.nanmean(o3[m:m+30,:,:],axis=0)
            aod_gl_mon.append(aod_globe)
        #    aod_mon.append(aod)
        #    rf_1850.append(swrf+lwrf)
        print('--------------------------------')
        no_season_globe.append(aod_gl_mon)
        dmsx=CMIP6(aod_gl_mon,aod_gl_mon,1,reg,globe)
        test.append(globe)
   #     test1.append(reg)
        globe=[]
        reg=[]
    if hi == 2 !=-1:
        total=np.concatenate((test[0],test[1]))
     #   total1=np.concatenate((test1[0],test1[1]))
        total2=np.concatenate((no_season_globe[0],no_season_globe[1]))
    #    total3=np.concatenate((no_season_SH[0],no_season_SH[1]))
        Ensamble_Global_RF.append(total)
   #    Ensamble_Regional_RF.append(total1)
        RF_Global_array.append(total2)
       # Wind_SH.append(total3)
        test1=[]
        test=[]
        no_season_SH=[]
        no_season_globe=[]
        hi=0 

#         hi=0 

# In[168]:



# Saves ALL arrays into the correct directory. Run this cell after each time you change the season.
# The first two lines for each variable will only need to be run ONCE as it saves ALL the months in the 165 year timescale.
#seasons='SON' # Already defined above.
end_sentence=''
file_output='{}/{}/Historic/Numpy_Array/Ozone_Region/{}/'.format(first_part,model,seasons)  
global_file_output='{}/{}/Historic/Numpy_Array/Global/{}/'.format(first_part,model,seasons)     
Ensamble_output='Ensamble/'

#np.save('/home/ybh10/CMIP6/{}/Historic/Numpy_Array/Ozone_Region/Ensamble/DMS_60_90S.npy'.format(model),DMS_SH)
#np.save('/home/ybh10/CMIP6/{}/Historic/Numpy_Array/Global/Ensamble/DMS_Global.npy'.format(model),DMS_Global_array)
#np.save('{}DMS_60_90S{}.npy'.format(file_output,end_sentence),np.nanmean(Ensamble_Regional_DMS,axis=0))
#np.save('{}{}DMS_60_90S{}.npy'.format(file_output,Ensamble_output,end_sentence),Ensamble_Regional_DMS)
#np.save('{}DMS_Global{}.npy'.format(global_file_output,end_sentence),np.nanmean(Ensamble_Global_DMS,axis=0))              
#np.save('{}{}DMS_Global{}.npy'.format(global_file_output,Ensamble_output,end_sentence),Ensamble_Global_DMS)              
# ######################################################################################

# # #np.save('/home/ybh10/CMIP6/{}/Historic/Numpy_Array/Ozone_Region/Ensamble/SO2_60_90S.npy'.format(model),SO2_SH)
# np.save('/home/ybh10/CMIP6/{}/Historic/Numpy_Array/Global/Ensamble/SO2_Global.npy'.format(model),SO2_Global_array)
# np.save('{}SO2_60_90S{}.npy'.format(file_output,end_sentence),np.nanmean(Ensamble_Regional_SO2,axis=0))
# np.save('{}{}SO2_60_90S{}.npy'.format(file_output,Ensamble_output,end_sentence),Ensamble_Regional_SO2) 
# np.save('{}SO2_Global{}.npy'.format(global_file_output,end_sentence),np.nanmean(Ensamble_Global_SO2,axis=0))
# np.save('{}{}SO2_Global{}.npy'.format(global_file_output,Ensamble_output,end_sentence),Ensamble_Global_SO2)                                   ######################################################################################

# ######################################################################################

# #np.save('/home/ybh10/CMIP6/{}/Historic/Numpy_Array/Ozone_Region/Ensamble/Ozone_Column_60_90S.npy'.format(model),Ozone_SH)
 #np.save('/home/ybh10/CMIP6/{}/Historic/Numpy_Array/Global/Ensamble/Ozone_Column_Global.npy'.format(model),Ozone_Global_array)
# np.save('{}Ozone_Column_60_90S{}.npy'.format(file_output,end_sentence),np.nanmean(Ensamble_Regional_Ozone,axis=0))
# np.save('{}{}Ozone_Column_60_90S{}.npy'.format(file_output,Ensamble_output,end_sentence),Ensamble_Regional_Ozone) 
# np.save('{}Ozone_Column_Global{}.npy'.format(global_file_output,end_sentence),np.nanmean(Ensamble_Global_Ozone,axis=0))
# np.save('{}{}Ozone_Column_Global{}.npy'.format(global_file_output,Ensamble_output,end_sentence),Ensamble_Global_Ozone) 

# #######################################################################################

# #np.save('/home/ybh10/CMIP6/{}/Historic/Numpy_Array/Ozone_Region/Ensamble/AOD_60_90S.npy'.format(model),AOD_SH)
np.save('/home/ybh10/CMIP6/{}/Historic/Numpy_Array/Global/Ensamble/AOD_Global.npy'.format(model),AOD_Global_array)
# np.save('{}AOD_60_90S{}.npy'.format(file_output,end_sentence),np.nanmean(Ensamble_Regional_AOD,axis=0))
# np.save('{}{}AOD_60_90S{}.npy'.format(file_output,Ensamble_output,end_sentence),Ensamble_Regional_AOD) 
np.save('{}AOD_Global{}.npy'.format(global_file_output,end_sentence),np.nanmean(Ensamble_Global_AOD,axis=0))
np.save('{}{}AOD_Global{}.npy'.format(global_file_output,Ensamble_output,end_sentence),Ensamble_Global_AOD)                      
# #######################################################################################

# #np.save('/home/ybh10/CMIP6/{}/Historic/Numpy_Array/Ozone_Region/Ensamble/Wind_60_90S.npy'.format(model),Wind_SH)
# np.save('/home/ybh10/CMIP6/{}/Historic/Numpy_Array/Global/Ensamble/Wind_Global.npy'.format(model),Wind_Global_array)
# np.save('{}Wind_60_90S{}.npy'.format(file_output,end_sentence),np.nanmean(Ensamble_Regional_Wind,axis=0))
# np.save('{}{}Wind_60_90S{}.npy'.format(file_output,Ensamble_output,end_sentence),Ensamble_Regional_Wind) 
# np.save('{}Wind_Global{}.npy'.format(global_file_output,end_sentence),np.nanmean(Ensamble_Global_Wind,axis=0))
# np.save('{}{}Wind_Global{}.npy'.format(global_file_output,Ensamble_output,end_sentence),Ensamble_Global_Wind)   

# #######################################################################################
#Oceanic_DMS
# # #np.save('/home/ybh10/CMIP6/{}/Historic/Numpy_Array/Ozone_Region/Ensamble/Oceanic_DMS_60_90S.npy'.format(model),ODMS_SH)
# np.save('/home/ybh10/CMIP6/{}/{}/Numpy_Array/Global/Ensamble/Oceanic_DMS_Global.npy'.format(model),ODMS_Global_array)      
# # np.save('{}Oceanic_DMS_60_90S{}.npy'.format(file_output,end_sentence),np.nanmean(Ensamble_Regional_Ocean_DMS,axis=0))
# # np.save('{}{}Oceanic_DMS_60_90S{}.npy'.format(file_output,Ensamble_output,end_sentence),Ensamble_Regional_Ocean_DMS)
# np.save('{}Oceanic_DMS_Global{}.npy'.format(global_file_output,end_sentence),np.nanmean(Ensamble_Global_Ocean_DMS,axis=0))
# np.save('{}{}Oceanic_DMS_Global{}.npy'.format(global_file_output,Ensamble_output,end_sentence),Ensamble_Global_Ocean_DMS)  

# #######################################################################################
#RF SW CS
# np.save('/home/ybh10/CMIP6/{}/Historic/Numpy_Array/Global/Ensamble/RFSW_CS_Global.npy',RFSW_Global_array)      
# np.save('{}RFSW_CS_Global{}.npy'.format(global_file_output,end_sentence),np.nanmean(Ensamble_Global_RFSW,axis=0))              
# np.save('{}{}RFSW_CS_Global{}.npy'.format(global_file_output,Ensamble_output,end_sentence),Ensamble_Global_RFSW) 

# #######################################################################################
#Surface Downward SW RF
# np.save('/home/ybh10/CMIP6/{}/Historic/Numpy_Array/Global/Ensamble/Surface_RFSW_CS_Global.npy'.format(model),Surface_RFSW_Global_array)      
# np.save('{}Surface_RFSW_CS_Global{}.npy'.format(global_file_output,end_sentence),np.nanmean(Ensamble_Global_Surface_RFSW,axis=0))              
# np.save('{}{}Surface_RFSW_CS_Global{}.npy'.format(global_file_output,Ensamble_output,end_sentence),Ensamble_Global_Surface_RFSW) 

# #######################################################################################
# TOA RF
np.save('/home/ybh10/CMIP6/{}/Historic/Numpy_Array/Global/Ensamble/TOA_RF_CS_Global.npy'.format(model),RF_Global_array)      
np.save('{}TOA_RF_CS_Global{}.npy'.format(global_file_output,end_sentence),np.nanmean(Ensamble_Global_RF,axis=0))              
np.save('{}{}TOA_RF_CS_Global{}.npy'.format(global_file_output,Ensamble_output,end_sentence),Ensamble_Global_RF) 
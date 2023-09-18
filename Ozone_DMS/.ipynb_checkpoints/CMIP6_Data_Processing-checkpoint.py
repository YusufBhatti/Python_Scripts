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

#from my_functions import *
os.chdir("/nesi/project/niwa02757/ybh10/CMIP6/UKESM1/Historic/")
Oz='Ozone_loss'
SO2='SO2'
DMS='DMS'
Historic='Historic'
ssp245='SSP245'
model='UKESM1'
dms_file = sorted(glob.glob('/nesi/project/niwa02757/ybh10/CMIP6/{}/Historic/Raw_Data/DMS/*.nc'.format(model)))
ODMS_file = sorted(glob.glob('/nesi/project/niwa02757/ybh10/CMIP6/Ocean_DMS/*{}*.nc'.format(model)))
# ODMS_file = sorted(glob.glob('/nesi/project/niwa02757/ybh10/CMIP6/Ocean_DMS/*{}*.nc'.format(model)))
# ODMS_file = sorted(glob.glob('/nesi/project/niwa02757/ybh10/CMIP6/Ocean_DMS/*{}*.nc'.format(model)))
# ODMS_file = sorted(glob.glob('/nesi/project/niwa02757/ybh10/CMIP6/Ocean_DMS/*{}*.nc'.format(model)))
# ODMS_file = sorted(glob.glob('/nesi/project/niwa02757/ybh10/CMIP6/Ocean_DMS/*{}*.nc'.format(model)))

files=[ODMS_file]

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

vn1=xr.open_dataset('/nesi/project/niwa02757/ybh10/CMIP6/UKESM1/Historic/Raw_Data/Ocean_DMS/dmsos_Omon_UKESM1-0-LL_historical_r10i1p1f2_gn_185001-194912.nc')
Olat=vn1.latitude
# Dates=['January', 'February','March','April','May','June','July','August','September','October','November','December']


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

def CMIP6(globe,dms_PPT,global_array):  # This function will give me all seasons when specified the first and last!
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
    


# PROCESSING THE RAW DATA IF NEW DATA IS AVAILABLE #######
########## 
## This code 

seasons='JJA' ## Choose the season you want, then the code and function will collate the correct 
month_data=seas(seasons) # information for the script below.
if seasons == 'DJF' !=-1: # The first = first month in slice ## IF NOT DJF
    First=month_data[0]
    first1=month_data[1] ### ONLY USE IF DJF - It should be 0
    Last=month_data[2] # The last month in the slice. Should be the last month.
else:
    First=month_data[0]
    Last=month_data[1]

end_sentence=''
file_output='/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Ozone_Region/{}/Ensamble/'.format(seasons)  
global_file_output='/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/{}/Ensamble/'.format(seasons)     

ssp_files=[ssp_dms_file,ssp_SO2_file,ssp_ODMS_file,ssp_O3_file,ssp_AOD_file,ssp_WIND_file]
monthy_1850=[]; monthy_1900=[]; monthy_1950=[]; monthy_2014=[]
monthy_reg_1850=[]; monthy_reg_1900=[]; monthy_reg_1950=[]; monthy_reg_2014=[]
name=['dms_','so2_','dmsos_Omon','o3_','toz','od550','dmsos_60','uas','rsutcs','rlutcs','rsdscs','rsdo','rss','rsds_','siconca','mmrso4','mmrss','cdnc','siconc_','chegpso4','loadss','loadso4']
    #DMS, SO2, Sea-DMS, Ozone, Column Ozone, AOD, N/A,         Wind,  TOA SWCS, TOA LWCS, Surface SW,


var=[]
var_reg=[]
ensamble=[]
ensamble_reg=[]

DMS_SH=[];     DMS_Global_array=[]; 
SO2_SH=[];     SO2_Global_array=[];
Ozone_SH=[];   Ozone_Global_array=[]
ODMS_SH=[];    ODMS_Global_array=[]
OZ_SH=[];      OZ_Global_array=[]
AOD_SH=[];     AOD_Global_array=[]
Wind_SH=[];    Wind_Global_array=[]
RFSW_SH=[];    RFSW_Global_array=[]
RSS_SH=[];     RSS_Global_array=[]
RSDS_SH=[];    RSDS_Global_array=[]
RSDO_SH=[];    RSDO_Global_array=[]
Sea_Ice_SH=[]; Sea_Ice_Global_array=[]
SSA_MMR_Global_array=[]
Sulfate_MMR_Global_array=[]
CDNC_Global_array=[]                
Ocean_Sea_Ice_Global_array=[]
CHEG_Global_array=[]                
SS_Global_array=[]                
SO4_Global_array=[]                

reg=[]; globe=[]; test=[]; test1=[]
Ensamble_Global_DMS=[]; Ensamble_Regional_DMS=[]
Ensamble_Global_SO2=[]; Ensamble_Regional_SO2=[]
Ensamble_Global_Ozone=[]; Ensamble_Regional_Ozone=[]
Ensamble_Global_AOD=[]; Ensamble_Regional_AOD=[]
Ensamble_Global_Wind=[]; Ensamble_Regional_Wind=[]
Ensamble_Global_Ocean_DMS=[]; Ensamble_Regional_Ocean_DMS=[]
Ensamble_Global_RFSW=[]; Ensamble_Regional_RFSW=[]
Ensamble_Global_Surface_RFSW=[]; Ensamble_Regional_Surface_RFSW=[]
Ensamble_Global_RSS=[]; Ensamble_Regional_RSS=[]
Ensamble_Global_RSDS=[]; Ensamble_Regional_RSDS=[]
Ensamble_Global_RSDO=[]; Ensamble_Regional_RSDO=[]
Ensamble_Global_Sea_Ice=[]; Ensamble_Regional_Sea_Ice=[]
Ensamble_Global_SSA_MMR=[]
Ensamble_Global_Sulfate_MMR=[]
Ensamble_Global_CDNC=[]
Ensamble_Global_Ocean_Sea_Ice=[]
test=[]; test1=[]; no_season_SH=[]; no_season_globe=[]
Ensamble_Global_CHEG=[]
Ensamble_Global_SS=[]
Ensamble_Global_SO4=[]

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
# #            no_season_SH.append(dms) ## These two lines just get ALL the months of the year
# #            no_season_globe.append(dms_globe) # They only need to be saved ONCE
# #             dmsx=CMIP6(dms,dms_globe,1,reg,globe)
#             test.append(dms_globe)
#             test1.append(reg)
#             globe=[]
#             reg=[]
#             hi=hi+1
#             print(hi)
#             if hi == 4 !=-1:
#                 total=np.concatenate((test[0],test[1],test[2],test[3]))
#                 total1=np.concatenate((test1[0],test1[1],test1[2],test1[3]))
# #                total2=np.concatenate((no_season_globe[0],no_season_globe[1],no_season_globe[2],no_season_globe[3]))
# #                total3=np.concatenate((no_season_SH[0],no_season_SH[1],no_season_SH[2],no_season_SH[3]))
#                 Ensamble_Global_DMS.append(total)
#                 Ensamble_Regional_DMS.append(total1)
# #               DMS_Global_array.append(total2)
# #                DMS_SH.append(total3)
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
# #            no_season_SH.append(so2)
# #            no_season_globe.append(so2_globe)
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
# #                total2=np.concatenate((no_season_globe[0],no_season_globe[1],no_season_globe[2],no_season_globe[3]))
# #                total3=np.concatenate((no_season_SH[0],no_season_SH[1],no_season_2_gn_18500101-19491230.nc
# #                SO2_Global_array.append(total2)
# #                SO2_SH.append(total3)
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
# #            no_season_SH.append(aod_mon)
# #            no_season_globe.append(aod_gl_mon)
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
# #                total2=np.concatenate((no_season_globe[0],no_season_globe[1]))
# #                total3=np.concatenate((no_season_SH[0],no_season_SH[1]))
#                 Ensamble_Global_Ozone.append(total)
#                 Ensamble_Regional_Ozone.append(total1)
# #                Ozone_Global_array.append(total2)
# #                Ozone_SH.append(total3)
#                 test1=[]
#                 test=[]
#                 no_season_SH=[]
#                 no_season_globe=[]
#                 hi=0             
#         if x.find (name[5]) !=-1:   #AOD
#             print(x)
#   #          print(f)
#             vn_o3=xr.open_dataset(x)
#             o3_globe=vn_o3.od550aer[:,:,:]
#             o3=vn_o3.od550aer[:,lat_inds[0]:lat_inds[(len(lat_inds)-1)]+1,:]
#             months=np.arange(0,len(o3_globe),30)
#             aod_mon=[]
#             aod_gl_mon=[]
#             for m in (months):
#                 aod_globe=np.nanmean(o3_globe[m:m+30,:,:],axis=0)
#                 aod=np.nanmean(o3[m:m+30,:,:],axis=0)
#                 aod_gl_mon.append(aod_globe)
#                 aod_mon.append(aod)
# #            no_season_SH.append(aod_mon)
# #            no_season_globe.append(aod_gl_mon)
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
# #                total2=np.concatenate((no_season_globe[0],no_season_globe[1]))
# #                total3=np.concatenate((no_season_SH[0],no_season_SH[1]))
#                 Ensamble_Global_AOD.append(total)
#                 Ensamble_Regional_AOD.append(total1)
# #                AOD_Global_array.append(total2)
# #                AOD_SH.append(total3)
#                 test1=[]
#                 test=[]
#                 no_season_SH=[]
#                 no_season_globe=[]
#                 hi=0  
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
# #            no_season_SH.append(aod_mon)
# #            no_season_globe.append(aod_gl_mon)
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
# #                total2=np.concatenate((no_season_globe[0],no_season_globe[1]))
# #                total3=np.concatenate((no_season_SH[0],no_season_SH[1]))
#                 Ensamble_Global_Wind.append(total)
#                 Ensamble_Regional_Wind.append(total1)
# #                Wind_Global_array.append(total2)
# #                Wind_SH.append(total3)
#                 test1=[]
#                 test=[]
#                 no_season_SH=[]
#                 no_season_globe=[]
#                 hi=0 
#         if x.find (name[2]) !=-1: # FOR OCEANIC DMS --> NEMO MODEL 
#             vn_o3=xr.open_dataset(x)
#             abc=vn_o3.dmsos[:]
#             abc_reg=vn_o3.dmsos[:,0:94,:] # Requires changing IF changing latitude constraint
# #            no_season_SH.append(abc_reg)
# #            no_season_globe.append(abc)
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
# #                total2=np.concatenate((no_season_globe[0],no_season_globe[1]))
# #                total3=np.concatenate((no_season_SH[0],no_season_SH[1]))
#                 Ensamble_Global_Ocean_DMS.append(total)
#                 Ensamble_Regional_Ocean_DMS.append(total1)
# #                ODMS_Global_array.append(total2)
# #                ODMS_SH.append(total3)
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
#         if x.find (name[10]) !=-1:   #Surface SW good
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
#         if x.find (name[11]) !=-1:  # Seawater SW  GOOOOOD
#             print(x)
#             vn_dms=xr.open_dataset(x)
#             dms_globe=vn_dms.rsdo[:,0,:,:]
#             no_season_globe.append(dms_globe) # They only need to be saved ONCE
#             dmsx=CMIP6(dms_globe,dms_globe,1,reg,globe)
#             test.append(globe)
#             globe=[]
#             reg=[]
#             hi=hi+1
#             print(hi)
#             if hi == 4 !=-1:
#                 total=np.concatenate((test[0],test[1],test[2],test[3]))
#                 total2=np.concatenate((no_season_globe[0],no_season_globe[1],no_season_globe[2],no_season_globe[3]))
#                 Ensamble_Global_RSDO.append(total)
#                 RSDO_Global_array.append(total2)
#                 test1=[]
#                 test=[]
#                 no_season_SH=[]
#                 no_season_globe=[]
#                 hi=0
#         if x.find (name[12]) !=-1: # Surface Net Down SW flux 
#             vn_o3=xr.open_dataset(x)
#             abc=vn_o3.rss[:]
#             no_season_globe.append(abc)
#             dmsx=CMIP6(abc,abc,1,reg,globe)
#             test.append(globe)
#             globe=[]
#             reg=[]
#             hi=hi+1
#             print(hi)
#             if hi == 2 !=-1:
#                 total=np.concatenate((test[0],test[1]))
#                 total2=np.concatenate((no_season_globe[0],no_season_globe[1]))
#                 Ensamble_Global_RSS.append(total)
#                 ODMS_Global_array.append(total2)
#                 test1=[]
#                 test=[]
#                 no_season_SH=[]
#                 no_season_globe=[]
#                 hi=0    
#         if x.find (name[13]) !=-1:   #Surface down SW in air all sky
#             print(x)
#             vn_o3=xr.open_dataset(x)
#             o3_globe=vn_o3.rsds[:]
#             months=np.arange(0,len(o3_globe),30)
#             aod_mon=[]
#             aod_gl_mon=[]
#             for m in (months):
#                 aod_globe=np.nanmean(o3_globe[m:m+30,:,:],axis=0)
#                 aod_gl_mon.append(aod_globe)
#             no_season_globe.append(aod_gl_mon)
#             dmsx=CMIP6(aod_gl_mon,1,globe)
#             test.append(globe)
#             globe=[]
#             reg=[]
#             hi=hi+1
#             print(hi)
#             if hi == 2 !=-1:
#                 total=np.concatenate((test[0],test[1]))
#                 total2=np.concatenate((no_season_globe[0],no_season_globe[1]))
#                 Ensamble_Global_RSDS.append(total)
#                 RSDS_Global_array.append(total2)
#                 test1=[]
#                 test=[]
#                 no_season_SH=[]
#                 no_season_globe=[]
#                 hi=0 
#         if x.find (name[14]) !=-1:   #Sea Ice Concentration
#             print(x)
#             vn_o3=xr.open_dataset(x)
#             o3_globe=vn_o3.siconca[:]
#         #    o3=vn_o3.uas[:,lat_inds[0]:lat_inds[(len(lat_inds)-1)]+1,:]
#             months=np.arange(0,len(o3_globe),30)
#             no_season_globe.append(o3_globe)
#             dmsx=CMIP6(o3_globe,1,globe)
#             test.append(globe)
#          #   test1.append(reg)
#             globe=[]
#             reg=[]
#             hi=hi+1
#             print(hi)
#             if hi == 2 !=-1:
#                 total=np.concatenate((test[0],test[1]))
#                 #total1=np.concatenate((test1[0],test1[1]))
#                 total2=np.concatenate((no_season_globe[0],no_season_globe[1]))
#                # total3=np.concatenate((no_season_SH[0],no_season_SH[1]))
#                 Ensamble_Global_Sea_Ice.append(total)
#               #  Ensamble_Regional_Wind.append(total1)
#                 Sea_Ice_Global_array.append(total2)
#               #  Wind_SH.append(total3)
#                 test1=[]
#                 test=[]
#                 no_season_SH=[]
#                 no_season_globe=[]
#                 hi=0 
#         if x.find (name[15]) !=-1: # Sulfate MMR
#             print(x)
#             vn_so2=xr.open_dataset(x)
#             so2_globe=vn_so2.mmrso4[:,0,:,:]
#  #           no_season_globe.append(so2_globe)
#             dmsx=CMIP6(so2_globe,1,globe)
#             test.append(globe)
#             globe=[]
#             reg=[]
#             hi=hi+1
#             print(hi)
#             if hi == 4 !=-1:
#                 total=np.concatenate((test[0],test[1],test[2],test[3]))
#                 Ensamble_Global_Sulfate_MMR.append(total)
# #                total2=np.concatenate((no_season_globe[0],no_season_globe[1],no_season_globe[2],no_season_globe[3]))
# #                Sulfate_MMR_Global_array.append(total2)
#                 test1=[]
#                 test=[]
#                 no_season_SH=[]
#                 no_season_globe=[]
#                 hi=0                
#         if x.find (name[16]) !=-1: # SSA MMR
#             print(x)
#             vn_so2=xr.open_dataset(x)
#             so2_globe=vn_so2.mmrss[:,0,:,:]
#         #    no_season_globe.append(so2_globe)
#             dmsx=CMIP6(so2_globe,1,globe)
#             test.append(globe)
#             globe=[]
#             reg=[]
#             hi=hi+1
#             print(hi)
#             if hi == 4 !=-1:
#                 total=np.concatenate((test[0],test[1],test[2],test[3]))
#                 Ensamble_Global_SSA_MMR.append(total)
# #                 total2=np.concatenate((no_season_globe[0],no_season_globe[1],no_season_globe[2],no_season_globe[3]))
# #                 SSA_MMR_Global_array.append(total2)
#                 test1=[]
#                 test=[]
#                 no_season_SH=[]
#                 no_season_globe=[]
#                 hi=0                      
#         if x.find (name[17]) !=-1:  # CDNC
#             print(x)
#             vn_dms=xr.open_dataset(x)
#             dms_globe=vn_dms.cdnc[:,7,:,:]
# #            no_season_globe.append(dms_globe) # They only need to be saved ONCE
#             dmsx=CMIP6(dms_globe,1,globe)
#             test.append(globe)
#             globe=[]
#             hi=hi+1
#             print(hi)
#             if hi == 4 !=-1:
#                 total=np.concatenate((test[0],test[1],test[2],test[3]))
#                # total2=np.concatenate((no_season_globe[0],no_season_globe[1],no_season_globe[2],no_season_globe[3]))
#                 Ensamble_Global_DMS.append(total)
# #                CDNC_Global_array.append(total2)
#                 test=[]
#                 no_season_globe=[]
#                 hi=0
#         if x.find (name[18]) !=-1:   #Sea Ice Concentration
#             print(x)
#             vn_o3=xr.open_dataset(x)
#             o3_globe=vn_o3.siconc[:]
#             months=np.arange(0,len(o3_globe),30)
#             no_season_globe.append(o3_globe)
#             dmsx=CMIP6(o3_globe,1,globe)
#             test.append(globe)
#             globe=[]
#             reg=[]
#             hi=hi+1
#             print(hi)
#             if hi == 2 !=-1:
#                 total=np.concatenate((test[0],test[1]))
#                 total2=np.concatenate((no_season_globe[0],no_season_globe[1]))
#                 Ensamble_Global_Ocean_Sea_Ice.append(total)
#                 Ocean_Sea_Ice_Global_array.append(total2)
#                 test=[]
#                 no_season_globe=[]
#                 hi=0 
        if x.find (name[19]) !=-1:  # cheg
            print(x)
            vn_dms=xr.open_dataset(x)
            dms_globe=vn_dms.chegpso4[:,0,:,:]
            no_season_globe.append(dms_globe) # They only need to be saved ONCE
            dmsx=CMIP6(dms_globe,1,globe)
            test.append(globe)
            globe=[]
            hi=hi+1
            print(hi)
            if hi == 4 !=-1:
                total=np.concatenate((test[0],test[1],test[2],test[3]))
                total2=np.concatenate((no_season_globe[0],no_season_globe[1],no_season_globe[2],no_season_globe[3]))
                Ensamble_Global_CHEG.append(total)
                CHEG_Global_array.append(total2)
                test=[]
                no_season_globe=[]
                hi=0
#                 hi=0 
        if x.find (name[20]) !=-1:   #SS
            print(x)
            vn_o3=xr.open_dataset(x)
            o3_globe=vn_o3.loadss[:]
            no_season_globe.append(o3_globe)
            dmsx=CMIP6(o3_globe,1,globe)
            test.append(globe)
            globe=[]
            reg=[]
            hi=hi+1
            print(hi)
            if hi == 2 !=-1:
                total=np.concatenate((test[0],test[1]))
                total2=np.concatenate((no_season_globe[0],no_season_globe[1]))
                Ensamble_Global_SS.append(total)
                SS_Global_array.append(total2)
                test=[]
                no_season_globe=[]
                hi=0 
        if x.find (name[21]) !=-1:   #SO4
            print(x)
            vn_o3=xr.open_dataset(x)
            o3_globe=vn_o3.loadso4[:]
            no_season_globe.append(o3_globe)
            dmsx=CMIP6(o3_globe,1,globe)
            test.append(globe)
            globe=[]
            reg=[]
            hi=hi+1
            print(hi)
            if hi == 2 !=-1:
                total=np.concatenate((test[0],test[1]))
                total2=np.concatenate((no_season_globe[0],no_season_globe[1]))
                Ensamble_Global_SO4.append(total)
                SO4_Global_array.append(total2)
                test=[]
                no_season_globe=[]
                hi=0 
# sw=sorted(glob.glob('/home/ybh10/CMIP6/UKESM1/Historic/Radiative_Forcing/Clear_Sky/SW/*'))
# lw=sorted(glob.glob('/home/ybh10/CMIP6/UKESM1/Historic/Radiative_Forcing/Clear_Sky/LW/*'))
# hi=0
# rf_1850=[]
# globe=[]
# reg=[]
# test=[]; test1=[]; no_season_SH=[]; no_season_globe=[]
# Ensamble_Global_RF=[]; Ensamble_Regional_RF=[];RF_Global_array=[]
# for s,l in zip(sw,lw):
#     hi=hi+1
#     print(hi)
#     if s.find ('rsutcs_C') !=-1:
#         vn_rf=xr.open_dataset(s)
#         swrf=vn_rf.rsutcs
#         print(s)
#     if l.find ('rlutcs_C') !=-1:
#         vn_rf=xr.open_dataset(l)
#         lwrf=vn_rf.rlutcs
#         rf=swrf+lwrf
#         print(l)

#         months=np.arange(0,len(rf),30)
#         aod_mon=[]
#         aod_gl_mon=[]
#         for m in (months):
#             aod_globe=np.nanmean(rf[m:m+30,:,:],axis=0)
#       #      aod=np.nanmean(o3[m:m+30,:,:],axis=0)
#             aod_gl_mon.append(aod_globe)
#        #     aod_mon.append(aod)
#             #rf_1850.append(swrf+lwrf)
#         print('--------------------------------')
#      #   no_season_globe.append(aod_gl_mon)
#         dmsx=CMIP6(aod_gl_mon,aod_gl_mon,1,reg,globe)
#         test.append(globe)
#     #    test1.append(reg)
#         globe=[]
#         reg=[]
#     if hi == 2 !=-1:
#         total=np.concatenate((test[0],test[1]))
#      #   total1=np.concatenate((test1[0],test1[1]))
#    #     total2=np.concatenate((no_season_globe[0],no_season_globe[1]))
# #                total3=np.concatenate((no_season_SH[0],no_season_SH[1]))
#         Ensamble_Global_RF.append(total)
#         #Ensamble_Regional_RF.append(total1)
#         RF_Global_array.append(total2)
# #                Wind_SH.append(total3)
#      #   test1=[]
#         test=[]
#       #  no_season_SH=[]
#         no_season_globe=[]
#         hi=0 

# Saves ALL arrays into the correct directory. Run this cell after each time you change the season.
# The first two lines for each variable will only need to be run ONCE as it saves ALL the months in the 165 year timescale.
#seasons='SON' # Already defined above.
end_sentence=''
file_output='/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Ozone_Region/{}/'.format(seasons)  
global_file_output='/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/{}/'.format(seasons)     
Ensamble_output='Ensamble/'

# #np.save('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Ozone_Region/Ensamble/DMS_60_90S.npy',DMS_SH)
# #np.save('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/Ensamble/DMS_Global.npy',DMS_Global_array)
# np.save('{}DMS_60_90S{}.npy'.format(file_output,end_sentence),np.nanmean(Ensamble_Regional_DMS,axis=0))
# np.save('{}{}DMS_60_90S{}.npy'.format(file_output,Ensamble_output,end_sentence),Ensamble_Regional_DMS)
# np.save('{}DMS_Global{}.npy'.format(global_file_output,end_sentence),np.nanmean(Ensamble_Global_DMS,axis=0))              
# np.save('{}{}DMS_Global{}.npy'.format(global_file_output,Ensamble_output,end_sentence),Ensamble_Global_DMS)              
# ######################################################################################

# # #np.save('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Ozone_Region/Ensamble/SO2_60_90S.npy',SO2_SH)
# # #np.save('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/Ensamble/SO2_Global.npy',SO2_Global_array)
# np.save('{}SO2_60_90S{}.npy'.format(file_output,end_sentence),np.nanmean(Ensamble_Regional_SO2,axis=0))
# np.save('{}{}SO2_60_90S{}.npy'.format(file_output,Ensamble_output,end_sentence),Ensamble_Regional_SO2) 
# np.save('{}SO2_Global{}.npy'.format(global_file_output,end_sentence),np.nanmean(Ensamble_Global_SO2,axis=0))
# np.save('{}{}SO2_Global{}.npy'.format(global_file_output,Ensamble_output,end_sentence),Ensamble_Global_SO2)                                   ######################################################################################

# ######################################################################################

# #np.save('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Ozone_Region/Ensamble/Ozone_Column_60_90S.npy',Ozone_SH)
# #np.save('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/Ensamble/Ozone_Column_Global.npy',Ozone_Global_array)
# np.save('{}Ozone_Column_60_90S{}.npy'.format(file_output,end_sentence),np.nanmean(Ensamble_Regional_Ozone,axis=0))
# np.save('{}{}Ozone_Column_60_90S{}.npy'.format(file_output,Ensamble_output,end_sentence),Ensamble_Regional_Ozone) 
# np.save('{}Ozone_Column_Global{}.npy'.format(global_file_output,end_sentence),np.nanmean(Ensamble_Global_Ozone,axis=0))
# np.save('{}{}Ozone_Column_Global{}.npy'.format(global_file_output,Ensamble_output,end_sentence),Ensamble_Global_Ozone) 

# #######################################################################################

# #np.save('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Ozone_Region/Ensamble/AOD_60_90S.npy',AOD_SH)
# #np.save('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/Ensamble/AOD_Global.npy',AOD_Global_array)
# np.save('{}AOD_60_90S{}.npy'.format(file_output,end_sentence),np.nanmean(Ensamble_Regional_AOD,axis=0))
# np.save('{}{}AOD_60_90S{}.npy'.format(file_output,Ensamble_output,end_sentence),Ensamble_Regional_AOD) 
# np.save('{}AOD_Global{}.npy'.format(global_file_output,end_sentence),np.nanmean(Ensamble_Global_AOD,axis=0))
# np.save('{}{}AOD_Global{}.npy'.format(global_file_output,Ensamble_output,end_sentence),Ensamble_Global_AOD)                      
# #######################################################################################

# #np.save('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Ozone_Region/Ensamble/Wind_60_90S.npy',Wind_SH)
# #np.save('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/Ensamble/Wind_Global.npy',Wind_Global_array)
# np.save('{}Wind_60_90S{}.npy'.format(file_output,end_sentence),np.nanmean(Ensamble_Regional_Wind,axis=0))
# np.save('{}{}Wind_60_90S{}.npy'.format(file_output,Ensamble_output,end_sentence),Ensamble_Regional_Wind) 
# np.save('{}Wind_Global{}.npy'.format(global_file_output,end_sentence),np.nanmean(Ensamble_Global_Wind,axis=0))
# np.save('{}{}Wind_Global{}.npy'.format(global_file_output,Ensamble_output,end_sentence),Ensamble_Global_Wind)   

# #######################################################################################
#Oceanic_DMS
# #np.save('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Ozone_Region/Ensamble/Oceanic_DMS_60_90S.npy',ODMS_SH)
# #np.save('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/Ensamble/Oceanic_DMS_Global.npy',ODMS_Global_array)      
# np.save('{}Oceanic_DMS_60_90S{}.npy'.format(file_output,end_sentence),np.nanmean(Ensamble_Regional_Ocean_DMS,axis=0))
# np.save('{}{}Oceanic_DMS_60_90S{}.npy'.format(file_output,Ensamble_output,end_sentence),Ensamble_Regional_Ocean_DMS)
# np.save('{}Oceanic_DMS_Global{}.npy'.format(global_file_output,end_sentence),np.nanmean(Ensamble_Global_Ocean_DMS,axis=0))
# np.save('{}{}Oceanic_DMS_Global{}.npy'.format(global_file_output,Ensamble_output,end_sentence),Ensamble_Global_Ocean_DMS)  

# #######################################################################################
#RF SW CS
# np.save('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/Ensamble/RFSW_CS_Global.npy',RFSW_Global_array)      
# np.save('{}RFSW_CS_Global{}.npy'.format(global_file_output,end_sentence),np.nanmean(Ensamble_Global_RFSW,axis=0))              
# np.save('{}{}RFSW_CS_Global{}.npy'.format(global_file_output,Ensamble_output,end_sentence),Ensamble_Global_RFSW) 

# # #######################################################################################
# #Surface Downward SW RF
# np.save('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/Ensamble/Surface_RFSW_CS_Global.npy',Surface_RFSW_Global_array)      
# np.save('{}Surface_RFSW_CS_Global{}.npy'.format(global_file_output,end_sentence),np.nanmean(Ensamble_Global_Surface_RFSW,axis=0))              
# np.save('{}{}Surface_RFSW_CS_Global{}.npy'.format(global_file_output,Ensamble_output,end_sentence),Ensamble_Global_Surface_RFSW)

# # #######################################################################################
# #Seawater Downward SW
# np.save('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/Ensamble/RSDO_Global.npy',RSDO_Global_array)      
# np.save('{}RSDO_Global{}.npy'.format(global_file_output,end_sentence),np.nanmean(Ensamble_Global_RSDO,axis=0))              
# np.save('{}{}RSDO_Global{}.npy'.format(global_file_output,Ensamble_output,end_sentence),Ensamble_Global_RSDO)

# # #######################################################################################
# #Surface Net Down SW Flux
# np.save('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/Ensamble/RSS_Global.npy',RSS_Global_array)      
# np.save('{}RSS_Global{}.npy'.format(global_file_output,end_sentence),np.nanmean(Ensamble_Global_RSS,axis=0))              
# np.save('{}{}RSS_Global{}.npy'.format(global_file_output,Ensamble_output,end_sentence),Ensamble_Global_RSS)

# # #######################################################################################
# #Surface Downward SW in air All Sky
# np.save('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/Ensamble/RSDS_Global.npy',RSDS_Global_array)      
# np.save('{}RSDS_Global{}.npy'.format(global_file_output,end_sentence),np.nanmean(Ensamble_Global_RSDS,axis=0))              
# np.save('{}{}RSDS_Global{}.npy'.format(global_file_output,Ensamble_output,end_sentence),Ensamble_Global_RSDS)

# #######################################################################################
# Sea Ice Area In Atmosphere Model
#np.save('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/Ensamble/Sea_Ice_Global.npy',Sea_Ice_Global_array)      
# np.save('{}Sea_Ice_Global{}.npy'.format(global_file_output,end_sentence),np.nanmean(Ensamble_Global_Sea_Ice,axis=0))
# np.save('{}{}Sea_Ice_Global{}.npy'.format(global_file_output,Ensamble_output,end_sentence),Ensamble_Global_Sea_Ice)

# #######################################################################################
# # Sulfate MMR 
#np.save('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/Ensamble/Sulfate_MMR_Global.npy',Sulfate_MMR_Global_array)   
# np.save('{}Sulfate_MMR_Global{}.npy'.format(global_file_output,end_sentence),np.nanmean(Ensamble_Global_Sulfate_MMR,axis=0))
# np.save('{}{}Sulfate_MMR_Global{}.npy'.format(global_file_output,Ensamble_output,end_sentence),Ensamble_Global_Sulfate_MMR)

# # #######################################################################################
# # SSA MMR 
# #np.save('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/Ensamble/SSA_MMR_Global.npy',SSA_MMR_Global_array)   
# np.save('{}SSA_MMR_Global{}.npy'.format(global_file_output,end_sentence),np.nanmean(Ensamble_Global_SSA_MMR,axis=0))
# np.save('{}{}SSA_MMR_Global{}.npy'.format(global_file_output,Ensamble_output,end_sentence),Ensamble_Global_SSA_MMR)
# # #######################################################################################
# # CDNC 
#np.save('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/Ensamble/CDNC_Global.npy',CDNC_Global_array)   
# np.save('{}CDNC_Global{}.npy'.format(global_file_output,end_sentence),np.nanmean(Ensamble_Global_CDNC,axis=0))
# np.save('{}{}CDNC_Global{}.npy'.format(global_file_output,Ensamble_output,end_sentence),Ensamble_Global_CDNC)

# #######################################################################################
# Sea Ice Area In OCEAN Model
# np.save('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/Ensamble/Ocean_Sea_Ice_Global.npy',Ocean_Sea_Ice_Global_array)      
# np.save('{}Ocean_Sea_Ice_Global{}.npy'.format(global_file_output,end_sentence),np.nanmean(Ensamble_Global_Ocean_Sea_Ice,axis=0))
# np.save('{}{}Ocean_Sea_Ice_Global{}.npy'.format(global_file_output,Ensamble_output,end_sentence),Ensamble_Global_Ocean_Sea_Ice)
# #######################################################################################
# CHEG
np.save('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/Ensamble/CHEG_Global.npy',CHEG_Global_array)      
np.save('{}CHEG_Global{}.npy'.format(global_file_output,end_sentence),np.nanmean(Ensamble_Global_CHEG,axis=0))
np.save('{}{}CHEG_Global{}.npy'.format(global_file_output,Ensamble_output,end_sentence),Ensamble_Global_CHEG)
# #######################################################################################
# SS
np.save('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/Ensamble/SS_Global.npy',SS_Global_array)      
np.save('{}SS_Global{}.npy'.format(global_file_output,end_sentence),np.nanmean(Ensamble_Global_SS,axis=0))
np.save('{}{}SS_Global{}.npy'.format(global_file_output,Ensamble_output,end_sentence),Ensamble_Global_SS)
# #######################################################################################
# SO4
np.save('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/Ensamble/SO4_Global.npy',SO4_Global_array)      
np.save('{}SO4_Global{}.npy'.format(global_file_output,end_sentence),np.nanmean(Ensamble_Global_SO4,axis=0))
np.save('{}{}SO4_Global{}.npy'.format(global_file_output,Ensamble_output,end_sentence),Ensamble_Global_SO4)
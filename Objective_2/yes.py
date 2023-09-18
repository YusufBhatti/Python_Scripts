#!/usr/bin/env python
# coding: utf-8

# In[ ]:


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
from scipy.stats import sem
import gc
AOD=[]
import xarray as xr
from my_functions import *
os.chdir("/nesi/project/niwa02757/ybh10/CMIP6/UKESM1/Historic/")

ODMS_file = sorted(glob.glob('/nesi/project/niwa02757/ybh10/CMIP6/UKESM1/Historic/Raw_Data/Ocean_DMS/*.nc'))
# O3_file = sorted(glob.glob('/nesi/project/niwa02757/ybh10/CMIP6/UKESM1/Historic/Ozone_loss//*.nc'))

months = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
vn2=xr.open_dataset(ODMS_file[3])
olat=vn2.latitude
olon=vn2.longitude

vn1=xr.open_dataset('/nesi/project/niwa02757/ybh10/CMIP6/UKESM1/Historic/Raw_Data/DMS/dms_AERmon_UKESM1-0-LL_historical_r10i1p1f2_gn_185001-189912.nc')
lats=vn1.lat
lons=vn1.lon
time=vn1.time
# Dates=['January', 'February','March','April','May','June','July','August','September','October','November','December']
lat_bnds = [-90, -60]
inpath_modis=sorted(glob.glob('/nesi/nobackup/niwa02757/ybh10/Observational_Data/MODIS/AOD/*.nc'))
latty_o=xr.open_dataset(inpath_modis[0])
lat_o=latty_o['lat']
lon_o=latty_o['lon']
vn1=xr.open_dataset('/nesi/project/niwa02757/ybh10/CMIP6/UKESM1/Historic/Raw_Data/Ocean_DMS/dmsos_Omon_UKESM1-0-LL_historical_r10i1p1f2_gn_185001-194912.nc')
Olat=vn1.latitude

DMS_ppt=((29/62.13)*1e12)
ODMS_ppm=1e6
ODMS_ppt=(((18/62.13)*1e12))

SO2_ppt=((29/64.06)*1e12)
O3_ppm=((29/48)*1e6)
season=['DJF','MAM','JJA','SON']
 
lat=lats[19:32] # 50 - 65s IN UKESM1:
files=[]
files_global=[]
AOD_ozone=[]; DMS_ozone=[]; Oceanic_DMS_ozone=[]; Ozone_Column_ozone=[]; RF_CS_ozone=[]; SO2_ozone=[];Surface_RF_CS_ozone=[]; Wind_ozone=[];RFSW_CS_ozone=[]
variables=[AOD_ozone,DMS_ozone,Oceanic_DMS_ozone,Ozone_Column_ozone,RF_CS_ozone,SO2_ozone,Surface_RF_CS_ozone,Wind_ozone]

DJF_var=[]; MAM_var=[]; JJA_var=[]; SON_var=[]
season_var=[DJF_var,MAM_var,JJA_var,SON_var]
DJF_global_var=[]; MAM_global_var=[]; JJA_global_var=[]; SON_global_var=[]
season_global_var=[DJF_global_var,MAM_global_var,JJA_global_var,SON_global_var]
last=[]
#ad=0; dm=1; od=2; oz=3; r=4; s=5; w=6; sw=7; sdsw=8
ad=0; dm=1; od=2; oz=3;  sw=4; r=5; RSDO=6; RSDS=7; RSS=8;  s=9;SSA_H2SO4=10;SSA=11; seaice=12;SO4MMR=13; RSDSCS=14; w=15

names=['AOD_','DMS_','Oceanic_DMS','Ozone','RF_CS','SO2','Wind','RFSW','Surface_RFSW']
for s,glvar in zip(season,season_global_var):
    Globa=sorted(glob.glob("/nesi/project/niwa02757/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/{}/Ensamble/*".format(s)))
    files_global.append(Globa)
    for g in (Globa):
#         if g.find (Globa[4]) !=-1: #RFSW_CS_Global
#             second_last=(np.load(g))
#         else:
#             if g.find (Globa[7]) !=-1: # Surface_RFSW
#                 last=(np.load(g))
        
        glo_filey=np.load(g)
        #print(g)
        glvar.append(glo_filey)
#             if g.find ('/Ensamble/Wind_Global.npy') !=-1:
#               #  print('end of season')
#                 glvar.append(second_last)
#                 glvar.append(last)
                    
                
#season_global_var.append(last)

# for s,var in zip(season,season_var):
#     Regi=sorted(glob.glob("/nesi/project/niwa02757/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Ozone_Region/{}/Ensamble/*".format(s)))
#     files.append(Regi)
#     for f in (Regi):
#         reg_filey=np.load(f)
#         var.append(reg_filey)
# #         glo_filey=np.load(g)
# #         glvar.append(glo_filey)
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


# In[3]:


# numbers=[ad,dm,od,oz,RSDO,RSDS,RSS,so2,SSA,seaice,SO4MMR,RSDSCS,w]

ad=0; dm=1; od=2; oz=3;  sw=4; r=5; RSDO=6; RSDS=7; RSS=8;  so2=9;SSA_H2SO4=10;SSA=11; seaice=12;SO4MMR=13; RSDSCS=14; w=15


# #numbers=[ad,dm,od,oz,RSDO,RSDS,RSS,s,SSA,seaice,SO4MMR,RSDSCS,w]
numbers=[so2]
maybe=[]
# large_data=[SON_global_var,SON_global_var,SON_global_var,SON_global_var,SON_global_var,SON_global_var,
#            SON_global_var,SON_global_var,SON_global_var,SON_global_var,SON_global_var,SON_global_var,SON_global_var]
# #large_data=[SON_global_var,SON_global_var,SON_global_var,SON_global_var,SON_global_var]
# olaty=np.nanmean(olat[0:138],axis=1)
# ocean_sw=large_data[1][RSDO]
# sea=np.nanmean(ocean_sw[:,100:165],axis=(0))
# #sea=sea.transpose(1,0)
# old_sea=np.nanmean(ocean_sw[:,90:111,:,:],axis=(0,1))
# sea_sw=np.empty((65,330,360)); sea_sw[:]=np.nan
# for i in range(0,65):
#     for n in range(0,330):
#         for l in range(0,360):
#             sea_sw[i,n,l]=((sea[i,n,l]/old_sea[n,l])-1)*100
# sea_sw=np.array(sea_sw)
# #sea_sw=sea_sw.transpose(1,0)

# ocean_dms=large_data[2][od]
# sea=np.nanmean(ocean_dms[:,100:165],axis=(0))
# #sea=sea.transpose(1,0)
# old_sea=np.nanmean(ocean_dms[:,90:111,:,:],axis=(0,1))
# sea_dms=np.empty((65,330,360)); sea_dms[:]=np.nan
# for i in range(0,65):
#     for n in range(0,330):
#         for l in range(0,360):
#             sea_dms[i,n,l]=((sea[i,n,l]/old_sea[n,l])-1)*100
# sea_dms=np.array(sea_dms)
# #sea_dms=sea_dms.transpose(1,0)
# #sea_dms=sea_dms*ODMS_ppm

# #large_data=[SON_global_var,SON_global_var,SON_global_var,sea_sw,sea_dms]
large_data=[DJF_global_var]
name='SSA'
# large_data=[SON_global_var,SON_global_var,sea_dms,SON_global_var,sea_sw,SON_global_var,SON_global_var,SON_global_var,
#            SON_global_var,SON_global_var,SON_global_var,SON_global_var,SON_global_var]
# #numbers=[oz,RSS,od,w,dm,ad]#,r] # ad chosen for ODMS due to the array having a 'one' at the start - easy to put in loop
# # R2_Abs_Diff=[]
#for e in range(0,13):

for num,datas in zip(numbers,large_data):
    print(num)
    if num == RSDO !=-1:
        maybe.append(sea_sw)
    else:
        if num == od !=-1:
            maybe.append(sea_dms)
            
        else:
            old_vari=np.nanmean(datas[num][:,0:91,:,:],axis=(1))
            vari=np.empty((13, 65, 144,192)); vari[:]=np.nan
          #  data=np.nanmean(datas[num],axis=(3))
            for t in range(0,datas[num].shape[0]):
                for i in range(0,vari.shape[1]):
                    for la in range(0,datas[num].shape[2]):
                        for lo in range(0,datas[num].shape[3]):
                        #for x in range(0,datas[num].shape[3]):
                              #  print('hi')
                            if num == w !=-1:
                                vari[t,i,la,lo]=(datas[num][t,i+100,la,lo]-old_vari[t,la,lo])
                            else:
                                vari[t,i,la,lo]=((datas[num][t,i+100,la,lo]/old_vari[t,la,lo])-1)*100
            print(np.shape(vari))
            print(np.nanmean(vari))
            see=np.nanmean(vari,axis=(0))
            maybe.append(see)
            maybe=np.array(maybe)
          #  see=see.transpose(1,0)      
    
gc.collect()

file = '/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/Difference/Rel_Diff/'
np.save('{}DJF_{}.npy'.format(file,name),maybe[0])
maybe=[]

for num,datas in zip(numbers,large_data):
    print(num)
    if num == RSDO !=-1:
        maybe.append(sea_sw)
    else:
        if num == od !=-1:
            maybe.append(sea_dms)
            
        else:
            old_vari=np.nanmean(datas[num][:,90:111,:,:],axis=(1))
            vari=np.empty((13, 65, 144,192)); vari[:]=np.nan
          #  data=np.nanmean(datas[num],axis=(3))
            for t in range(0,datas[num].shape[0]):
                for i in range(0,vari.shape[1]):
                    for la in range(0,datas[num].shape[2]):
                        for lo in range(0,datas[num].shape[3]):
                        #for x in range(0,datas[num].shape[3]):
                              #  print('hi')
                            if num == w !=-1:
                                vari[t,i,la,lo]=(datas[num][t,i+100,la,lo]-old_vari[t,la,lo])
                            else:
                                vari[t,i,la,lo]=((datas[num][t,i+100,la,lo]/old_vari[t,la,lo])-1)*100
            print(np.shape(vari))
            print(np.nanmean(vari))
            see=np.nanmean(vari,axis=(0))
            maybe.append(see)
            maybe=np.array(maybe)

          #  see=see.transpose(1,0)      
    
gc.collect()

file = '/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/Difference/1940_1960_mean/Rel_Diff/DJF/'
np.save('{}{}.npy'.format(file,name),maybe[0])

# saves=['AOD','DMS','ODMS','Ozone','RSDO','RSDS','RSS','SO2','SSA','Sea_Ice','H2SO4MMR','RSDSCS','Wind']
# os.chdir(file+'/Rel_Diff/SON/')
# for i,names in zip(range(0,len(maybe)),saves):
#     np.save((names),(maybe[i]))
# large_data=[SON_global_var]
# maybe=[]

# for num,datas in zip(numbers,large_data):
#     print(num)
#     if num == RSDO !=-1:
#         maybe.append(sea_sw)
#     else:
#         if num == od !=-1:
#             maybe.append(sea_dms)
            
#         else:
#             old_vari=np.nanmean(datas[num][:,0:91,:,:],axis=(1))
#             vari=np.empty((13, 65, 144,192)); vari[:]=np.nan
#           #  data=np.nanmean(datas[num],axis=(3))
#             for t in range(0,datas[num].shape[0]):
#                 for i in range(0,vari.shape[1]):
#                     for la in range(0,datas[num].shape[2]):
#                         for lo in range(0,datas[num].shape[3]):
#                         #for x in range(0,datas[num].shape[3]):
#                               #  print('hi')
#                             if num == w !=-1:
#                                 vari[t,i,la,lo]=(datas[num][t,i+100,la,lo]-old_vari[t,la,lo])
#                             else:
#                                 vari[t,i,la,lo]=((datas[num][t,i+100,la,lo]/old_vari[t,la,lo])-1)*100
#             print(np.shape(vari))
#             print(np.nanmean(vari))
#             see=np.nanmean(vari,axis=(0))
#             maybe.append(see)
#             maybe=np.array(maybe)
    
# gc.collect()

# file = '/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/Difference/Rel_Diff/'
# np.save('{}SON_{}.npy'.format(file,name),maybe[0])
# maybe=[]

# for num,datas in zip(numbers,large_data):
#     print(num)
#     if num == RSDO !=-1:
#         maybe.append(sea_sw)
#     else:
#         if num == od !=-1:
#             maybe.append(sea_dms)
            
#         else:
#             old_vari=np.nanmean(datas[num][:,90:111,:,:],axis=(1))
#             vari=np.empty((13, 65, 144,192)); vari[:]=np.nan
#           #  data=np.nanmean(datas[num],axis=(3))
#             for t in range(0,datas[num].shape[0]):
#                 for i in range(0,vari.shape[1]):
#                     for la in range(0,datas[num].shape[2]):
#                         for lo in range(0,datas[num].shape[3]):
#                         #for x in range(0,datas[num].shape[3]):
#                               #  print('hi')
#                             if num == w !=-1:
#                                 vari[t,i,la,lo]=(datas[num][t,i+100,la,lo]-old_vari[t,la,lo])
#                             else:
#                                 vari[t,i,la,lo]=((datas[num][t,i+100,la,lo]/old_vari[t,la,lo])-1)*100
#             print(np.shape(vari))
#             print(np.nanmean(vari))
#             see=np.nanmean(vari,axis=(0))
#             maybe.append(see)
#             maybe=np.array(maybe)
    
# gc.collect()

# file = '/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/Difference/1940_1960_mean/Rel_Diff/SON/'
# np.save('{}{}.npy'.format(file,name),maybe[0])


# h2so4=np.load('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/DJF/Ensamble/Sulfate_MMR_Global.npy')
# ssa=np.load('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/DJF/Ensamble/SSA_MMR_Global.npy')
# vari=np.empty((13, 165, 144,192)); vari[:]=np.nan
# for e in range (0,13):
#     for i in range(0,165):
#         for lo in range(0,144):
#             for la in range(0,192):
#                 vari[e,i,lo,la]=ssa[e,i,lo,la]/h2so4[e,i,lo,la]
# print(np.nanmean(vari))
# np.save('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/SON/Ensamble/SSA_H2SO4.npy',vari)

# h2so4=np.load('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/SON/Ensamble/Sulfate_MMR_Global.npy')
# ssa=np.load('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/SON/Ensamble/SSA_MMR_Global.npy')
# vari=np.empty((13, 165, 144,192)); vari[:]=np.nan
# for e in range (0,13):
#     for i in range(0,165):
#         for lo in range(0,144):
#             for la in range(0,192):
#                 vari[e,i,lo,la]=ssa[e,i,lo,la]/h2so4[e,i,lo,la]
# print(np.nanmean(vari))
# np.save('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/SON/Ensamble/SSA_H2SO4.npy',vari)

# h2so4=np.load('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/MAM/Ensamble/Sulfate_MMR_Global.npy')
# ssa=np.load('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/MAM/Ensamble/SSA_MMR_Global.npy')
# vari=np.empty((13, 165, 144,192)); vari[:]=np.nan
# for e in range (0,13):
#     for i in range(0,165):
#         for lo in range(0,144):
#             for la in range(0,192):
#                 vari[e,i,lo,la]=ssa[e,i,lo,la]/h2so4[e,i,lo,la]
# print(np.nanmean(vari))
# np.save('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/MAM/Ensamble/SSA_H2SO4.npy',vari)    

# h2so4=np.load('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/JJA/Ensamble/Sulfate_MMR_Global.npy')
# ssa=np.load('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/JJA/Ensamble/SSA_MMR_Global.npy')
# vari=np.empty((13, 165, 144,192)); vari[:]=np.nan
# for e in range (0,13):
#     for i in range(0,165):
#         for lo in range(0,144):
#             for la in range(0,192):
#                 vari[e,i,lo,la]=ssa[e,i,lo,la]/h2so4[e,i,lo,la]
# print(np.nanmean(vari))
# np.save('/home/ybh10/CMIP6/UKESM1/Historic/Numpy_Array/Global/JJA/Ensamble/SSA_H2SO4.npy',vari)    


# In[ ]:





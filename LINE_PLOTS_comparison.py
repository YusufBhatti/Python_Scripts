# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 15:59:49 2020

@author: ybh10
"""

import os
os.chdir("../Scripts")
from my_functions import *
os.chdir("../Aerosol_AOD")
import glob
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
from netCDF4 import Dataset
def topsy(data):
    top=[]
    top_RMSE=[]
    bot=[]
    bot_RMSE=[]
    for i in (1,3,5):
        top.append(data[i])
        # top_RMSE.append(RMSE[i])
    for i in (2,4,6):
        bot.append(data[i])
        # bot_RMSE.append(RMSE[i])
    Diff=np.concatenate((top,bot))
    # RMSE_Tot=np.concatenate((top_RMSE,bot_RMSE))
    return Diff
Chems=['CHEM1','CHEM2','CHEM3','CHEM1-SSF','CHEM2-SSF','CHEM3-SSF']
Dates=['January', 'February','March','April','May','June','July','August','September','October','November','December']

inpath=sorted(glob.glob('//macdiarmid/PHYS381/Observations/MODIS_Aqua_c6.1_monthly_AOD/*.nc'))
inpath=inpath[0:60]
vn2 = xr.open_dataset(inpath[0])
lat_o = vn2.lat

MODIS=np.load('AOD_MODIS_DATA.npy')
MODIS=np.concatenate((MODIS[6:12],MODIS[0:6]))
# MODIS2=MODIS[6:0]
# MODIS=np.concatenate((MODIS1,MODIS2))
global_aod=np.load('CHEM_vs_Modis_Difference.npy') # (7,12,180,360)
global_aod=np.concatenate((global_aod[:,6:12],global_aod[:,0:6]),axis=1)
# global_aod1=global_aod[:,6:]
# global_aod2=global_aod[:,:6]
# global_aod=np.concatenate((global_aod1,global_aod2),axis=1)
#### (['REF','CHEM1','CHEM1-SSF','CHEM2','CHEM2-SSF','CHEM3','CHEM3-SSF'])
MODEL_DATA=np.load('Chem_Model_Data.npy')
MODEL_DATA=np.concatenate((MODEL_DATA[:,6:12],MODEL_DATA[:,0:6]),axis=1)
# MODEL_DATA1=MODEL_DATA[:,6:]
# MODEL_DATA2=MODEL_DATA[:,:6]
# MODEL_DATA=np.concatenate((MODEL_DATA1,MODEL_DATA2),axis=1)
global_mean=areaweight(global_aod,lat_o)
lat_bnds = [-60, -40]

lat_inds = np.where((lat_o > lat_bnds[0]) & (lat_o < lat_bnds[1]))
lat_inds = np.squeeze(lat_inds)
AOD_SO=global_aod[:,:,lat_inds[0]:lat_inds[19]+1,:]
lat_o_SO=vn2.lat[lat_inds[0]:lat_inds[19]+1]
MODIS_SO=MODIS[:,lat_inds[0]:lat_inds[19]+1,:]
MODEL_SO=MODEL_DATA[:,:,lat_inds[0]:lat_inds[19]+1,:]
MODEL_SO=np.nanmean(MODEL_DATA,axis=1)
MODIS_SO=np.nanmean(MODIS,axis=0)
AOD_SO_Mean=areaweight(AOD_SO,lat_o_SO) # Weighting the lon and lat, instead of nan meaning - accurate.
AOD_SO_Diff=topsy(AOD_SO_Mean)
# AOD_SO_Diff=AOD_SO_Diff[:,6:]+AOD_SO_Diff[:,:6]

AOD_global_mean=(global_mean)
AOD_global_meany=np.nanmean(global_mean,axis=1)
AOD_global=topsy(AOD_global_meany)
AOD_global_Diff=topsy(AOD_global_mean)
Dates=['January', 'February','March','April','May','June','July','August','September','October','November','December']
Dates=Dates[6:]+Dates[:6]
SO_Mean=np.nanmean(AOD_SO_Diff,axis=1)
RMSEy=[]
from math import sqrt
from sklearn.metrics import mean_squared_error
for i in range(0,7):
    data1= MODIS_SO; data2 = MODEL_SO[i]
    mask1 = np.isnan(data1)
    rms_gong = sqrt(mean_squared_error(data1[~mask1],data2[~mask1]))
    RMSEy.append(rms_gong)
RMSE=topsy(RMSEy)

fig, axes = plt.subplots(nrows=2, ncols=3, sharex=True, figsize=(20,15))
fig.subplots_adjust(hspace=0.2, wspace=0.1)
for ax1,DD,GD,ch,meany,rmse in zip(axes.flat,AOD_SO_Diff,AOD_global_Diff,Chems,AOD_global,RMSE):
    ax1.grid(linestyle='--',alpha=0.5)
    clevs=np.arange(-0.08,0.15+0.05,0.05)    # clevs1=np.arange(-100,100,2)
    ax1.set_ylim([-0.08,0.15])
    ax1.plot(Dates,AOD_SO_Mean[0],color='red',linestyle=':',label='REF Southern Ocean')
    ax1.plot(Dates,AOD_global_mean[0],linestyle=':',color='blue',label='REF Global ')

    cs = ax1.plot(Dates,DD,color='red',label='Southern Ocean')
    ax1.plot(Dates,GD,color='blue',label='Global')
    x=np.arange(0,12,1)
    ax1.set_xticks(x)
    ax1.set_xticklabels(Dates,fontsize=15,rotation=45)
    ax1.set_title(ch+"\n Mean and RMSE = {} and {}".format(round(meany,3),round(rmse,3)),fontsize=20)
    plt.suptitle('Global Mean Difference and RMSE from MODIS',y=.955,fontsize=25,fontweight='bold')
    if ch == 'CHEM1-SSF':
        ax1.legend(fontsize=20)
        ax1.set_ylabel('AOD',fontsize=20)
    if ch == 'CHEM1':
        ax1.legend(fontsize=20)
        ax1.set_ylabel('AOD',fontsize=20)

# plt.show()
plt.savefig("Chem_Line_Plot_Lauras.png",dpi=600,bbox_inches = 'tight')

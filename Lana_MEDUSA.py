# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 14:45:13 2020

@author: ybh10
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 14:36:09 2020

@author: ybh10
"""

import os
os.chdir("../Scripts")
from my_functions import *
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
from netCDF4 import Dataset
Dates=['January', 'February','March','April','May','June','July','August','September','October','November','December']
Date=['J','A','S','O','N','D','J','F','M','A','M','J']
os.chdir("../DMS_Emissions/Lana/2003")


os.chdir("../DMS_Emissions/MEDUSA/2003")

model=['00011','00039','00038','00041','00042','00043','00044']
Chemical=['REF','CHEM1','CHEM1-SSF','CHEM2','CHEM2-SSF','CHEM3','CHEM3-SSF']
inpath=sorted(glob.glob('//macdiarmid/PHYS381/Observations/MODIS_Aqua_c6.1_monthly_AOD/*.nc'))
inpath=inpath[0:60]
vn2 = xr.open_dataset(inpath[0])
lat_o = vn2.lat
lon_o = vn2.lon
filename='//macdiarmid/PHYS381/Aerosols/run_00011/00011_atmosphere_optical_thickness_due_to_insoluble_aitken_mode_sulphate_aerosol_climatology.nc'
vn1 = xr.open_dataset(filename)
lat_m = vn1.latitude
lon_m = vn1.longitude
time_m = vn1.time
lat_ob_SO = vn2.lat
lat_mo_SO = vn1.latitude
AOD_Obv=[]
AOD_Obby=[]
count=-1
d=np.arange(0,60,12)
upper= np.arange(5,61,5)
AOD_O=[]
AOD_Obv_Lon=[]
AOD_Obby_Lon=[]
AOD_O_Lon=[]
for pathway in (inpath): ##OBSERVATIONAL DATA
    lat_bnds = [-60, -40]

    lat_inds = np.where((lat_o > lat_bnds[0]) & (lat_o < lat_bnds[1]))
    lat_inds = np.squeeze(lat_inds)
    Dataset_o=xr.open_dataset(pathway) 
    Aerosol_Optical_O=Dataset_o.MYD08_M3_6_1_AOD_550_Dark_Target_Deep_Blue_Combined_Mean_Mean[:]
    aerosol_o=np.mean(Aerosol_Optical_O,axis=(1))  # Mean in the longitude
    AOD_Obv.append(aerosol_o.data)
    AOD_Obv_Lon.append(Aerosol_Optical_O.data)
    if pathway.find ('20071231.MONTH_12') !=-1: # Gets the last file in the loop, 
        for x in range(0,12): # 60 files produced, and sorts them to make every 5 file the
            count =count+1 # Same 5 months. E.G. It WAS J,A,S,O. But now its JJJJJ,AAAAA,SSSSS...
            print(count)
            for i in (d): # The d arange just sorts it every 5 files
                print(i+count)
                print(x)
                a=AOD_Obv[i+count]
                a_lon=AOD_Obv_Lon[i+count]
                AOD_Obby.append(a)
                AOD_Obby_Lon.append(a_lon)
                if i+count == 59 !=-1: ### The 59 is the last file for it to iterate though, resulting
                    for u in (upper): ## in the loop to go through this bit when it's done above.
                        a=AOD_Obby[u-5:u]  
                        a_lon=AOD_Obby_Lon[u-5:u]  
                        b=np.nanmean(a,axis=0) ## This loop gets the mean of all the same months - creating
                        b_lon=np.nanmean(a_lon,axis=0) ## This loop gets the mean of all the same months - creating
                        AOD_O.append(b) ## a CLIMATOLOGY of 5 years, for each month, with 20 lat data (12,20)
                        AOD_O_Lon.append(b_lon) ## a CLIMATOLOGY

lat_o_SO=Dataset_o.lat[lat_inds[0]:lat_inds[19]+1]
AOD_O_Lon=np.array(AOD_O_Lon)
mask_test = np.isnan(AOD_O_Lon)
AOD_Diff=[]
# lat_o_SO=Dataset_o.lat[lat_inds[0]:lat_inds[19]+1]
AOD_M_Lon=[]
AOD_M=[]
AOD_Test_WO_Mod=[]
for run in (model): #MODEL dATA
    Dataset=xr.open_dataset(run+'_AOD.nc')
    lat_inds = np.where((lat_m > lat_bnds[0]) & (lat_m < lat_bnds[1]))
    lat_inds = np.squeeze(lat_inds)
    lat_m_SO=Dataset.latitude[lat_inds[0]:lat_inds[15]+1] ## Create the latitude variable for the model run
    Aerosol_Optical=Dataset.Aerosol_Optical_Depth[:,2]
    Aerosol_oppy=check_lat_lon(lat_m,lon_m,Aerosol_Optical.data)
    Aerosol_model=Aerosol_oppy[2]
    lats_model=Aerosol_oppy[0]
    lons_model=Aerosol_oppy[1]
    AOD_M_Modified=my_interpolate(lats_model,lons_model,Aerosol_model,lat_o,lon_o)
    # Aerosol_Optical=AOD_M_Modified[:,lat_inds[0]:lat_inds[15]+1,:]
    # AOD_M_Modifiedy=np.flip(AOD_M_Modified,axis=(1))
    AOD_M_Modified[mask_test] = np.nan
    AOD_M_Lon.append(AOD_M_Modified)

for i in range(0,7): ### Gets the difference for all SEVEN model runs compared to MODIS data.
    aod_difference=AOD_M_Lon[i]-AOD_O_Lon
    AOD_Diff.append(aod_difference)
AOD_M_Lon=np.array(AOD_M_Lon)
AOD_Diff=np.array(AOD_Diff)
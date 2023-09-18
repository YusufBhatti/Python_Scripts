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
#from my_functions import *
os.chdir("/nesi/project/niwa02757/ybh10/CMIP6/UKESM1/Historic/")
def calc_pdf(data, left, right, nbins):
    '''
    This is a function to calculate pdf
data: data to calculate PDF must be 1D
left: left bound of data
right: right bound of data
    nbins: number of bins
    '''
    data = data[~np.isnan(data)]
    right = right-0.1
    
    bins = np.linspace(left, right, nbins)
    
    frequency = np.zeros(len(bins))
    
    for i in range(0,nbins):
        if np.isnan(np.nanmean(data)) == False:
            if i < nbins-1:
                mask = (bins[i]<= data) & (data < bins[i+1])
                frequency[i] = len(data[mask])/len(data)
            else:
                mask = data >= bins[i]
                frequency[i] = len(data[mask])/len(data)
        else:
            frequency[i] = np.nan
    return(bins,frequency)

def calc_cdf(perd):
    cpd = np.zeros(len(perd))
    for i in range(0,len(perd)):
        cpd[i] = perd[i]+ cpd[i-1]
    return(cpd)


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]
# bins, freq = calc_pdf(wspd_cha, left=0, right=15, nbins=20)
# cdf = calc_cdf(freq)
def check_lev_lat_lon(lat,lon,array):
    #like check_lat_lon but ignores level dimension in 2nd position

    #ensure lat runs N-S and lon runs 180W-180E
    if lat[0]<1:
        lat = lat[::-1]
        array = array[:,:,::-1,:]
    if lon[-1]>180:
        array = t_lon_array(lon,array)
        lon = t_lon(lon)
def check_lat_lon_model(lat,lon,array):
    #ensure lat runs N-S and lon runs 180W-180E
#     if lat[0]<1:
#         lat = lat[::-1]
#         array = array[:,::-1,:]
    if lon[-1]<180:
        array = t_lon_array(lon,array)
        lon = t_lon_360(lon)
        
    return lat, lon, array
def check_lev_lat_lon_model(lat,lon,array):
    #like check_lat_lon but ignores level dimension in 2nd position

    #ensure lat runs N-S and lon runs 180W-180E
#     if lat[0]<1:
#         lat = lat[::-1]
#         array = array[:,:,::-1,:]
    if lon[-1]<180:
        array = t_lon_array(lon,array)
        lon = t_lon_360(lon)        
    return lat, lon, array
###############################################################################
def t_lon(input_lon):
    #Reorders longitude from 0:360 to -180:180
    h = np.int(np.size(input_lon)/2)
    t_lon = np.mod(input_lon+180,360)-180
    output_lon = np.concatenate((t_lon[h:],t_lon[:h]),axis=0)
    
    return output_lon

def t_lon_360(input_lon):
    #Reorders longitude from -180:180 to 0:360 
    h = np.int(np.size(input_lon)/2)
    t_lon = np.mod(input_lon+360,360)
    output_lon = np.concatenate((t_lon[h:],t_lon[:h]),axis=0)
    
    return output_lon
###############################################################################    
def t_lon_array(input_lon,input_array):
    #Reorders an array so longitudes from 0:360 go to -180:180 (assumes that lon is the last dimension)
    h = np.int(np.size(input_lon)/2)
    
    if input_array.ndim==2:
        output_array = np.concatenate((input_array[:,h:],input_array[:,:h]),axis=1)
    if input_array.ndim==3:
        output_array = np.concatenate((input_array[:,:,h:],input_array[:,:,:h]),axis=2)
    if input_array.ndim==4:
        output_array = np.concatenate((input_array[:,:,:,h:],input_array[:,:,:,:h]),axis=3)
    
    return output_array

def t_lon_array_360(input_lon,input_array):
    #Reorders an array so longitudes from -180:180 go to 0:360 (assumes that lon is the last dimension)
    h = np.int(np.size(input_lon))
    
    if input_array.ndim==2:
        output_array = np.concatenate((input_array[:,h:],input_array[:,:h]),axis=1)
    if input_array.ndim==3:
        output_array = np.concatenate((input_array[:,:,h:],input_array[:,:,:h]),axis=2)
    if input_array.ndim==4:
        output_array = np.concatenate((input_array[:,:,:,h:],input_array[:,:,:,:h]),axis=3)
    
    return output_array        



### CHLOROPHYLL MONTH DATA MODIS

########### DAILY MODIS DATA  ################
# year_count=0
# modis=[]
# month_chl=np.empty((12,4,4320,8640)); month_chl[:]=np.nan


# year=['2020']
# month=['01','02','03','04','05','06','07','08','09','10','11','12']
# day=['01','10','20','30']
# mon_count=-1
# day_count=0
# for m in (month):
#     mon_count=mon_count+1
#     day_count=0
#     for d in (day):
#         modis_data=xr.open_dataset('/home/ybh10/Observational_Data/MODIS/CHL/Daily/2020/g4.timeAvgMap.MODISA_L3m_CHL_8d_4km_2018_chlor_a.2020{}{}-2020{}{}.180W_90S_180E_90N.nc'.format(m,d,m,d)).MODISA_L3m_CHL_8d_4km_2018_chlor_a
#         month_chl[mon_count,day_count,:,:]=modis_data[:,:]
#         modis.append(modis_data[:,:])
#          #   print('day_count={}'.format(day_count))
#         day_count=day_count+1
# print('mean of chl 2020 dec {}'.format(np.nanmean(month_chl[11])))
# print('mean of chl 2020 dec, last day {}'.format(np.nanmean(month_chl[11,3])))


#month_chl=np.load('/home/ybh10/Objective_2/Chlorophyll_Climatology/2020/Daily_2020_MODIS_GRID.npy')
month_chl=np.load('/home/ybh10/Objective_2/Chlorophyll_Climatology/Monthly_CHL_MODIS_2003_2020.npy')
month_chl=np.nanmedian(month_chl,axis=0)
print(np.shape(month_chl))

modis_data_2=xr.open_dataset('/home/ybh10/Observational_Data/MODIS/CHL/Daily/2009_2021_Jan/g4.timeAvgMap.MODISA_L3m_CHL_8d_4km_2018_chlor_a.20180101-20180131.180W_90S_180E_90N.nc')
olat=modis_data_2.lat.data
olon=modis_data_2.lon.data
            


### CALCULATES EACH GRID CELL ########################
modis_check=check_lat_lon_model(olat,olon,month_chl) # for DAILY DATA
#modis_check=check_lev_lat_lon_model(olat,olon,month_chl) # for DAILY DATA
modis_new_lat=modis_check[0]
modis_new_lon=modis_check[1]
modis_new_data=modis_check[2]
lat_num=1080
lon_num=2158
interval=4
first_dim=18
sec_dim=12

latitude=np.empty((lat_num)); latitude[:]=np.nan
longitude=np.empty((lon_num)); longitude[:]=np.nan
mean_chl=np.empty((sec_dim,lat_num,lon_num))
count=0
modis_bins=[]
mean_chl=np.empty((sec_dim,lat_num,lon_num)); mean_chl[:]=np.nan
std_chl=np.empty((sec_dim,lat_num,lon_num)); std_chl[:]=np.nan
median_chl=np.empty((sec_dim,lat_num,lon_num)); median_chl[:]=np.nan
IQR=np.empty((sec_dim,lat_num,lon_num)); IQR[:]=np.nan
laty=np.arange(0,4320,interval)
lony=np.arange(0,8629,interval)
lo_co=-1
#for m in range(0,first_dim):
for d in range(0,sec_dim):
    la_co=-1
    for la in (laty):
        la_co=la_co+1
        for lo in (lony):
            lo_co=lo_co+1
            data=(modis_new_data[d,la:la+interval,lo:lo+interval])
            latitude[la_co]=np.nanmean(modis_new_lat[la:la+interval])
            longitude[lo_co]=np.nanmean(modis_new_lon[lo:lo+interval])
            median_chl[d,la_co,lo_co]=np.nanmedian(data)
            std_chl[d,la_co,lo_co]=np.nanstd(data)
            mean_chl[d,la_co,lo_co]=np.nanmean(data)
            Q1=np.nanquantile(data,.25); Q3=np.nanquantile(data,.75)
            IQR[d,la_co,lo_co]=Q3-Q1

            if lo_co == lon_num-1:
                lo_co=-1
              #  print(la_co)

            
File_Path='/home/ybh10/Objective_2/Chlorophyll_Climatology/2020'
Data_Type='ANNUAL_CLIMATOLOGY_1080x2158_Median'

np.save('/home/ybh10/Objective_2/Chlorophyll_Climatology/ANNUAL_CLIMATOLOGY_MODIS_GRID.npy',month_chl)   
np.save('{}/Median_{}.npy'.format(File_Path,Data_Type),median_chl)   
np.save('{}/Mean_{}.npy'.format(File_Path,Data_Type),mean_chl) 
np.save('{}/STD_{}.npy'.format(File_Path,Data_Type),std_chl) 
np.save('{}/IQR_{}.npy'.format(File_Path,Data_Type),IQR)
np.save('{}/{}_MODIS_GRID.npy'.format(File_Path,Data_Type),month_chl)
np.save('{}/{}_ORCA1_Latitude.npy'.format(File_Path,Data_Type),latitude)
np.save('{}/{}_ORCA1_Longitude.npy'.format(File_Path,Data_Type),longitude)
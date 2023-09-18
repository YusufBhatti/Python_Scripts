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
#os.chdir("/nesi/project/niwa02757/ybh10/CMIP6/UKESM1/Historic/")
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


chl_medusa=xr.open_dataset('/home/ybh10/Objective_2/MEDUSA/coupled/Processed/daily_Ocean_Surface_Chlor_coupled_jan.nc')
medusa_chl=chl_medusa.m01s00i096
model_mean_chl=np.nanmean(medusa_chl,axis=0)*1e6
lat=medusa_chl.latitude.data
lon=medusa_chl.longitude.data
#obs_check=check_lat_lon(lat,lon,medusa_chlviridis.data)

# ##data=interpolate(modis_data_2,olat,olon,lat,lon,model_mean_chl)
# a=my_interpolate(obs_check[0],obs_check[1],obs_check[2],olat,olon)
# model_int=np.nanmean(a,axis=0)[:1200]*1e6
#bin_number=90
bin_number='10000_bins_90'
Model='MODEL_Grid/Max_Concentration_{}/MODEL_GRID_'.format(bin_number)

cdf_modis_g=np.load('/home/ybh10/Objective_2/Numpy_Array/CDF/{}MODIS_CDF_SO_MAX_{}.npy'.format(Model,bin_number))
mean_chl=np.load('/home/ybh10/Objective_2/Numpy_Array/CDF/{}interpolated_MODIS_CHL_Spatially_Mean_Data_{}.npy'.format(Model,bin_number))
std_chl=np.load('/home/ybh10/Objective_2/Numpy_Array/CDF/{}interpolated_MODIS_CHL_Spatially_STD_Data_{}.npy'.format(Model,bin_number))
median_chl=np.load('/home/ybh10/Objective_2/Numpy_Array/CDF/{}interpolated_MODIS_CHL_Spatially_Median_Data_{}.npy'.format(Model,bin_number))
modis_freq=np.load('/home/ybh10/Objective_2/Numpy_Array/CDF/{}MODIS_modis_freq_SO_MAX_{}.npy'.format(Model,bin_number))
latitude=np.load('/home/ybh10/Objective_2/Numpy_Array/CDF/{}MODIS_CDF_SO_Latitude.npy'.format(Model))
longitude=np.load('/home/ybh10/Objective_2/Numpy_Array/CDF/{}MODIS_modis_freq_SO_Longitude.npy'.format(Model))
obins=np.load('/home/ybh10/Objective_2/Numpy_Array/CDF/{}MODIS_BINS_{}.npy'.format(Model,bin_number))


import random

count=0
lat=[]
lon=[]
day=[]
data_model=np.empty((30,144,192)); data_model[:]=np.nan
data_modis=np.empty((30,144,192)); data_modis[:]=np.nan
for t in range(0,30):
    for la in range(0,144):
        for lo in range(0,192):
            n = round(random.uniform(0, 1.00000000), 8)# outputs a RANDOM number between 0.01 - 1 
            occ=find_nearest(cdf_modis_g[:,la,lo], n) # This points to the nearest value in the cdf to n MODIS
            if np.isnan(occ) == True:
                data_modis[t,la,lo]=np.nan
            else:
                bins=np.where(cdf_modis_g[:,la,lo]==occ)[0][0] # The CHL bin which contains the that cdf occ value MODIS
                CHL_CON=obins[bins]
                data_modis[t,la,lo]=CHL_CON

np.save('/home/ybh10/Objective_2/Numpy_Array/CDF/MODEL_Grid/Max_Concentration_{}/MODEL_GRID_RNG_STOCHASTIC_Data_{}'.format(bin_number,bin_number),data_modis)  
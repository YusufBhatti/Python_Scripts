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
    from my_functions import t_lon
    from my_functions import t_lon_array
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
    from my_functions import t_lon_360
    from my_functions import t_lon_array
    #ensure lat runs N-S and lon runs 180W-180E
    if lat[0]<1:
        lat = lat[::-1]
        array = array[:,:,::-1,:]
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
modis_data_2=xr.open_dataset('/home/ybh10/Objective_2/MODIS_CHL/g4.timeAvgMap.MODISA_L3m_CHL_8d_4km_2018_chlor_a.20180101-20180131.180W_90S_180E_90N.nc').MODISA_L3m_CHL_8d_4km_2018_chlor_a
olat=modis_data_2.lat.data
olon=modis_data_2.lon.data
year = ['2009','2010','2011','2012','2013','2014','2015','2016','2017','2018']
month='01'
day=[1,10,20,30]
month_chl=np.empty((10,4,4320,8640)); month_chl[:]=np.nan
year_count=0
day_count=0
for y in range(2009,2019):
    for d in (day):
        modis_data=xr.open_dataset('/home/ybh10/Objective_2/MODIS_CHL/g4.timeAvgMap.MODISA_L3m_CHL_8d_4km_2018_chlor_a.{:02d}01{:02d}-{:02d}01{:02d}.180W_90S_180E_90N.nc'.format(y,d,y,d)).MODISA_L3m_CHL_8d_4km_2018_chlor_a
    #     print(np.nanmean(modis_data))
         #print('/home/ybh10/Objective_2/MODIS_CHL/g4.timeAvgMap.MODISA_L3m_CHL_8d_4km_2018_chlor_a.{}01{}-{}01{}.180W_90S_180E_90N.nc'.format(y,d,y,d))
        month_chl[year_count,day_count,:,:]=modis_data[:,:]  
        day_count=day_count+1
     #   print('day_count={}'.format(day_count))
    day_count=0
    year_count=year_count+1
    #print('year_count={}'.format(year_count))


# model_mean_chl=np.nanmean(medusa_chl,axis=0)*1e6
# lat=medusa_chl.latitude.data
# lon=medusa_chl.longitude.data

### CALCULATES CDF FOR EACH 1.25 (120km) X 1.875 (180km) GRID CELL ########################
modis_check=check_lat_lon_model(olat,olon,month_chl)
modis_new_lat=modis_check[0]
modis_new_lon=modis_check[1]
modis_new_data=modis_check[2]
number_of_bins=10000
max_concentration=50

cdf_modis_grid=np.empty((number_of_bins,173,346)); cdf_modis_grid[:]=np.nan
cdf_modis_grid=np.empty((number_of_bins,173,346)); cdf_modis_grid[:]=np.nan
latitude=np.empty((173)); latitude[:]=np.nan
longitude=np.empty((346)); longitude[:]=np.nan
mean_chl=np.empty((173,346))
count=0
modis_bins=[]
modis_freq=np.empty((59858,number_of_bins)); modis_freq[:]=np.nan
mean_chl=np.empty((173,346)); mean_chl[:]=np.nan
std_chl=np.empty((173,346)); std_chl[:]=np.nan
median_chl=np.empty((173,346)); median_chl[:]=np.nan
IQR=np.empty((173,346)); IQR[:]=np.nan


laty=np.arange(0,4320,25)
lony=np.arange(0,8629,25)
la_co=-1
lo_co=-1
for la in (laty):
    la_co=la_co+1
    for lo in (lony):
        lo_co=lo_co+1
        data=(modis_new_data[:,:,la:la+25,lo:lo+25])
        latitude[la_co]=np.nanmean(modis_new_lat[la:la+25])
        longitude[lo_co]=np.nanmean(modis_new_lon[lo:lo+25])
        median_chl[la_co,lo_co]=np.nanmedian(data)
        std_chl[la_co,lo_co]=np.nanstd(data)
        Q1=np.nanquantile(data,.25); Q3=np.nanquantile(data,.75)
        IQR[la_co,lo_co]=Q3-Q1
        mean_chl[la_co,lo_co]=np.nanmean(data)
        modis_stats=data.flatten()
        if np.isnan(np.nanmean(modis_stats)) == True:
            obins, ofreq = calc_pdf(modis_stats, left=0, right=max_concentration, nbins=number_of_bins)
            cdf_modis_grid[:,la_co,lo_co] = np.nan            
        else:
            obins, ofreq = calc_pdf(modis_stats, left=0, right=max_concentration, nbins=number_of_bins)
            cdf_modis_grid[:,la_co,lo_co] = calc_cdf(ofreq)
     #   modis_bins.append(obins)
        modis_freq[count,:]=ofreq
        count=count+1
        if lo_co == 345:
            lo_co=-1
            print(la_co)
        
        gc.collect()
import random

count=0
lat=[]
lon=[]
day=[]
data_model=np.empty((30,173,346)); data_model[:]=np.nan
data_modis=np.empty((30,173,346)); data_modis[:]=np.nan
for t in range(0,30):
    for la in range(0,173):
        for lo in range(0,346):
            n = round(random.uniform(0, 1.00000000), 8)# outputs a RANDOM number between 0.01 - 1 
            occ=find_nearest(cdf_modis_grid[:,la,lo], n) # This points to the nearest value in the cdf to n MODIS
            if np.isnan(occ) == True:
                data_modis[t,la,lo]=np.nan
            else:
                bins=np.where(cdf_modis_grid[:,la,lo]==occ)[0][0] # The CHL bin which contains the that cdf occ value MODIS
                CHL_CON=obins[bins]
                data_modis[t,la,lo]=CHL_CON
                
print('Stochastic data shape should be 30,144,192, it is: {}'.format(np.shape(data_modis)))
print('Stochastic data mean should be >1, it is: {}'.format(np.nanmean(data_modis)))
print('Stochastic Max CHL Concentration is: {}'.format((max_concentration)))
print('Stochastic number of bins is: {}'.format((number_of_bins)))

bin_number='{}_bins_{}'.format(number_of_bins,max_concentration)
grid='Max_Concentration_{}/'.format(bin_number)
np.save('/home/ybh10/Objective_2/Numpy_Array/CDF/ORCA1_Grid/{}MODEL_GRID_MODIS_BINS_{}'.format(grid,bin_number),obins)  
np.save('/home/ybh10/Objective_2/Numpy_Array/CDF/ORCA1_Grid/{}MODEL_GRID_MODIS_CDF_SO_MAX_{}'.format(grid,bin_number),cdf_modis_grid)  
np.save('/home/ybh10/Objective_2/Numpy_Array/CDF/ORCA1_Grid/{}MODEL_GRID_MODIS_modis_freq_SO_MAX_{}'.format(grid,bin_number),modis_freq)      
np.save('/home/ybh10/Objective_2/Numpy_Array/CDF/ORCA1_Grid/{}MODEL_GRID_MODIS_CDF_SO_Latitude'.format(grid),latitude)  
np.save('/home/ybh10/Objective_2/Numpy_Array/CDF/ORCA1_Grid/{}MODEL_GRID_MODIS_modis_freq_SO_Longitude'.format(grid),longitude)   
np.save('/home/ybh10/Objective_2/Numpy_Array/CDF/ORCA1_Grid/{}MODEL_GRID_interpolated_MODIS_CHL_Spatially_Median_Data_{}'.format(grid,bin_number),median_chl)   
np.save('/home/ybh10/Objective_2/Numpy_Array/CDF/ORCA1_Grid/{}MODEL_GRID_interpolated_MODIS_CHL_Spatially_STD_Data_{}'.format(grid,bin_number),std_chl)  
np.save('/home/ybh10/Objective_2/Numpy_Array/CDF/ORCA1_Grid/{}MODEL_GRID_interpolated_MODIS_CHL_Spatially_IQR_Data_{}'.format(grid,bin_number),IQR)  
np.save('/home/ybh10/Objective_2/Numpy_Array/CDF/ORCA1_Grid/{}MODEL_GRID_interpolated_MODIS_CHL_Spatially_Mean_Data_{}'.format(grid,bin_number),mean_chl)  
np.save('/home/ybh10/Objective_2/Numpy_Array/CDF/ORCA1_Grid/{}MODEL_GRID_RNG_STOCHASTIC_Data_{}'.format(grid,bin_number),data_modis)  

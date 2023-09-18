
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 13:28:58 2019

@author: ee18yb - Yusuf Bhatti
"""

"""
Read in netcdf files and read in AOD

"""
### Import modules ###

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.colors as mcolors
import cartopy.feature as cfeature
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap, cm
from matplotlib.dates import date2num
from datetime import datetime, timedelta
from netCDF4 import num2date, date2num, date2index


### Read in and define netcdf files of monthly 10 Tg SO2 progression ###
### by defining the directory pathway ###

expname='xnsbi'

### Dates in year-month-day ###
mnames=['01014','01028','01111','01125','01209','01223','10106','10120','10203']

### Directory pathway ###
inpath = 'E://Data Files (NetCDF)/'+expname+'/Daily_AOD/'


### Create the subplots (nrows, ncolns) and set the size of plots ###
fig, axes = plt.subplots(nrows=4, ncols=2, sharex=True, figsize=(15.745,15))
fig.subplots_adjust(hspace=-0.033, wspace=0.2)

###################################################################################
######################### Extracting the data for the plots #######################
###################################################################################

### Looping the variables needed to be extracted for the different months ###
for ax,name in zip(axes.flat, mnames):
    infile= inpath+expname+'a_pa200'+name+'sAODtracersMDSJamesB.nc' # Complete the pathway to loop all files # 
    vn=Dataset(infile)
    lon       = vn.variables['longitude'][:]
    lat      = vn.variables['latitude'][:]
    time      = vn.variables['t'][:]
    Ait_Ins0 = vn.variables['AITKEN_MODE__INSOL__STRATO_AOD'][:]
    Ait_Sol0 = vn.variables['AITKEN_MODE__SOLUBLE__STRATO__AOD'][:]
    Coa_Ins0 = vn.variables['COARSE_MODE__INSOL__STRATO_AOD'][:]
    Coa_Sol0 = vn.variables['COARSE_MODE__SOLUBLE__STRATO_AOD'][:]
    Acc_Ins0 = vn.variables['ACCUM_MODE__INSOL__STRATO_AOD'][:]
    Acc_Sol0 = vn.variables['ACCUM_MODE__SOLUBLE__STRATO_AOD'][:]
    
### Add up all components [above] to get an overall aod ###
    aod       = Ait_Ins0+Ait_Sol0+Coa_Ins0+Coa_Sol0+Acc_Ins0+Acc_Sol0
    print(np.shape(aod))
    
### Extract the second dimension to obtain the ###
    aod550=np.squeeze(aod[:,2,:,:])
    print(np.shape(aod550))
   
### Set the time and date of each plots ####
    t_unit="days since 2000-09-01T00:00:00Z"
    #t_cal =  vn1.variables[tname].calendar
    t_cal =  "360_day"

    dates =num2date(time,units=t_unit,calendar=t_cal)
    print(dates)
    ntime=len(time)
    
### Creating map variables ###
    import  matplotlib
    
### set the font size globally to get the ticklabels big ###

    matplotlib.rcParams['font.sans-serif'] = "Arial"
    matplotlib.rcParams['font.family'] = "sans-serif"

    label_size=15
    mpl.rcParams['xtick.labelsize'] =label_size
    mpl.rcParams['ytick.labelsize'] =label_size

    plt.rc('grid', color = 'black')
    plt.rc('grid', alpha = 0.3) # alpha is percentage of transparency
    
    plt.rcParams['xtick.labelsize'] = 10
    plt.rcParams['ytick.labelsize'] = 10
    
    
###################################################################################
###################### Ploting models for the Months ##############################
###################################################################################


#### Setting up a [Mercator] projection with the focus at 180 W and ####
################ coarse resoltuion coastlines ##########################
    map_ax = Basemap(llcrnrlon=0.,llcrnrlat=-80,urcrnrlon=358.,urcrnrlat=85,
            projection='cyl',resolution ='c',ax=ax)
    
### Add coastlines and countries ###
    map_ax.drawcoastlines()
    map_ax.drawcountries()
    
### organise the longitude and latitude using meshgrid ###
    lons, lats = np.meshgrid(lon, lat)
### Define x and y as Longitude and Latitude ###
    x, y = map_ax(lons, lats) 
    
### Plotting Latitude and Longitude grid every 30 and 60 degrees, respectively ###
    meridians = np.arange(0.,360,60.)    
    parallels = np.arange(-90.,90,30.)
    map_ax.drawmeridians(meridians,labels=[0,0,0,1],fontsize=12)  
    map_ax.drawparallels(parallels,labels=[1,0,0,0],fontsize=12)


### Set the axis ###
    clevs=np.arange(0+0.025,0.7+0.025,0.025)
### Change values outside contour limits to min and max contours limits
    aod550[aod550 > 0.7] = 0.7
    #aod550[aod550 < 0.025] = -10



###################################################################################
################################## Construct the plot #############################
###################################################################################

    cs = ax.contourf(lon,lat,aod550,clevs,cmap='jet',linewidths=1.)
    
### Add colourbar ###    
    cbar = map_ax.colorbar(cs,location='right',pad="5%")
    cbar.set_label('stratospheric AOD',fontsize=12)
    sdate=str(dates[0])
    
### Add title with the looped dates ###
    plt.title("10 Tg SO$_2$ @ 21-23 km height"+" ("+sdate[0:10]+")",fontsize=14)
    
### Save and display the plot as a .png and save in high quality (dpi=600)###
    plt.savefig("E://Model_Plots/xnsbi_20/Daily_14days_8_" +expname+".png",dpi=600,bbox_inches = 'tight')
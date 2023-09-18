        # -*- coding: utf-8 -*-
"""
Created on Wed 12/02/2020 10:34:50

@author:  Yusuf Bhatti
"""

"""
Read in netcdf files and read in AOD

"""
### Import modules ###


import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
#import cartopy.crs as ccrs
import matplotlib.colors as mcolors
import cartopy.feature as cfeature
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap, cm
from matplotlib.dates import date2num 
from datetime import datetime, timedelta
from netCDF4 import num2date, date2num, date2index
import datetime as dt
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import matplotlib.gridspec as gridspec
from matplotlib import cm
import matplotlib.colors as colors

#'00060','00062','00056','00058','00063','00064','00065','00066'
model=['00056','00058']

months = ['jan']
years = [97,98,99]
date=[]
for y in range(0, len(years)):
    for m in range(0, len(months)):
        date.append(str(years[y]) + months[m])

def areaweight(x,latitude):
    #Area weights variable x by the cosine of the latitude.
    cos_lat = np.cos(latitude*np.pi/180)
   
    #assume two dimensions for now, LR, 1/3/17
               #add for three dimensions, LR, 8/1/18
    s = x.shape
    n = x.ndim
              
    if n==2 and s[0]==np.size(latitude):
        cost3m = np.empty((s[0],s[1])); cost3m[:,:] = np.nan
        for i in range(0,s[1]):
            cost3m[:,i]= cos_lat
        m = np.nansum(np.nansum(x*cost3m,1),0)/np.nansum(np.nansum(cost3m,1),0) #a lat-lon array        print('a')
                             
    elif n==2 and s[0]!=np.size(latitude):
        cost3m = np.empty((s[0],s[1])); cost3m[:,:] = np.nan
        for i in range(0,s[0]):
            cost3m[i,:]= np.squeeze(cos_lat)
        m = np.nansum(x*cost3m,1)/np.nansum(cost3m,1) #a lev-lat or time-lat array
        print('b')

       
    elif n==3:
        cost3m = np.empty((s[0],s[1],s[2])); cost3m[:,:,:] = np.nan
        for i in range(0,s[0]):
            for j in range(0,s[2]):
                cost3m[i,:,j] = np.squeeze(cos_lat)
        m = np.nansum(np.nansum(x*cost3m,2),1)/np.nansum(np.nansum(cost3m,2),1)        #lat is the middle dimension, e.g. time-lat-lon or lev-lat-lon.     
       
    elif n==4:
        cost3m = np.empty((s[0],s[1],s[2],s[3])); cost3m[:,:,:,:] = np.nan
        for i in range(0,s[0]):
            for j in range(0,s[1]):
                for k in range(0,s[3]):
                    cost3m[i,j,:,k] = np.squeeze(cos_lat)
        m = np.nansum(np.nansum(x*cost3m,3),2)/np.nansum(np.nansum(cost3m,3),2)               #time-lev-lat-lon             
              
    output = m
    return output

filename = '/nesi/nobackup/niwa02757/ybh10/nudged_2018/2018_jan_test/'
dms= filename+'bs798_DMS_2018jan.nc'





lons, lats = np.meshgrid(lon, lat)

import  matplotlib
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"

label_size=5
mpl.rcParams['xtick.labelsize'] =label_size
mpl.rcParams['ytick.labelsize'] =label_size

plt.rc('grid', color = 'black')
plt.rc('grid', alpha = 0.3) # alpha is percentage of transparency
    
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
    
lons, lats = np.meshgrid(lon, lat)
clevs = np.arange(-24,24+4,4)
fig = plt.figure(figsize=(16,9), edgecolor='w')
ax = fig.add_subplot(2,2,1)

ax1 = fig.add_subplot(2,2,2)
ax2 = fig.add_subplot(2,2,3)
ax3 = fig.add_subplot(2,2,4)


fig.subplots_adjust(hspace=0.263, wspace=0.2)

#map_ax = Basemap(projection='spstere', boundinglat=-20., lon_0=90,                            #llcrnrlon=0.,llcrnrlat=-90,urcrnrlon=360.,urcrnrlat=90
#                         resolution ='c',ax=ax)

weightall = areaweight(gaall,lat)
weightcs = areaweight(gacs,lat)
weightcloud = areaweight(cloudy,lat)

weightSSA = areaweight(SSA,lat)

print("weightall =",(weightall))

print("weightcs=",(weightcs))
print("weightcloud=",(weightcloud))

print("weightSSA=",(weightSSA))

fig.suptitle('Gong Atmos',fontsize=20) # or plt.suptitle('Main title')

        ### Plotting Latitude and Longitude grid every 30 and 60 degrees, respectively ###
map_ax = Basemap(llcrnrlon=0.,llcrnrlat=-85,urcrnrlon=360.,urcrnrlat=85,
            projection='cyl',resolution ='c',ax=ax)
map_ax.drawcoastlines()
map_ax.drawcountries() 
x, y = map_ax(lons, lats) 
lonpt, latpt = map_ax(x,y,inverse=True)
meridians = np.arange(0.,360,60.)  

  
parallels = np.arange(-90.,90,30.)
map_ax.drawmeridians(meridians,labels=[0,0,0,1],fontsize=15)  
map_ax.drawparallels(parallels,labels=[1,0,0,0],fontsize=15)
gaall[gaall > 24] = 24
gaall[gaall < -24] = -24

cs= ax.contourf(lon,lat,gaall,clevs,cmap='seismic')
####################################################################################################

map_ax1 = Basemap(llcrnrlon=0.,llcrnrlat=-85,urcrnrlon=360.,urcrnrlat=85,
            projection='cyl',resolution ='c',ax=ax1)
x, y = map_ax1(lons, lats) 
lonpt, latpt = map_ax1(x,y,inverse=True)
map_ax1.drawcoastlines()
map_ax1.drawcountries()
gacs [gacs > 24] = 23
gacs[gacs < -24] = -24

map_ax1.drawmeridians(meridians,labels=[0,0,0,1],fontsize=15)  
map_ax1.drawparallels(parallels,labels=[1,0,0,0],fontsize=15)
ax1.contourf(lon,lat,gacs,clevs,cmap='seismic')

####################################################################################################

map_ax2 = Basemap(llcrnrlon=0.,llcrnrlat=-85,urcrnrlon=360.,urcrnrlat=85,
            projection='cyl',resolution ='c',ax=ax2)
x, y = map_ax2(lons, lats) 
lonpt, latpt = map_ax2(x,y,inverse=True)
map_ax2.drawcoastlines()
map_ax2.drawcountries()

map_ax2.drawmeridians(meridians,labels=[0,0,0,1],fontsize=15)  
map_ax2.drawparallels(parallels,labels=[1,0,0,0],fontsize=15)
clevs1=np.arange(-0.20,0.20+0.01,0.01)
cloudy[cloudy > 0.2] = 0.2

cs1=ax2.contourf(lon,lat,cloudy,clevs1,cmap='seismic')
ax2.title.set_text("Cloud Volume Fraction Diff. mean ="+ str(round((weightcloud),6)))
cbar1 = map_ax2.colorbar(cs1,location='bottom',pad="15%")
cbar1.set_label('Cloud Volume Fraction',fontsize=10)

####################################################################################################

map_ax3 = Basemap(llcrnrlon=0.,llcrnrlat=-85,urcrnrlon=360.,urcrnrlat=85,
            projection='cyl',resolution ='c',ax=ax3)
x, y = map_ax2(lons, lats) 
lonpt, latpt = map_ax2(x,y,inverse=True)
map_ax3.drawcoastlines()
#SSA[SSA>1.5e-07]=1.5e-07
#SSA[SSA>1.5e-07]=1.5e-07

map_ax3.drawcountries()
map_ax3.drawmeridians(meridians,labels=[0,0,0,1],fontsize=15)  
map_ax3.drawparallels(parallels,labels=[1,0,0,0],fontsize=15)
clevs2=np.arange(0.8e-07,2.4e-07+0.1e-7,0.05e-7)
SSA[SSA>2.4e-07]=2.4e-07

cs2=ax3.contourf(lon,lat,SSA,clevs2,cmap='Reds')
map_ax3.fillcontinents()

ax3.title.set_text("SSA Concentration Diff")
cbar2 = map_ax3.colorbar(cs2,location='bottom',pad="15%")
cbar2.set_label('kg kg-1',fontsize=10)


#plt.title("00060 (CS) RF JJA Average: 1995-1999 ",fontsize=10)
ax.title.set_text("RF Diff (all sky). mean ="+ str(round((weightall),2))+ 'W/m$^{2}$')
ax1.title.set_text("RF Diff (clear sky). mean ="+ str(round((weightcs),2))+ 'W/m$^{2}$')

                                            
########## Colourbar ##########
# Set the location and size of the colorbar 
cax = fig.add_axes([0.11, 0.5, 0.8, 0.025])
# Add the colorbar
cbar = fig.colorbar(cs, cax=cax, extend='both',orientation='horizontal')
#cbar1 = plt.colorbar(cs1,  extend='both',orientation='horizontal')

#ax.text(5,-78,'mean = '+ str(round(np.mean(SWasga),2))+ 'W/m$^{2}$',fontsize=6)
#cax = fig.add_axes()
#cbar = map_ax.colorbar(cs,location='bottom')
cbar.set_label('W/m$^2$',fontsize=10)

#plt.show
### Save and display the plot as a .png and save in high quality (dpi=600) ###
plt.savefig("Gong_Jaegle_Plots/Gong_Coupled_weighted_Cloud_vs_RF_vs_SSA.png",dpi=600,bbox_inches = 'tight')

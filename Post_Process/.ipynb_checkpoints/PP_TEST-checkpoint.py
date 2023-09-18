#!/usr/bin/env python
# coding: utf-8

# In[ ]:


"""
Created on Mon Feb 17 15:23:46 2020

@author: ybh10
"""

import iris
import os
import numpy as np

months = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
year = '2003'

#for yr in range(0,np.size(year)):
for i in range(0,np.size(months)):
    #  os.chdir('/home/ybh10/cylc-run/u-bx043/share/data/History_Data/')
     #   filename = 'bx043a.pm2003'+months[i]+'.pp'
        #print vars

    
        #cloud_area_fraction  = vars[346]; iris.save(cloud_area_fraction ,"bx_cloud_area_fraction_2003_{}.nc".format(months[i]))
        #mole_concentration_of_dimethyl_sulfide_in_sea_water = vars[401]; iris.save(mole_concentration_of_dimethyl_sulfide_in_sea_water ,"bx_mole_concentration_of_dimethyl_sulfide_in_sea_water_2003_{}.nc".format(months[i]))   
     #   constraint=iris.AttributeConstraint(STASH='m01s01i298')
      #  cube = iris.load(filename,constraint)
  #      os.chdir('/nesi/nobackup/niwa02757/ybh10/DMS/Lana/2003')        

       # iris.save(cube, '/nesi/nobackup/niwa02757/ybh10/DMS/Lana/2003/bx_CNDC_at_cloud_top_2003_{}.nc'.format(months[i]))
       # constraint=iris.AttributeConstraint(STASH='m01s01i299')
       # cube = iris.load(filename,constraint)
       # iris.save(cube, '/nesi/nobackup/niwa02757/ybh10/DMS/Lana/2003/bx_CNDC_weight_cloud_top_2003_{}.nc'.format(months[i]))
  #      CNDC_at_cloud_top= vars[474] ; iris.save(CNDC_at_cloud_top,"bx_CNDC_at_cloud_top_2003_{}.nc".format(months[i]))
   #     CNDC_weight_cloud_top= vars[474] ; iris.save(CNDC_weight_cloud_top,"bx_CNDC_weight_cloud_top_2003_{}.nc".format(months[i]))

      #  toa_outgoing_LW_flux = vars[474] ; iris.save(toa_outgoing_LW_flux,"bx_toa_outgoing_LW_flux_all_sky_2003_{}.nc".format(months[i]))
        os.chdir('/home/ybh10/cylc-run/u-bw945/share/data/History_Data/')       
        filename = 'bw945a.pm2003'+months[i]+'.pp'
        constraint=iris.AttributeConstraint(STASH='m01s01i298')
        cube = iris.load(filename,constraint)
       # os.chdir('/nesi/nobackup/niwa02757/ybh10/DMS/MEDUSA/2003')       

        iris.save(cube, '/nesi/nobackup/niwa02757/ybh10/DMS/MEDUSA/2003/bw_CNDC_at_cloud_top_2003_{}.nc'.format(months[i]))
        constraint=iris.AttributeConstraint(STASH='m01s01i299')
        cube = iris.load(filename,constraint)
        iris.save(cube, '/nesi/nobackup/niwa02757/ybh10/DMS/MEDUSA/2003/bw_CNDC_weight_cloud_top_2003_{}.nc'.format(months[i]))
       # cloud_area_fraction  = vars[332]; iris.save(cloud_area_fraction ,"bw_cloud_area_fraction_2003_{}.nc".format(months[i]))
       # mass_fraction_of_dimethyl_sulfide_in_air = vars[358]; iris.save(mass_fraction_of_dimethyl_sulfide_in_air,"bw_mass_fraction_of_dimethyl_sulfide_in_air_2003_{}.nc".format(months[i]))         
#        mole_concentration_of_dimethyl_sulfide_in_sea_water = vars[387]; iris.save(mole_concentration_of_dimethyl_sulfide_in_sea_water ,"bw_mole_concentration_of_dimethyl_sulfide_in_sea_water_2003_{}.nc".format(months[i]))   
       # toa_outgoing_SW_flux = vars[464] ; iris.save(toa_outgoing_SW_flux,"bw_toa_outgoing_SW_flux_all_sky_2003_{}.nc".format(months[i]))
      #  toa_outgoing_LW_flux = vars[459] ; iris.save(toa_outgoing_LW_flux,"bw_toa_outgoing_LW_flux_all_sky_2003_{}.nc".format(months[i]))
      #  print(toa_outgoing_LW_flux)


# In[ ]:





# In[ ]:





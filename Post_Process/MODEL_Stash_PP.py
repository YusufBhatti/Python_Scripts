#!/usr/bin/env python
# coding: utf-8

# In[13]:


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 17:27:10 2020

@author: ybh10
"""
from six.moves import cPickle as pickle
import iris
import numpy as np
import os
os.chdir("/home/ybh10/Scripts")
### THINGS NEEDED TO CHANGE WHEN CHANGING SUITE ID:
suite='cb134' # the suite ID BUT NOT INCLUDING THE u-
Daily_save_file='/home/ybh10/DMS_Emissions/Chemistry_Scheme/CHEM2.1/Daily' # Directory to save daily variables
Monthly_save_file='/home/ybh10/DMS_Emissions/Chemistry_Scheme/CHEM2.1/Monthly' # Directory to save monthly variables
## Change the months depending on if the model starts on Jan - Dec, or if just one month... etc.
months = ['1989jan','1989feb','1989mar','1989apr','1989may','1989jun','1989jul','1989aug','1988sep','1988oct','1988nov','1988dec']
#months = ['sep','oct'] # for one month. But must get rid of the [i] in the actual code: filename=(('{}a.p{}{}{}').format(suite,use,year,months[i]))
year='1988' ## Important - this is what year your suite-id states. e.g.  cb134a.pm1988sep so here 1988.
#year=np.arange(2003,2010,1) # for multiple years.
## This is to seperate the daily and monthly variables and directories appropriately. The letter here
## is what the suite-id states. e.g. cb134a.pm1988sep here is 'm' - for monthly. cb134a.pd1988sep - 'd' for daily
daily_file_type='d'
monthly_file_type='m'

#os.chdir("/nesi/nobackup/niwa02757/ybh10/cylc-run/u-"+suite+"/share/data/History_Data/")
section = (0,1,2,3,34,38,50) # This gathers the SECTIONS AND ITEMS for the stash requests - in ORDER.
## Change as appropriate, but MUST correspond to the 'Daily_Variables' section below. E.g. 
## DMS_CONCENTRATION_IN_SEAWATER = 132. Ditem = Daily variables. Mitem = Monthly variables. Some are the 
## same, some are different. It depends on what variable you want to extract for daily/monthly.
Ditem = ['096','132'],['208','209','241','298','299'],['205','206','284','285','300','301','302','303'],['209','210','236'],['071','072','073','074','075'], ['201','202','203','439','485','486','487','488'],['214','216','217','219']
Ditems=np.array(Ditem)
Mitem = ['096','132'],['208','209','241','298','299'],['205','206','284','300','301','302','303'],['209','210','236'],['072','073','074','075','102','104','108','114','967','968'], ['201','202','203','294','295','296','297','298','299','300','439','485','486','487','488'],['214','215','216','217','219']
Mitems=np.array(Mitem)
x=-1
M_msi=[]
## This section compiles the stash requests in order to create a fully formatted stash (e.g. m01s00i096)
for i,s in zip(Mitems,section):
    x=x+1
    for it in (i):
        a=('m01s{:02d}i{}').format(section[x],it)
        M_msi.append(a)
D_msi=[]
x=0
for i,s in zip(Ditems,section):
    x=x+1
    for it in (i):
        a=('m01s{:02d}i{}').format(section[x],it)
        D_msi.append(a)
## This section labels what the stash requests ACTUAL name is. This section is NOT FUNDEMENTAL, but very
## useful if you want to know what each variable extracted is named.If changing one of these names, make sure
## to change the Stash Request (m01s##i###) too! They correspond to each other!

Daily_Variables = [
             'All_Sky_Outgoing_SW_Flux_TOA','Clear_Sky_Upward_SW_Flux_TOA ',
'Cloud_Droplet_Number_Concentration','CDNC_Cloud_Top','Weight_for_CDNC',
             'All_Sky_Outgoing_LW_Flux_TOA','Clear_Sky_Outgoing_LW_Flux_TOA','Sulfate_Optical_Depth','atmosphere_optical_thickness_due_to_dust_ambient_aerosol', 
'atmosphere_optical_thickness_due_to_soluble_aitken_mode','atmosphere_optical_thickness_due_to_soluble_accumulation_mode','atmosphere_optical_thickness_due_to_soluble_coarse_mode', 
'atmosphere_optical_thickness_due_to_insoluble_aitken_mode',
             '10m_Wind_u','10m_wind_v','1.5m_Temp',
             'DMS_mmr','SO2_mmr','H2SO4_mmr','MSA_mmr','DMSO_MMR',
             'H2SO4_to_aitken', 'H2SO4_to_accum','H2SO4_to_coarse','CCN_number_concentration', 
'NUCLEATION_MODE_SOL_H2SO4', 'AITKEN_MODE_SOL_H2SO4','ACCUMULATION_MODE_SOL_H2SO4','COARSE_MODE_SOL_H2SO4',
             'DMS_surface_emissions','so2_high_lev_emissions','so2_natural_emissions', 'ozone_column']

Monthly_Variables = ['Ocean_surface_chlorophyll','Seawater_DMS_Concentration',
             'All_Sky_Outgoing_SW_Flux_TOA','Clear_Sky_Upward_SW_Flux_TOA ',
'Cloud_Droplet_Number_Concentration','CDNC_Cloud_Top','Weight_for_CDNC','All_Sky_Outgoing_LW_Flux_TOA','Clear_Sky_Outgoing_LW_Flux_TOA','Sulfate_Optical_Depth',
'atmosphere_optical_thickness_due_to_soluble_aitken_mode','atmosphere_optical_thickness_due_to_soluble_accumulation_mode','atmosphere_optical_thickness_due_to_soluble_coarse_mode', 'atmosphere_optical_thickness_due_to_insoluble_aitken_mode',
             '10m_Wind_u','10m_wind_v','1.5m_Temp',
             'SO2_mmr','H2SO4_MMR','MSA_mmr','DMSO_MMR','NUCLEATION_MODE_SOL_H2SO4_MMR','AITKEN_MODE_SOL_H2SO4_MMR',
'ACCUMULATION_MODE_SOL_H2SO4_MMR','rm ','CDNC_Third','CDNC ',
             'H2SO4_to_aitken', 'H2SO4_to_accum','H2SO4_to_coarse',
'Cond_H2SO4_to_nucleation_soluble','Cond_H2SO4_to_aitken_soluble','Cond_H2SO4_to_accum_soluble','Cond_H2SO4_to_coarse_soluble',
'Cond_H2SO4_to_aitken_insoluble','Ccon_H2SO4_to_accum_insoluble','Cond_H2SO4_to_coarse_insoluble','CCN_number_concentration', 
'NUCLEATION_MODE_SOL_H2SO4','AITKEN_MODE_SOL_H2SO4','ACCUMULATION_MODE_SOL_H2SO4', 'COARSE_MODE_SOL_H2SO4',
             'DMS_surface_emissions','so2_surface_emissions','so2_high_lev_emissions','so2_natural_emissions',
'ozone_column']
## Compiling it into a dictionary. so: OCEAN_SURFACE_CHLOR = m01s00i096. So when extracting from the
## raw model data, it will take that variable (m01s00i096) and name it to OCEAN_SURFACE_CHLOR.
Monthly = {}
for v,m in zip(Monthly_Variables,M_msi):
    Monthly[v]=m
Daily = {}
for v,m in zip(Daily_Variables,D_msi):
    Daily[v]=m


Usage_Profile = [monthly_file_type,daily_file_type]
# Magic.
def monthly_iris(x): # Creates NC files from iris.
    constraint=iris.AttributeConstraint(STASH=Monthly[x])
    cube = iris.load(filename,constraint)
    iris.save(cube, '{}/{}_{}_{}.nc'.format(Monthly_save_file,monthly_file_type,x,months[i]))
    print(('{}a.p{}{}{} & {}').format(suite,use,year,months[i],Monthly[x]))    
    print(cube)
    print('{}/{}_{}_{}.nc'.format(Monthly_save_file,monthly_file_type,x,months[i]))
    return monthly_iris     

def daily_iris(x): # Creates NC files from iris.
    constraint=iris.AttributeConstraint(STASH=Daily[x])
    cube = iris.load(filename,constraint)
    iris.save(cube,'{}/{}_{}_{}.nc'.format(Daily_save_file,daily_file_type,x,months[i]))
    print(('{}a.p{}{}{} & {}').format(suite,use,year,months[i],Daily[x]))    
    print(cube)
    print('{}/{}_{}_{}.nc'.format(Daily_save_file,daily_file_type,x,months[i]))

    return daily_iris    

## FOR THE FIRST TIME RUNNING, IT IS RECOMMENDED TO UNCOMMENT BELOW. This will save the Stash and names out!! for good reference
#os.chdir("/home/ybh10/Numpy_Array/")
#save_dict(Daily,'Daily_Stash_Output.npy')
#save_dict(Monthly,'Monthly_Stash_Output.npy')

#os.chdir("/nesi/nobackup/niwa02757/ybh10/cylc-run/u-"+suite+"/share/data/History_Data/")


# In[14]:



for use in (Usage_Profile):
    #     #for yr in range(0,np.size(year))
    os.chdir('/home/ybh10/cylc-run/u-{}/share/data/History_Data'.format(suite))
    for i in range(0,np.size(months)):
        filename=(('{}a.p{}{}').format(suite,use,months[i]))
        print(filename)
#        if use == 'm':
#            for x in (Monthly_Variables):
#                monthly_iris(x)
        if use == 'd':
            for x in (Daily_Variables):
                daily_iris(x)
            # function for daily.






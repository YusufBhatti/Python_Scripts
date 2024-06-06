"""
Created on Wed Oct  7 17:27:10 2020

@author: ybh10

Converts the RAW simulation output data from .pp files into netcdf files. 

"""
from six.moves import cPickle as pickle
import iris
import numpy as np
import os
import glob 
import xarray as xr
import time
#time.sleep(14400)
os.chdir("/home/ybh10/Scripts")
#from my_functions import *
###!! THINGS NEEDED TO CHANGE WHEN CHANGING SUITE ID:
#suites=['cx359','cx381','cx740','cx579','cx741','cx361','ct130'] # the suite ID BUT NOT INCLUDING THE u-
suites=['co939'] # the suite ID BUT NOT INCLUDING THE u-

#suite='cv821'
 # Directory to save daily variables
for suite in (suites):
# Monthly_save_file='/home/ybh10/DMS_Emissions/Chemistry_Scheme/CHEM2.1/Monthly' # Directory to save monthly variables
##!! Change the months depending on if the model starts on Jan - Dec, or if just one month... etc.
#months = 
    months = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']

    control_suites=['co872','co939','cp001']
    perturbed_suites=['cp003','cp010','cp009','co803','co925','co940','co966','co967','cp002','cp006','cw469']
    CHEM3_MOD=['cp888','cq392','cq418','cq419','cq420','cq421','cq446','cq502','cq540','cq572','cr037','cr038','cr039','cr040','cr041','cr042']
    Lana_suites=['cq362','cq448']
    Lana_Chem3_suites=['cr147','cr152']
    CHEM3_MOD_New=['cr233','cr248','cr249','cr251','cr252','cr253','cr244','cr250','cr254','cr256']
    Lana_Chem3_new_suites=['cr298','cr437']
    Lana_Chem3_MSIA_new_suites=['cr356']
    REV3_suites=['cr687','cr721']
    MEDUSA_B17=['cw911']
    REV3_B17=['cw909']
    REV3_B17_SOCOL=['cy056']

    Mod_N00=['cr799','cr803','cr826','cr827','cr828','cr854','cr855','cr856','cr860','cr876']
    Mod_W14=['cr877','cr878']
    Mod_GM16=['cr982','cs138','cs061']
    Mod_B17=['ct129','ct130','cx362','cx363','cx364']

    Lana_W14=['cs715','']
    Lana_GM16=['cs716','']
    Lana_GM16_test=['cs877']
    Lana_B17=['ct131','ct132']
    Lana_LM86=['cq362','cq448','cw432']

    MODIS_B17_CHEM3=['cv395','cx358','cx360','cx361']
    MODIS_B17_CAMS6=['cv396']
    MODIS_B17_SOCOL=['cw941','cx269','cx381']
    MODIS_B17_NORESM=['cw987','cx359','cx366','cx366','cx367']
    MODIS_B17_MIROC=['cx566','cx741'] # 741 = real one
    MODIS_B17_GDFL=['cx578','cx740'] # 740 = real one
    MODIS_B17_GEOS=['cx579']

    MODIS_B17_Jaegle=['cv821']

    MODIS_B17_CLIM=['cx515','cx517'] #  2015 - 18, 2009-15

    if suite in control_suites:
        file_type='MEDUSA_LM86'
        print(file_type)
    if suite in perturbed_suites:
        file_type='MODIS_LM86'
        print(file_type)
    if suite in CHEM3_MOD:
        file_type='CHEM3_MOD'
        print(file_type)
    if suite in Lana_suites:
        file_type='Lana_LM86'
        print(file_type)
    if suite in Lana_Chem3_suites:
        file_type='Lana_CHEM3'
        print(file_type)
    if suite in CHEM3_MOD_New:
        file_type='CHEM3_MOD_New'
        print(file_type)
    if suite in Lana_Chem3_new_suites:
        file_type='Lana_CHEM3_NEW'
        print(file_type)
    if suite in Lana_Chem3_MSIA_new_suites:
        file_type='Lana_CHEM3_MSIA'
        print(file_type)
    if suite in REV3_suites:
        file_type='Hulswar_LM86'
    if suite in REV3_B17:
        file_type='Hulswar_B17'
    if suite in MEDUSA_B17:
        file_type='MEDUSA_B17'
    if suite in REV3_B17_SOCOL:
        file_type='Hulswar_B17_SOCOL'


        print(file_type)
    if suite in Mod_N00:
        file_type='MODIS_N00'
        print(file_type)
    if suite in Mod_W14:
        file_type='MODIS_W14'
        print(file_type)
    if suite in Mod_GM16:
        file_type='MODIS_GM16'
        print(file_type)
    if suite in Mod_B17:
        file_type='MODIS_B17'
        print(file_type)
    if suite in MODIS_B17_CHEM3:
        file_type='MODIS_B17_CHEM3'
        print(file_type)
    if suite in MODIS_B17_CAMS6:
        file_type='CAM6'
        print(file_type)
    if suite in MODIS_B17_Jaegle:
        file_type='MODIS_B17_Jaegle'
        print(file_type)
    if suite in MODIS_B17_SOCOL:
        file_type='MODIS_B17_SOCOL'
    if suite in MODIS_B17_NORESM:
        file_type='MODIS_B17_NORESM'
    if suite in MODIS_B17_MIROC:
        file_type='MODIS_B17_MIROC'
    if suite in MODIS_B17_GDFL:
        file_type='MODIS_B17_GFDL'
    if suite in MODIS_B17_GEOS:
        file_type='MODIS_B17_GEOS_CHEM'



    if suite in MODIS_B17_CLIM:
        file_type='MODIS_B17_CLIM'
        print(file_type)



    if suite in Lana_W14:
        file_type='Lana_W14'
        print(file_type)
    if suite in Lana_GM16:
        file_type='Lana_GM16'
        print(file_type)
    if suite in Lana_B17:
        file_type='Lana_B17'
        print(file_type)
    if suite in Lana_LM86:
        file_type='Lana_LM86'
        print(file_type)

    # if suite in MODIS_B17_CAMS6:
    #     file_type='MODIS_B17_CAMS6'
    #     print(file_type)



    run_length='Years'

    year=np.arange(2019,2020,1) 
    dates=[]
    for y in (year):
        for m in (months):
            dates.append('{}{}'.format(y,m))
    
    dates=dates[3:]
    #months = ['sep','oct'] # for one month. But must get rid of the [i] in the actual code: filename=(('{}a.p{}{}{}').format(suite,use,year,months[i]))
    year='2009' ## Important - this is what year your suite-id states. e.g.  cb134a.pm1988sep so here 1988.
    #year=np.arange(2003,2010,1) # for multiple years.
    ## This is to seperate the daily and monthly variables and directories appropriately. The letter here
    ## is what the suite-id states. e.g. cb134a.pm1988sep here is 'm' - for monthly. cb134a.pd1988sep - 'd' for daily
    daily_file_type='d'
    monthly_file_type='e'

    #os.chdir("/nesi/nobackup/niwa02757/ybh10/cylc-run/u-"+suite+"/share/data/History_Data/")
    section = (1,2,34,38,50) # This gathers the SECTIONS AND ITEMS for the stash requests - in ORDER.

    ## Change as appropriate, but MUST correspond to the 'Daily_Variables' section below. E.g. 
    ## DMS_CONCENTRATION_IN_SEAWATER = 132. Ditem = Daily variables. Mitem = Monthly variables. Some are the 
    ## same, some are different. It depends on what variable you want to extract for daily/monthly.
    Ditem =['209','208','241','298','299'],\
    ['205','206','285','300','301','302','303'],\
    ['001','071','072','073','102','104','108','114','111','117','101','103','107','113','119'],\
    ['439','485','486','487','488','516','517','518','519','520','539','504','505','506','507','508','509','510','401','402','403','404','405'],\
    ['214','215']
    #['135','136','137','138','139','140','141','142','143','144','145','146','150','151','152','153','214','215','217','219','331','332','333','334','335','336','337','338','339'],\
                #['135','136','137','138','139','140','141','142','143','144','145','146','150','151','152','153','331','332','333','334','335','336','337','338','339']

    Mitem =['031'],\
    ['075','145'],\
    ['217','219'],\
    ['003']


    msection = (0,34,50,34) # This gathers the SECTIONS AND ITEMS for the stash requests - in ORDER.

    # Mitem = ['096','132']
    # 'm01s34i003','m01s34i041','m01s34i081'
    Ditems=np.array(Ditem)
    Mitems=np.array(Mitem)

    x=-1
    ## This section compiles the stash requests in order to create a fully formatted stash (e.g. m01s00i096)
    D_msi=[]
    M_msi=[]
    x=-1
    for item,s in zip(Ditems,section):
        x=x+1
        for it in (item):
            a=('m01s{:02d}i{}').format(section[x],it)
            D_msi.append(a)
    x=-1

    for mmitem,ms in zip(Mitems,msection):
        x=x+1
        for its in (mmitem):
            a=('m01s{:02d}i{}').format(msection[x],its)
            M_msi.append(a)
    ###############################################################################################################################
    ## This section labels what the stash requests ACTUAL lsname is. This section is NOT FUNDEMENTAL, but very
    ## useful if you want to know what each variable extracted is named.If changing one of these names, make sure
    ## to change the Stash Request (m01s##i###) too! They correspond to each other!
    ###############################################################################################################################
    # --------------------------------------------->   # m01s00                                                                                                  # m01s00
    Daily_Variables = [                                
     # --------------------------------------------->   # m01s01                                                                                                 # m01s01
                 'CLEAR_SKY_OUTGOING_SW_FLUX_TOA','ALL_SKY_OUTGOING_SW_FLUX_TOA',                              
    'DROPLET_NUMBER_CONCENTRATION','CDNC_Cloud_Top','WEIGHT_FOR_CDNC',
     # --------------------------------------------->   # m01s02                                                                                                  # m01s02
                 'ALL_SKY_OUTGOING_LW_FLUX_TOA','CLEAR_SKY_OUTGOING_LW_FLUX_TOA',                           

    'atmosphere_optical_thickness_due_to_dust_ambient_aerosol','atmosphere_optical_thickness_due_to_soluble_aitken_mode_sulphate_aerosol',

    'atmosphere_optical_thickness_due_to_soluble_accumulation_mode_sulphate_aerosol','atmosphere_optical_thickness_due_to_soluble_coarse_mode_sulphate_aerosol', 

    'atmosphere_optical_thickness_due_to_insoluble_aitken_mode_sulphate_aerosol',                                                                    
    # --------------------------------------------->   # m01s03                                                                                                   # m01s03

    # --------------------------------------------->   # m01s34                                                                                                   # m01s34 
                 'O3_Mass_Mixing_Ratio','DMS_Mass_Mixing_Ratio','SO2_MASS_MIXING_RATIO','H2SO4_Mass_Mixing_Ratio','NUCLEATION_MODE_SOL_H2SO4_MMR','AITKEN_MODE_SOL_H2SO4_MMR',                                   # m01s34
    'ACCUMULATION_MODE_SOL_H2SO4_MMR','COARSE_MODE_SOL_H2SO4_MMR','ACCUMULATION_MODE_SOL_SSA_MMR','COARSE_MODE_SOL_SSA_MMR',    'number_mixing_ratio_nucleation_sol_mode','number_mixing_ratio_aitken_sol_mode','number_mixing_ratio_accumulation_sol_mode','number_mixing_ratio_coarse_sol_mode','number_mixing_ratio_aitken_insol_mode',
                                                                                                    # m01s34
    # --------------------------------------------->   # m01s38
                 'ccn_concentration_25r',
    'H2SO4_nucleation_sol','H2SO4_aitken_sol','H2SO4_accum_sol','H2SO4_coarse_sol','H2SO4_nucleation_load','H2SO4_aitken_load','H2SO4_accum_load','H2SO4_coarse_load','Total_H2SO4_load','Total_SSA_load',     'Number_nucleation_sol_mode','Number_aitken_sol_mode','Number_accumulation_sol_mode','Number_coarse_sol_mode','Number_aitken_insol_mode',
                          'Number_accumulation_insol_mode','Number_coarse_insol_mode','DRY_PARTICLE_DIAMETER_nucleation_sol', 'DRY_PARTICLE_DIAMETER_aitken_sol', 'DRY_PARTICLE_DIAMETER_accumulation_sol',
              'DRY_PARTICLE_DIAMETER_coarse_sol', 'DRY_PARTICLE_DIAMETER_aitken_insol',
        # --------------------------------------------->   # m01s50

            'DMS_SURF_EMISSIONS','SO2_SURF_EMISSIONS']

                       # --------------------------------------------->   # m01s50

    ###############################################################################################################################
    ###############################################################################################################################
    ## Compiling it into a dictionary. so: OCEAN_SURFACE_CHLOR = m01s00i096. So when extracting from the
    ## raw model data, it will take that variable (m01s00i096) and name it to OCEAN_SURFACE_CHLOR.
    Monthly_Variables = ['Sea_Ice_Fraction',                                
     # --------------------------------------------->   # m01s01                                                                                                 # m01s01

     # --------------------------------------------->   # m01s02                                                                                                  # m01s02

    # --------------------------------------------->   # m01s03                                                                                                   # m01s03

    # --------------------------------------------->   # m01s34                                                                                                   # m01s34 
    'DMSO_Mass_Mixing_Ratio','MSIA_Mass_Mixing_Ratio',                                                                                                # m01s34
    # --------------------------------------------->   # m01s38
                 'SO2_natural_emissions','OZONE_COLUMN_IN_DOBSON_UNITS',
                        'NO3_Mass_Mixing_Ratio']


    Daily = {}
    Monthly={}
    for m,v in zip(Daily_Variables,D_msi):
        #print(v,m)
        Daily[m]=v
    for mo,vo in zip(Monthly_Variables,M_msi):
        #print(v,m)
        Monthly[mo]=vo

    Usage_Profile = [monthly_file_type,daily_file_type]
    # Magic.

    print(file_type)
    def monthly_iris(x): # Creates NC files from iris.
        constraint=iris.AttributeConstraint(STASH=Monthly[x])
        cube = iris.load(filename,constraint)
        iris.save(cube,'{}/{}_{}.nc'.format(Output_file_monthly,x,dates[i]))
        print(('{}a.p{}{}{} & {}').format(suite,use,year,dates[i],Monthly[x]))    
        print(cube)
        print('{}/{}_{}.nc'.format(Output_file_daily,x,dates[i]))

        return monthly_iris     

    def daily_iris(x): # Creates NC files from iris.
        constraint=iris.AttributeConstraint(STASH=Daily[x])
        cube = iris.load(filename,constraint)
        iris.save(cube,'{}/{}_{}.nc'.format(Output_file_daily,x,dates[i]))
        print(('{}a.p{}{}{} & {}').format(suite,use,year,dates[i],Daily[x]))    
        print(cube)
        print('{}/{}_{}.nc'.format(Output_file_daily,x,dates[i]))

        return daily_iris    

    def hourly_iris(x): # Creates NC files from iris.
        constraint=iris.AttributeConstraint(STASH=Daily[x])
        cube = iris.load(filename,constraint)
        iris.save(cube,'{}/{}_{}.nc'.format(Output_file_daily,x,count))
       # iris.save(cube,'{}/../{}_{}.nc'.format(Output_file_daily,x,count))

        #rint(('{}a.p{}{} & {}').format(suite,year,filename[-11:-3],Daily[x]))    
        #rint(cube)
        print('{}/{}_{}.nc'.format(Output_file_daily,x,filename[-11:-3]))

        return hourly_iris    
    ###############################################################################################################################
    ## FOR THE FIRST TIME RUNNING, IT IS RECOMMENDED TO UNCOMMENT BELOW. This will save the Stash and names out!! for good reference
    ###############################################################################################################################
    #os.chdir("/home/ybh10/Numpy_Array/")
    #save_dict(Daily,'Daily_Stash_Output.npy')
    #save_dict(Monthly,'Monthly_Stash_Output.npy')
    ###############################################################################################################################

    ###############################################################################################################################
    #################################### FILE SAVE AND OUTPUT #####################################################################
    ###############################################################################################################################
    Output_file_daily='/home/ybh10/Objective_3/Postprocess_Data/{}/Raw'.format(file_type)
    Output_file_monthly='/home/ybh10/Objective_3/Postprocess_Data/{}/Raw'.format(file_type)
    
    Monthly_Variables=Monthly_Variables[5:]
    end
    for i in range(0,np.size(dates)):
    #     if i == 1:
    #         end
        for use in (Usage_Profile):
            os.chdir('/home/ybh10/cylc-run/u-{}/share/data/History_Data'.format(suite))

            filename=(('{}a.p{}{}.pp').format(suite,use,dates[i]))

            print(filename)
            if use == 'e':
                for x in (Monthly_Variables[:]):

                    #print(x)
                    try:
                        xr.open_dataset('{}/{}_{}.nc'.format(Output_file_daily,x,dates[i]))
                        pass
                    except:
                        try:
                            monthly_iris(x)
                        except:
                            pass
            if use == 'd':
                #end
                #            if i == 1:
                #    end
                for x in (Daily_Variables[:]):

                    print(x)
                    try:
                        xr.open_dataset('{}/{}_{}.nc'.format(Output_file_daily,x,dates[i])) 
                        pass
                    except:
                        try:
                            daily_iris(x)

                        except:
                            pass
                    

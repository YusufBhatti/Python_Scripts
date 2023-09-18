#!/usr/bin/env python2
# -*- coding: utf-8 -*-
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
        os.chdir('/home/ybh10/cylc-run/u-bx043/share/data/History_Data/')
        filename = 'bx043a.pm2003'+months[i]+'.pp'
        vars = iris.load(filename)
        #print vars
        
        os.chdir('/nesi/nobackup/niwa02757/ybh10/DMS/Lana/2003')    
        toa_outgoing_SW_flux_clear_sky = vars[480] ; iris.save(toa_outgoing_SW_flux_clear_sky,"bx_toa_outgoing_SW_flux_Clear_Sky_2003_{}.nc".format(months[i]))
        toa_outgoing_SW_flux = vars[478] ; iris.save(toa_outgoing_SW_flux,"bx_toa_outgoing_SW_flux_2003_{}.nc".format(months[i]))
        toa_outgoing_LW_flux_clear_sky = vars[476] ; iris.save(toa_outgoing_LW_flux_clear_sky,"bx_toa_outgoing_LW_flux_Clear_Sky_2003_{}.nc".format(months[i]))
        toa_outgoing_LW_flux = vars[474] ; iris.save(toa_outgoing_LW_flux,"bx_toa_outgoing_LW_flux_2003_{}.nc".format(months[i]))
        atmosphere_optical_thickness_due_to_soluble_aitken_mode_sulphate_aerosol = vars[23]; iris.save(atmosphere_optical_thickness_due_to_soluble_aitken_mode_sulphate_aerosol,"bx_atmosphere_optical_thickness_due_to_soluble_aitken_mode_sulphate_aerosol_2003_{}.nc".format(months[i]))
        atmosphere_optical_thickness_due_to_soluble_accumulation_mode_sulphate_aerosol = vars[24]; iris.save(atmosphere_optical_thickness_due_to_soluble_accumulation_mode_sulphate_aerosol,"bx_atmosphere_optical_thickness_due_to_soluble_accumulation_mode_sulphate_aerosol_2003_{}.nc".format(months[i]))       
        atmosphere_optical_thickness_due_to_soluble_coarse_mode_sulphate_aerosol = vars[25]; iris.save(atmosphere_optical_thickness_due_to_soluble_coarse_mode_sulphate_aerosol,"bx_atmosphere_optical_thickness_due_to_soluble_coarse_mode_sulphate_aerosol_2003_{}.nc".format(months[i]))        
        atmosphere_optical_thickness_due_to_insoluble_aitken_mode_sulphate_aerosol = vars[26]; iris.save(atmosphere_optical_thickness_due_to_insoluble_aitken_mode_sulphate_aerosol,"bx_atmosphere_optical_thickness_due_to_insoluble_aitken_mode_sulphate_aerosol_2003_{}.nc".format(months[i]))
        atmosphere_optical_thickness_due_to_dust_ambient_aerosol = vars[342]; iris.save(atmosphere_optical_thickness_due_to_dust_ambient_aerosol,"bx_atmosphere_optical_thickness_due_to_dust_ambient_aerosol_2003_{}.nc".format(months[i]))
        
        atmosphere_cloud_liquid_water_content = vars[337]; iris.save(atmosphere_cloud_liquid_water_content,"bx_atmosphere_cloud_liquid_water_content_2003_{}.nc".format(months[i]))
        atmosphere_cloud_ice_content = vars[336]; iris.save(atmosphere_cloud_ice_content,"bx_atmosphere_cloud_ice_content_2003_{}.nc".format(months[i]))
        cloud_volume_fraction_in_atmosphere_layer = vars[197]; iris.save(cloud_volume_fraction_in_atmosphere_layer,"bx_cloud_volume_fraction_in_atmosphere_layer_2003_{}.nc".format(months[i]))
        cloud_area_fraction  = vars[346]; iris.save(cloud_area_fraction ,"bx_cloud_area_fraction _2003_{}.nc".format(months[i]))
        mass_fraction_of_cloud_liquid_water_in_air = vars[370]; iris.save(mass_fraction_of_cloud_liquid_water_in_air,"bx_mass_fraction_of_cloud_liquid_water_in_air_2003_{}.nc".format(months[i]))
        
        product_of_effective_radius_of_cloud_liquid_water_particle_and_cloud_liquid_water_area_fraction_exposed_to_space_and_sunlit_binary_mask  = vars[1]; iris.save(product_of_effective_radius_of_cloud_liquid_water_particle_and_cloud_liquid_water_area_fraction_exposed_to_space_and_sunlit_binary_mask ,"bx_product_of_effective_radius_of_cloud_liquid_water_particle_and_cloud_liquid_water_area_fraction_exposed_to_space_and_sunlit_binary_mask_2003_{}.nc".format(months[i]))
        product_of_cloud_liquid_water_area_fraction_exposed_to_space_and_sunlit_binary_mask  = vars[2]; iris.save(product_of_cloud_liquid_water_area_fraction_exposed_to_space_and_sunlit_binary_mask  ,"bx_product_of_cloud_liquid_water_area_fraction_exposed_to_space_and_sunlit_binary_mask_2003_{}.nc".format(months[i]))
        
        mass_fraction_of_sulfuric_acid_in_soluble_nucleation_mode_dry_aerosol_in_air = vars[243]; iris.save(mass_fraction_of_sulfuric_acid_in_soluble_nucleation_mode_dry_aerosol_in_air,"bx_mass_fraction_of_sulfuric_acid_in_soluble_nucleation_mode_dry_aerosol_in_air_2003_{}.nc".format(months[i]))   
        mass_fraction_of_sulfuric_acid_in_soluble_aitken_mode_dry_aerosol_in_air = vars[245]; iris.save(mass_fraction_of_sulfuric_acid_in_soluble_aitken_mode_dry_aerosol_in_air,"bx_mass_fraction_of_sulfuric_acid_in_soluble_aitken_mode_dry_aerosol_in_air_2003_{}.nc".format(months[i]))   
        mass_fraction_of_sulfuric_acid_in_soluble_accumulation_mode_dry_aerosol_in_air = vars[249]; iris.save(mass_fraction_of_sulfuric_acid_in_soluble_accumulation_mode_dry_aerosol_in_air,"bx_mass_fraction_of_sulfuric_acid_in_soluble_accumulation_mode_dry_aerosol_in_air_2003_{}.nc".format(months[i]))   
        mass_fraction_of_sulfuric_acid_in_soluble_coarse_mode_dry_aerosol_in_air = vars[254]; iris.save(mass_fraction_of_sulfuric_acid_in_soluble_coarse_mode_dry_aerosol_in_air,"bx_mass_fraction_of_sulfuric_acid_in_soluble_coarse_mode_dry_aerosol_in_air_2003_{}.nc".format(months[i]))   
        mass_fraction_of_dimethyl_sulfide_in_air = vars[372]; iris.save(mass_fraction_of_dimethyl_sulfide_in_air,"bx_mass_fraction_of_dimethyl_sulfide_in_air_2003_{}.nc".format(months[i]))   
        mass_fraction_of_ozone_in_air = vars[393]; iris.save(mass_fraction_of_ozone_in_air,"bx_mass_fraction_of_ozone_in_air_2003_{}.nc".format(months[i]))   
        mass_fraction_of_sulfur_dioxide_in_air = vars[397]; iris.save(mass_fraction_of_sulfur_dioxide_in_air,"bx_mass_fraction_of_sulfur_dioxide_in_air_2003_{}.nc".format(months[i]))   
        mass_fraction_of_sulfuric_acid_in_air = vars[398]; iris.save(mass_fraction_of_sulfuric_acid_in_air,"bx_mass_fraction_of_sulfuric_acid_in_air_2003_{}.nc".format(months[i]))   
        mole_concentration_of_dimethyl_sulfide_in_sea_water = vars[401]; iris.save(mole_concentration_of_dimethyl_sulfide_in_sea_water ,"bx_mole_concentration_of_dimethyl_sulfide_in_sea_water _2003_{}.nc".format(months[i]))   
        air_temperature = vars[334]; iris.save(air_temperature,"bx_air_temperature_2003_{}.nc".format(months[i]))   
        upward_air_velocity = vars[484]; iris.save(upward_air_velocity,"bx_upward_air_velocity_2003_{}.nc".format(months[i]))   
        wind_speed = vars[488]; iris.save(wind_speed,"bx_wind_speed_2003_{}.nc".format(months[i]))   
        y_wind  = vars[495]; iris.save(y_wind ,"bx_y_wind _2003_{}.nc".format(months[i]))   
        x_wind = vars[492]; iris.save(x_wind,"bx_x_wind_2003_{}.nc".format(months[i]))   

        os.chdir('/home/ybh10/cylc-run/u-bw945/share/data/History_Data/')       
        filename = 'bw945a.pm2003'+months[i]+'.pp'
        vars = iris.load(filename)
        os.chdir('/nesi/nobackup/niwa02757/ybh10/DMS/MEDUSA/2003')       
        
        product_of_effective_radius_of_cloud_liquid_water_particle_and_cloud_liquid_water_area_fraction_exposed_to_space_and_sunlit_binary_mask  = vars[1]; iris.save(product_of_effective_radius_of_cloud_liquid_water_particle_and_cloud_liquid_water_area_fraction_exposed_to_space_and_sunlit_binary_mask,"bw_product_of_effective_radius_of_cloud_liquid_water_particle_and_cloud_liquid_water_area_fraction_exposed_to_space_and_sunlit_binary_mask_2003_{}.nc".format(months[i]))
        product_of_cloud_liquid_water_area_fraction_exposed_to_space_and_sunlit_binary_mask  = vars[2]; iris.save(product_of_cloud_liquid_water_area_fraction_exposed_to_space_and_sunlit_binary_mask  ,"bw_product_of_cloud_liquid_water_area_fraction_exposed_to_space_and_sunlit_binary_mask_2003_{}.nc".format(months[i]))

        atmosphere_optical_thickness_due_to_soluble_aitken_mode_sulphate_aerosol = vars[22]; iris.save(atmosphere_optical_thickness_due_to_soluble_aitken_mode_sulphate_aerosol,"bw_atmosphere_optical_thickness_due_to_soluble_aitken_mode_sulphate_aerosol_2003_{}.nc".format(months[i]))
        atmosphere_optical_thickness_due_to_soluble_accumulation_mode_sulphate_aerosol = vars[23]; iris.save(atmosphere_optical_thickness_due_to_soluble_accumulation_mode_sulphate_aerosol,"bw_atmosphere_optical_thickness_due_to_soluble_accumulation_mode_sulphate_aerosol_2003_{}.nc".format(months[i]))  

        atmosphere_optical_thickness_due_to_soluble_coarse_mode_sulphate_aerosol = vars[24]; iris.save(atmosphere_optical_thickness_due_to_soluble_coarse_mode_sulphate_aerosol,"bw_atmosphere_optical_thickness_due_to_soluble_coarse_mode_sulphate_aerosol_2003_{}.nc".format(months[i]))        
        atmosphere_optical_thickness_due_to_insoluble_aitken_mode_sulphate_aerosol = vars[25]; iris.save(atmosphere_optical_thickness_due_to_insoluble_aitken_mode_sulphate_aerosol,"bw_atmosphere_optical_thickness_due_to_insoluble_aitken_mode_sulphate_aerosol_2003_{}.nc".format(months[i]))
        atmosphere_optical_thickness_due_to_dust_ambient_aerosol = vars[328]; iris.save(atmosphere_optical_thickness_due_to_dust_ambient_aerosol,"bw_atmosphere_optical_thickness_due_to_dust_ambient_aerosol_2003_{}.nc".format(months[i]))

        atmosphere_cloud_liquid_water_content = vars[325]; iris.save(atmosphere_cloud_liquid_water_content,"bw_atmosphere_cloud_liquid_water_content_2003_{}.nc".format(months[i]))
        atmosphere_cloud_ice_content = vars[324]; iris.save(atmosphere_cloud_ice_content,"bw_atmosphere_cloud_ice_content_2003_{}.nc".format(months[i]))
        cloud_area_fraction  = vars[332]; iris.save(cloud_area_fraction ,"bx_cloud_area_fraction _2003_{}.nc".format(months[i]))
        cloud_volume_fraction_in_atmosphere_layer = vars[195]; iris.save(cloud_volume_fraction_in_atmosphere_layer,"bw_cloud_volume_fraction_in_atmosphere_layer_2003_{}.nc".format(months[i]))

        mass_fraction_of_sulfuric_acid_in_soluble_nucleation_mode_dry_aerosol_in_air = vars[241]; iris.save(mass_fraction_of_sulfuric_acid_in_soluble_nucleation_mode_dry_aerosol_in_air,"bw_mass_fraction_of_sulfuric_acid_in_soluble_nucleation_mode_dry_aerosol_in_air_2003_{}.nc".format(months[i]))   
        mass_fraction_of_sulfuric_acid_in_soluble_aitken_mode_dry_aerosol_in_air = vars[243]; iris.save(mass_fraction_of_sulfuric_acid_in_soluble_aitken_mode_dry_aerosol_in_air,"bw_mass_fraction_of_sulfuric_acid_in_soluble_aitken_mode_dry_aerosol_in_air_2003_{}.nc".format(months[i]))   
        mass_fraction_of_sulfuric_acid_in_soluble_accumulation_mode_dry_aerosol_in_air = vars[247]; iris.save(mass_fraction_of_sulfuric_acid_in_soluble_accumulation_mode_dry_aerosol_in_air,"bw_mass_fraction_of_sulfuric_acid_in_soluble_accumulation_mode_dry_aerosol_in_air_2003_{}.nc".format(months[i]))   
        mass_fraction_of_sulfuric_acid_in_soluble_coarse_mode_dry_aerosol_in_air = vars[252]; iris.save(mass_fraction_of_sulfuric_acid_in_soluble_coarse_mode_dry_aerosol_in_air,"bw_mass_fraction_of_sulfuric_acid_in_soluble_coarse_mode_dry_aerosol_in_air_2003_{}.nc".format(months[i]))   
        mass_fraction_of_dimethyl_sulfide_in_air = vars[358]; iris.save(mass_fraction_of_dimethyl_sulfide_in_air,"bx_mass_fraction_of_dimethyl_sulfide_in_air_2003_{}.nc".format(months[i]))   
        mass_fraction_of_ozone_in_air = vars[379]; iris.save(mass_fraction_of_ozone_in_air,"bw_mass_fraction_of_ozone_in_air_2003_{}.nc".format(months[i]))   
        mass_fraction_of_sulfur_dioxide_in_air = vars[383]; iris.save(mass_fraction_of_sulfur_dioxide_in_air,"bw_mass_fraction_of_sulfur_dioxide_in_air_2003_{}.nc".format(months[i]))   
        mass_fraction_of_sulfuric_acid_in_air = vars[384]; iris.save(mass_fraction_of_sulfuric_acid_in_air,"bw_mass_fraction_of_sulfuric_acid_in_air_2003_{}.nc".format(months[i]))   
        mole_concentration_of_dimethyl_sulfide_in_sea_water = vars[387]; iris.save(mole_concentration_of_dimethyl_sulfide_in_sea_water ,"bw_mole_concentration_of_dimethyl_sulfide_in_sea_water _2003_{}.nc".format(months[i]))   
        
        toa_outgoing_SW_flux_clear_sky = vars[466] ; iris.save(toa_outgoing_SW_flux_clear_sky,"bw_toa_outgoing_SW_flux_Clear_Sky_2003_{}.nc".format(months[i]))
        toa_outgoing_SW_flux = vars[464] ; iris.save(toa_outgoing_SW_flux,"bx_toa_outgoing_SW_flux_2003_{}.nc".format(months[i]))
        toa_outgoing_LW_flux_clear_sky = vars[462] ; iris.save(toa_outgoing_LW_flux_clear_sky,"bw_toa_outgoing_LW_flux_Clear_Sky_2003_{}.nc".format(months[i]))
        toa_outgoing_LW_flux = vars[459] ; iris.save(toa_outgoing_LW_flux,"bx_toa_outgoing_LW_flux_2003_{}.nc".format(months[i]))

        air_temperature = vars[332]; iris.save(air_temperature,"bw_air_temperature_2003_{}.nc".format(months[i]))   
        upward_air_velocity = vars[470]; iris.save(upward_air_velocity,"bw_upward_air_velocity_2003_{}.nc".format(months[i]))   
        wind_speed = vars[474]; iris.save(wind_speed,"bw_wind_speed_2003_{}.nc".format(months[i]))   
        y_wind  = vars[478]; iris.save(y_wind ,"bw_y_wind _2003_{}.nc".format(months[i]))   
        x_wind = vars[481]; iris.save(x_wind,"bw_x_wind_2003_{}.nc".format(months[i]))   

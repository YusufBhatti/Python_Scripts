# -*- coding: utf-8 -*-
"""
This code reads in a years worth of data 
in this case we are donwloading mean sea level presure dara
every 6 hours for a specified year. 
A for loop goes through all the ERA5 from 1980 to 2019
to collect all the relevant data.


READS ERA5 data from teh Copernicus DATA STORE

The details of the cdsapi package are at https://pypi.org/project/cdsapi/

Created on Sun Apr 19 13:49:26 2020

@author: ajm226
"""

import cdsapi
import numpy as np

year_array=np.arange(2003,2019) # identifying a range of years
for year in year_array:
    c = cdsapi.Client()  # opening up a client to interact with the CDS
    c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'product_type': 'reanalysis',
            'variable': 'mean_sea_level_pressure',
            'year': str(year),  # this variable is coming from the for loop
            'month': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
           ],
"""            'day': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
                '13', '14', '15',
                '16', '17', '18',
                '19', '20', '21',
                '22', '23', '24',
                '25', '26', '27',
                '28', '29', '30',
                '31',
            ],
  
            'time': [
                '00:00', '06:00', '12:00',
                '18:00',     
            ],  # this is UTC time """
            'format': 'netcdf',
        },
        'mslp'+str(year)+'.nc')   # this is where the data gets written to on your system

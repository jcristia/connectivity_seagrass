
# small run for testing purposes
# test sensitivity of particle locations to changes in the time step
# now using a reader resolution of 0.005

timestep = 5

#import sys
#sys.path.append("/Linux/src/opendrift-master")
import numpy as np
from datetime import datetime
from datetime import timedelta
import ogr
import os

from opendrift.models.oceandrift import OceanDrift
from opendrift.readers import reader_netCDF_CF_unstructured_Salish005
from opendrift.readers import reader_basemap_landmask

# temporary fix for when running on anything other than mank01:
#sys.path.append("/Linux/src/opendrift-master")

#####################################################

# load readers outside of loop
file_salish_cluster = r'D:\Hakai\models\salishsea\salishseacast_20160101_20160301\forcing\SalishSea_1h_20160101_20160301_opendrift.nc'

reader_salish = reader_netCDF_CF_unstructured_Salish005.Reader(file_salish_cluster)

reader_basemap = reader_basemap_landmask.Reader(
                       llcrnrlon=-142.0, llcrnrlat=42.0,
                       urcrnrlon=-121.0, urcrnrlat=60.6,
                       resolution='f', projection='merc')

#####################################################

o = OceanDrift(loglevel=0)
o.add_reader([reader_basemap, reader_salish])

# seed 12 particles randomly throughout entire extent
# do it without stranding
o.seed_elements(lon=[-125.1749, -125.0553, -124.9804, -124.5178, -123.912, -123.6736, -123.2168, -123.0741, -122.5614, -122.9538, -123.7274, -122.6249], lat=[50.09507, 50.0416, 49.78544, 49.5716, 49.34444, 49.14915, 49.04325, 48.65438, 48.58271, 48.30416, 48.17336, 47.1926], number=12, time=reader_salish.start_time)

o.set_config('drift:current_uncertainty', 0)
o.set_config('general:coastline_action', 'previous')
#o.set_config('general:coastline_action', 'stranding')
o.set_config('drift:scheme', 'euler')

o.run(end_time=reader_salish.start_time + timedelta(days=3), time_step=timestep, time_step_output=1800,
    outfile=r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_dev_scratch\sensitivity_testing\time_step\working_simulation\outputs/seagrass12particles_' + str(timestep) + 's.nc', export_variables=["age_seconds", "land_binary_mask"])   
    
print(o)

#####################################################
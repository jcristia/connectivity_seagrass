
# small run for testing purposes
# changing window resolution in the unstructured.
# change it manually in the reader. The value here is just for properly naming the output.

# WINDOW RES:
window_res = 0.005

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

# multiple opendrift runs
# iterate over seagrass shapefiles

sg_path = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_dev_scratch\sensitivity_testing\window_res\working_simulation\seagrass_split'
sg_files = os.listdir(sg_path)
shapefiles = []
for file in sg_files:
    if file.endswith('.shp'):
        shapefiles.append(os.path.join(sg_path, file))

for shp in shapefiles:

    # get base number for output names
    base = os.path.splitext(os.path.basename(shp))[0]

    # get number of particles to seed
    shp = ogr.Open(shp)
    lyr = shp.GetLayer(0)
    for feature in lyr:
        particles = feature.GetField('particles')
        break

    # REMOVE THIS LATER:
    particles = int(particles / 100)

    o = OceanDrift(loglevel=0)
    o.add_reader([reader_basemap, reader_salish])

    #time_step = timedelta(hours=1)
    #num_steps = 1
    #for i in range(num_steps):
    #    o.seed_from_shapefile(shp, number=particles, time=reader_salish.start_time + i*time_step)

    # seed 10 particles randomly throughout entire extent
    # do it without stranding
    o.seed_elements(lon=[-125.1749, -125.0553, -124.9804, -124.5178, -123.912, -123.6736, -123.2168, -123.0741, -122.5614, -122.9538, -123.7274, -122.6249], lat=[50.09507, 50.0416, 49.78544, 49.5716, 49.34444, 49.14915, 49.04325, 48.65438, 48.58271, 48.30416, 48.17336, 47.1926], number=12, time=reader_salish.start_time)


    o.set_config('drift:current_uncertainty', 0)
    o.set_config('general:coastline_action', 'previous')
    #o.set_config('general:coastline_action', 'stranding')
    o.set_config('drift:scheme', 'euler')

    #o.run(end_time=reader_salish.start_time + timedelta(days=3), time_step=30, time_step_output=1800,
    #      outfile=r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_dev_scratch\sensitivity_testing\window_res\working_simulation\outputs/seagrass_' + base + "_" + str(window_res).split('.')[1] + '.nc', export_variables=["age_seconds", "land_binary_mask"])

    o.run(end_time=reader_salish.start_time + timedelta(days=3), time_step=30, time_step_output=1800,
        outfile=r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_dev_scratch\sensitivity_testing\window_res\working_simulation\outputs/seagrass12particles_' + str(window_res).split('.')[1] + '.nc', export_variables=["age_seconds", "land_binary_mask"])   
    
    print(o)

    o.animation(filename=r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_dev_scratch\sensitivity_testing\window_res\working_simulation\outputs\seagrass12particles_' + str(window_res).split('.')[1] + '.gif')

    #o2 = OceanDrift(loglevel=0)
    #o2.io_import_file(r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_dev_scratch\sensitivity_testing\window_res\working_simulation\outputs\seagrass_sg1_01.nc')
    #o2.animation(filename=r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_dev_scratch\sensitivity_testing\window_res\working_simulation\outputs\seagrass_sg1_01.gif')
    #o.animation(compare=o2, legend=['0.0001 deg', '0.01 deg'], filename=r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_dev_scratch\sensitivity_testing\window_res\working_simulation\outputs\seagrass_sg1_compare_01_0001.gif')


#####################################################
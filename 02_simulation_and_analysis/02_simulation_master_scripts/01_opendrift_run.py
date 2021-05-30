# Opendrift run

part_fact = 1


import sys
sys.path.append("/Linux/src/opendrift-master")
import numpy as np
from datetime import datetime, timedelta
import ogr
import os

from opendrift.models.oceandrift import OceanDrift
from opendrift.readers import reader_netCDF_CF_unstructured_Salish005
from opendrift.readers import reader_basemap_landmask


#####################################################

# load readers outside of loop
file_salish_cluster = r'/home/jcristia/models/salishsea/SalishSea_1h_20170101_20170316_opendrift.nc'

reader_salish = reader_netCDF_CF_unstructured_Salish005.Reader(file_salish_cluster)

reader_basemap = reader_basemap_landmask.Reader(
                       llcrnrlon=-126.3, llcrnrlat=46.8,
                       urcrnrlon=-122.0, urcrnrlat=51.1,
                       resolution='f', projection='merc')

#####################################################

sg_path = r'seagrass_split'
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

    particles = int(particles * part_fact)

    o = OceanDrift(loglevel=0)
    o.add_reader([reader_basemap, reader_salish])

    time_step = timedelta(hours=4)
    num_steps = 84
    for i in range(num_steps):
        o.seed_from_shapefile(shp, number=particles, time=reader_salish.start_time + i*time_step)

    np.save('outputs/lon_' + base + '.npy', o.elements_scheduled.lon)
    np.save('outputs/lat_' + base + '.npy', o.elements_scheduled.lat)

    o.set_config('drift:current_uncertainty', 0.22)
    o.set_config('general:coastline_action', 'stranding')
    o.set_config('drift:scheme', 'euler')

    o.run(end_time=reader_salish.end_time, time_step=60, time_step_output=1800,
          outfile=r'outputs/seagrass_' + base + '.nc', export_variables=["age_seconds", "land_binary_mask"])

#####################################################
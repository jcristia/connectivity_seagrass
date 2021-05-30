# Average connections for all PLDs within one simulation

import os
import pandas as pd
import geopandas as gp

# input
shp_folder = r'D:\Hakai\script_runs\seagrass\seagrass_20200327_SS201408\shp_merged'
file = 'connectivity_pld{}.shp'
plds = ['01', '03', '07', '21', '60']
shp_out = 'connectivity_average.shp'

# read in shapefiles within folder as a geopanadas dataframe (gdf) and append to overall gdf
conn = []
for p in plds:
    conn.append(os.path.join(shp_folder, file.format(p)))
gdf_all = gp.GeoDataFrame()
for shp in conn:
    gdf = gp.read_file(shp)
    gdf_all = gdf_all.append(gdf)
gdf_all = gdf_all.astype({'from_id':int, 'to_id':int})

# groupby
# on aggregation, use a custom function
def mean_cust_denom(x):
    s = x.sum()
    m = s/float(len(plds))
    return m
gdf_group = gdf_all.groupby(['from_id', 'to_id']).agg(
    prob_avg = ('prob', mean_cust_denom),
    time_int = ('time_int', 'first'),
    totalori = ('totalori', 'first'),
    date_start = ('date_start', 'first'),
    geometry = ('geometry', 'first')
    )
gdf_group = gdf_group.astype({'time_int':int, 'totalori':int})
gdf_group = gdf_group.reset_index()

# output
gdf_f = gp.GeoDataFrame(gdf_group, crs=gdf.crs)
gdf_f.to_file(filename=os.path.join(shp_folder, shp_out), driver='ESRI Shapefile')
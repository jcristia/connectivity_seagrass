# Outputs an average of averages for connections and an average for dPC for nodes


import os
import pandas as pd
import geopandas as gp
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from math import log10, floor, sqrt
import statsmodels



###########################
# configuration
###########################

dirs = [
    'seagrass_20200228_SS201701',
    'seagrass_20200309_SS201705',
    'seagrass_20200309_SS201708',
    'seagrass_20200310_SS201101',
    'seagrass_20200310_SS201105',
    'seagrass_20200310_SS201108',
    'seagrass_20200327_SS201401',
    'seagrass_20200327_SS201405',
    'seagrass_20200327_SS201408',
    ]

root = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_runs_cluster\seagrass'
shp_merged = 'shp_merged'
conn_avg = 'connectivity_average.shp'
conn_ind = 'connectivity_pld{}.shp'
conn_conefor = r'patch_centroids_metrics_commavg.shp'
plds = ['01', '03', '07', '21', '60']
subfolders = 9
subfolder = 'seagrass_{}'
output_ind = 'output_figs'
output_all = os.path.join(root, 'output_figs_SALISHSEA_ALL')
overall_indices = r'conefor\conefor_connectivity_{}\overall_indices.txt'
grand_average = r'output_figs_SALISHSEA_ALL\connectivity_average_ALL.shp'


# create directories and paths
if not os.path.exists(output_all):
    os.mkdir(output_all)



########################


# overall connectivity

df = gp.GeoDataFrame()
for dir in dirs:
    conn = (os.path.join(root, dir, shp_merged, conn_avg))
    conn = gp.read_file(conn)
    df = df.append(conn)
df['date'] = df.date_start.str[:7]
df['dateYR'] = df.date.str[:4]
df['dateSEA'] = df.date.str[5:8]

def customstd(x):
    var = np.var(x, ddof=0)
    std = sqrt(var)
    return std
ldirs = len(dirs)
def customstd9(x):
    mean = sum(x)/ldirs
    var = (sum((x-mean)**2))/ldirs
    std = sqrt(var)
    return std
def conn_mean(x):
    s = x.sum()
    m = s/ldirs
    return m
df_freqstd = df.groupby(['from_id', 'to_id']).agg(
    freq = ('from_id', 'count'),
    probavgm = ('prob_avg', conn_mean),
    prob_stdf0 = ('prob_avg', customstd),
    prob_std9 = ('prob_avg', customstd9),
    geometry = ('geometry', 'first'),
    date = ('date', 'first'),
    dateYR = ('dateYR', 'first'),
    dateSEA = ('dateSEA', 'first'),
    timeintavg = ('time_int', 'mean'),
    timeintmin = ('time_int', 'min'),
    timeintmax = ('time_int', 'max')
    ).reset_index()

gdf = gp.GeoDataFrame(df_freqstd)
gdf.crs = df.crs
gdf.to_file(filename=os.path.join(output_all, 'connectivity_average_ALL.shp'), driver='ESRI Shapefile')



#########################
# Conefor values

gdf = gp.GeoDataFrame()
for dir in dirs:
    df = gp.read_file(os.path.join(root, dir, shp_merged, conn_conefor))
    gdf = gdf.append(df)
gdf['date'] = gdf.date_start.str[:7]
gdf['dateYR'] = gdf.date.str[:4]
gdf['dateSEA'] = gdf.date.str[5:8]

gdf_avg = gdf.groupby(['uID']).agg(
    dPC = ('dPC', 'mean'),
    dPCintra = ('dPCintra', 'mean'),
    dPCflux = ('dPCflux', 'mean'),
    dPCconnect = ('dPCconnect', 'mean'),
    dBC_PC = ('dBC_PC', 'mean'),
    dPCv = ('dPC', 'var'),
    dPCintrav = ('dPCintra', 'var'),
    dPCfluxv = ('dPCflux', 'var'),
    dPCconnecv = ('dPCconnect', 'var'),
    dBC_PCv = ('dBC_PC', 'var'),
    dPCstd = ('dPC', 'std'),
    dPCintrast = ('dPCintra', 'std'),
    dPCfluxstd = ('dPCflux', 'std'),
    dPCconnecs = ('dPCconnect', 'std'),
    dBC_PCstd = ('dBC_PC', 'std'),
    #comid_mode = ('comidns', pd.Series.mode),
    #comid_uniq = ('comidns', 'nunique'),  # IDs change each time, so this doesn't work.
    geometry = ('geometry', 'first'),
    date = ('date', 'first'),
    dateYR = ('dateYR', 'first'),
    dateSEA = ('dateSEA', 'first')
    ).reset_index()
gdf_avg['dPCinter'] = gdf_avg.dPCflux + gdf_avg.dPCconnect
gdf_avg = gp.GeoDataFrame(gdf_avg)
gdf_avg.crs = gdf.crs
gdf_avg.to_file(filename=os.path.join(output_all, 'patch_centroids_metrics_ALL.shp'), driver='ESRI Shapefile')




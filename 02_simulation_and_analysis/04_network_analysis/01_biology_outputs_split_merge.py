# Merge outputs from biology_opendrift.
# For runs where seagrass datasets were split up, there will be multiple output shapefiles that need to be merged together.

# NOTE: merging the destination pts is optional. This result in millions of particles.


import os
import geopandas
import math
import pandas as pd
from shapely.geometry import LineString
import shutil

######################
# User inputs
######################

# folder structure
root = r'D:\Hakai\script_runs\seagrass\{}\seagrass_{}\outputs\shp'
project_folder = 'seagrass_20200327_SS201408'

# number of sub folders:
subfolders = 9

# PLDs
plds = ['01', '03', '07', '21', '60']

# output folder for merged shapefiles
out_folder = r'D:\Hakai\script_runs\seagrass\{}\shp_merged'

# root names of input/output shapefiles
conn = 'connectivity_sg{}_pld{}.shp'
pts = 'dest_biology_pts_sg{}.shp'
centroids = 'patch_centroids.shp'

# set seagrass_og, seagrass_crs
seagrass_og = r'D:\Hakai\script_runs\seagrass\seagrass_20200228_SS201701\seagrass_1\seagrass_og'
seagrass_crs = {'init' :'epsg:3005'}

time_step_output = 0.5 # from biology script


######################
# MERGE
######################

# create output folder
out = out_folder.format(project_folder)
if not os.path.exists(out):
    os.mkdir(out)

# copy in centroids
# just need to do this once since each sub folder has the centroids of all meadows
if not os.path.isfile(os.path.join(out, centroids)):
    for file in os.listdir(root.format(project_folder, 1)):
        if file.startswith(centroids.split('.')[0]):
            shutil.copyfile(os.path.join(root.format(project_folder, 1), file), os.path.join(out, file))

# merge destination points
#out_shp = os.path.join(out, 'dest_biology_pts.shp')
#if not os.path.isfile(out_shp):
#    merge_list_dest = []
#    for sub in range(1, subfolders + 1):
#        dest_pt = os.path.join(root.format(project_folder, sub), pts.format(sub))
#        merge_list_dest.append(dest_pt)
#    arcpy.Merge_management(merge_list_dest, out_shp)




# create connection line shapefile
# Issue: for some simulations, the particles strand/settle before a certain PLD level is reached (e.g. everything strands by day 29, so there will not be a PLD60 dataset created).
# However, I still need a PLD60 dataset or else it looks like no connectivity happened for the range between 21 and 60.
# Therefore, create a dataset of all connections up to the max time step reached.
def connection_lines(shp_out, seagrass_og, seagrass_crs, conn_lines_out, date_start, pld_int, pld):

    print("writing connection lines to shapefile")
    od = geopandas.read_file(shp_out)
    sg = geopandas.read_file(seagrass_og)

    ### on od, select particles where time_int_s minus time_int is less than or equal to PLD
    od_pld = od[(od.time_int - od.time_int_s <= pld_int)]

    # get each unique combination of originID and destID and get count of particles that survived
    od_unique = od_pld[(od_pld.dest_id != -1) & (od_pld.mortstep == -1)].groupby(['uID','dest_id']).size().reset_index(name='Freq')
    # how to read this:
    # first we select from od the ones that settled and survived
    # then we groupby unique combinations of uID and dest_id
    # then we get the count of those unique combinations
    # this normally makes uID the index and doesn't have a column name for count (the series we created), so we reset index and give the count a column name

    # df of time interval where first settlement occurred
    df_time_int = od_pld[(od_pld.dest_id != -1) & (od_pld.mortstep == -1)].groupby(['uID','dest_id'])['time_int'].min().reset_index(name='time_int')
    
    # set up for creating self connection lines. Size of circle lines based on amount settled and average area of all patches.
    def CircleCoords(xLeft, yCenter, r, n): # CREDIT: MGET - Roberts and Treml
        return [(xLeft + r - math.cos(2*math.pi/n*x)*r, math.sin(2*math.pi/n*x)*r + yCenter) for x in range(n+1)]
    # min and max quantities used for normalization
    quantity_min = od_unique[od_unique.dest_id == od_unique.uID].Freq.min()
    quantity_max = od_unique[od_unique.dest_id == od_unique.uID].Freq.max()
    # get average area
    area_mean = sg.area.mean()
    # get radius of a circle with this area
    radius = math.sqrt(area_mean/math.pi)
    
    # for each unique combinaton create line from centroids
    connection_lines = pd.DataFrame(columns=['from_id','to_id','quantity','totalori','prob','time_int', 'line', 'pld'])
    conn_i = 0
    for row in od_unique.itertuples(index=False):
        # get total amount of particles released from patch
        total = od.uID[od.uID ==  row[0]].value_counts().values[0]

        # time interval where first settlement occurred
        time_int = df_time_int[(df_time_int.uID == row[0]) & (df_time_int.dest_id == row[1])]['time_int'].values[0]
    
        # get centroid of from and to patches
        centroid_origin = sg[sg.uID == row[0]].centroid
        centroid_dest = sg[sg.uID == row[1]].centroid

        if row[0] != row[1]:
            geom_line = LineString([centroid_origin.tolist()[0], centroid_dest.tolist()[0]])
        else:
            # normalize the quantites to 0.5 - 1 range (or I can do 0-1 but then the smallest one won't show up)
            #quantity_norm = 0.5 * (row[2] - quantity_min) / float(quantity_max - quantity_min) + 0.5
            quantity_norm = (row[2] - quantity_min) / float(quantity_max - quantity_min)
            radius_adj = radius * quantity_norm
            geom_line = LineString(CircleCoords(centroid_origin.x.tolist()[0], centroid_origin.y.tolist()[0], radius_adj, 90))
    
        connection_lines.loc[conn_i] = [row[0],row[1],float(row[2]),float(total),row[2]/float(total), time_int,geom_line, pld]
        conn_i += 1
    
    connection_lines['date_start'] = date_start   
    connection_lines = geopandas.GeoDataFrame(connection_lines, geometry='line')
    connection_lines.crs = seagrass_crs
    connection_lines.to_file(filename=conn_lines_out, driver='ESRI Shapefile')





# merge connectivity lines
for pld in plds:
    merge_list_conn = []
    out_shp_conn = os.path.join(out, 'connectivity_pld' + pld + '.shp')
    if not os.path.isfile(out_shp_conn):
        for sub in range(1, subfolders + 1):
            shp_conn = os.path.join(root.format(project_folder, sub), conn.format(sub, pld))

            # not all 60 day pld files exist
            if not os.path.isfile(shp_conn):
                print("connections line do not exist for " + shp_conn)
                print("creating connection lines")
                # get dest_pts
                dest_pt = os.path.join(root.format(project_folder, sub), pts.format(sub))
                # set pld_int
                pld_int = (int(pld) * 24) / time_step_output
                # get date_start from an existing shapefile. Assume the pld01 will always exist.
                shp_conn_01 = os.path.join(root.format(project_folder, sub), conn.format(sub, '01'))
                c01 = geopandas.read_file(shp_conn_01)
                date_start = c01['date_start'][0]
                connection_lines(dest_pt , seagrass_og, seagrass_crs, shp_conn, date_start, pld_int, int(pld))

            merge_list_conn.append(shp_conn)

        # merge with geopandas
        gdf = geopandas.GeoDataFrame(pd.concat([geopandas.read_file(i) for i in merge_list_conn], ignore_index=True), crs=geopandas.read_file(merge_list_conn[0]).crs)
        gdf.to_file(filename=out_shp_conn, driver='ESRI Shapefile')
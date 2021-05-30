# find the distance between all seagrass points, using land as a barrier

# There is now an updated way to do this, since I last did it with the Cost Path tools:
# https://www.esri.com/arcgis-blog/products/spatial-analyst/analytics/doing-more-with-euclidean-distance-barriers-and-paths/
# https://community.esri.com/t5/arcgis-pro-questions/measure-distance-between-points-across-polygon-least-cost-path/td-p/66607



import arcpy
import os

default_gdb = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_dev_scratch\distance_analysis\distance_analysis_mapping\Default.gdb'
arcpy.env.workspace = default_gdb
seagrass = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_runs_localstage\seagrass\seagrass\seagrass_20200228_SS\seagrass_prep\seagrass.gdb\seagrass_all_19FINAL'
land = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\spatial\boundaries\boundaries.gdb\coastline_bc_ak_wa_or_cleaned_less10000'
clip_fc = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_dev_scratch\distance_analysis\distance_analysis_mapping\Default.gdb\clip_fc'
euc_lines_gdb = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_dev_scratch\distance_analysis\distance_analysis_mapping\euc_lines.gdb'


# get centroids - make sure points fall within polys
arcpy.FeatureToPoint_management(seagrass, 'sg_centroid_fc', 'INSIDE')

# copy land, clip, add attribute, convert to raster
arcpy.CopyFeatures_management(land, 'land_fc')
arcpy.Clip_analysis(land, clip_fc, 'land_fc_clip') # this will cut down on processing time
arcpy.AddField_management('land_fc_clip', 'land', 'SHORT')
arcpy.CalculateField_management('land_fc_clip', 'land', 1)
arcpy.FeatureToRaster_conversion('land_fc_clip', 'land', 'land_r_50', '50')

# check if there are any points that overlap with the land raster
# raster to poly
arcpy.RasterToPolygon_conversion('land_r_50', 'land_r_50_fc', 'NO_SIMPLIFY', 'Value', 'MULTIPLE_OUTER_PART')
# copy centroids
arcpy.CopyFeatures_management('sg_centroid_fc', 'sg_centroid_fc_manualedit')
# !!!MANUALLY!!!: select by location where they intersect, manually move these few points


# rasterize with Feature to Raster tool to make sure I can control cell size and nothing gets lost
# snap to land raster
arcpy.env.snapRaster = 'land_r_50'
arcpy.FeatureToRaster_conversion('sg_centroid_fc_manualedit', 'uID', 'sg_centroid_r', '50')




# go through each point and find distances
# the tools are set up to find all the distances between ONE starting point and multiple destination points
# First, in the Euclidean distance tool, find the distances of every cell back to the one starting point
# Then, in Coast Path as Polyline, input the Destination points to get all the distances from start to dest.

with arcpy.da.SearchCursor('sg_centroid_fc_manualedit', ['uID']) as cursor:
    for row in cursor:
        uID = row[0]
        if uID < 939:
            continue
        print('processing uID {}'.format(uID))
        centroid_sel = arcpy.SelectLayerByAttribute_management('sg_centroid_fc_manualedit', 'NEW_SELECTION', 'uID={}'.format(row[0]))
        arcpy.env.snapRaster = 'land_r_50'
        origin = 'sg_centroid_r_{}'.format(uID)
        arcpy.FeatureToRaster_conversion(centroid_sel, 'uID', origin, '50')

        # Euclidean distance with barriers
        in_rast = origin
        max_dist = ''
        cellsize = 'land_r_50'
        outDirectionRaster = '' # you would use this one if you don't have barriers
        distance_method = 'PLANAR'
        in_barrier_data = 'land_r_50'
        outBackDirectionRaster = 'eucbackdirect' # this direction raster is for when you have barriers

        outEucDistance = arcpy.sa.EucDistance(in_rast, max_dist, cellsize, outDirectionRaster, distance_method, in_barrier_data, outBackDirectionRaster)
        # Save the output 
        outEucDistance.save('eucdistance')

        # Cost Path as Polyline
        inputDestinationLayer = 'sg_centroid_r'
        inputCostLayer = 'eucdistance'
        inputDirectionLayer = 'eucbackdirect'
        outLines = os.path.join(euc_lines_gdb, 'euc_lines_{}'.format(uID))
        pathType = 'EACH_CELL'
        destfield = ''
        arcpy.sa.CostPathAsPolyline(inputDestinationLayer, inputCostLayer, inputDirectionLayer, outLines, pathType, destfield)

        # add origin_uid attribute
        arcpy.AddField_management(outLines, 'origin_id', 'SHORT')
        arcpy.CalculateField_management(outLines, 'origin_id', uID)

        arcpy.Delete_management(origin)
        arcpy.Delete_management('eucbackdirect')
        arcpy.Delete_management('eucdistance')



# merge FCs
arcpy.env.workspace = euc_lines_gdb
fcs = arcpy.ListFeatureClasses()
arcpy.Merge_management(fcs, os.path.join(default_gdb, 'euc_lines_ALL'))



# NOTE
# for some reason the meadow with uID __ could not be processed
# If you look at each line feature class, no distances are calculated to it.
# I tried moving it farther off shore, but it did not make a difference.
# Since all I want is a general trend, I can probably just leave this one out.
# No point wasting time on this.
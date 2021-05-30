# prepare seagrass data
# Basic summary: merge all the seagrass datasets, erase overlap with coastline, aggregate polygons that are within a certain distance, remove holes, delete some polys based on size and distance to other polys
# split into equal area parts

import arcpy
import os

arcpy.env.overwriteOutput = True

# reference individual seagrass feature classes
bcmca = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\1Data_BASE\Seagrass\Seagrass_bcmca\BCMCA_ECO_VascPlants_Eelgrass_Polygons_DATA.shp'
islandstrust = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\1Data_BASE\Seagrass\Seagrass_islandstrust\IslandsTrust_Eelgrass.shp'
misc = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\1Data_BASE\Seagrass\Seagrass_mymapping\misc.gdb\misc'
ross = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\1Data_BASE\Seagrass\Seagrass_mymapping\Ross.gdb\Eelgrass_RossSites'
tsawwassen = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\1Data_BASE\Seagrass\Seagrass_mymapping\tsawwassen.gdb\eelgrass'
hakaiCEC = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\1Data_BASE\Seagrass\Seagrass_polys_Hakai_CEC\BC_SHZN_ZOS_polys.shp'
hakaisurveyed2014 = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\1Data_BASE\Seagrass\Seagrass_polys_Hakai_mapped\seagrass_polys_2014.gdb\seagrass_2014'
hakaisurveyed2016 = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\1Data_BASE\Seagrass\Seagrass_polys_Hakai_mapped\seagrass_polys_2016.gdb\hakai_seagrass_2016'
washington = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\1Data_BASE\Seagrass\Seagrass_Washington\SVMP_distribution.gdb\generalized_eelgrass_poly'
mpatt = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\1Data_BASE\Seagrass\Seagrass_mpatt\mpatt_eco_plants_eelgrass_polygons_data.shp'

# workspace
sg_folder = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_runs\seagrass\seagrass\seagrass_20200228_SS\seagrass_prep'
sg_gdb = 'seagrass.gdb'  # this must exist
arcpy.env.workspace = os.path.join(sg_folder, sg_gdb)
arcpy.env.outputCoordinateSystem = arcpy.SpatialReference("NAD 1983 BC Environment Albers")

# remove meadows where the hydrodynamic model grid is not resolved
remove = True
grid = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\spatial\seagrass\seagrass_prep\seagrass_removeWhereNotResolved.gdb\grid_salishsea'
grid_outline = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\spatial\seagrass\seagrass_prep\seagrass_removeWhereNotResolved.gdb\grid_salishsea_outline'


# buffer Island Trusts data
# these are polygons made from linear features, so it is acceptable to buffer them
# buffer is based on needing these to get somewhat close to the coast so that the overall buffer done later will reach the coast
# first select just presence data (some absence data is included)
islandstrust_presence = arcpy.SelectLayerByAttribute_management(islandstrust, 'NEW_SELECTION', "Presence = 'Y'")
islandstrust_buff = arcpy.Buffer_analysis(islandstrust_presence, "sg_islandstrust_buff", 75 , dissolve_option='ALL')


# merge
sg_merge = arcpy.Merge_management([islandstrust_buff, bcmca, misc, ross, tsawwassen, hakaiCEC, hakaisurveyed2014, hakaisurveyed2016, washington, mpatt], "seagrass_all_01MERGE")


# delete all fields
fieldObjList = arcpy.ListFields(sg_merge)     
fieldNameList = []
for field in fieldObjList:
    if not field.required:
        fieldNameList.append(field.name)
arcpy.DeleteField_management(sg_merge, fieldNameList)


# dissolve
sg_diss = arcpy.Dissolve_management(sg_merge, "seagrass_all_02DISS", multi_part="SINGLE_PART")


# erase areas overlapping the coastline
coast = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\spatial\boundaries\boundaries.gdb\coastline_bc_ak_wa_or_cleaned_less10000'
sg_erase = arcpy.Erase_analysis(sg_diss, coast, "seagrass_all_03ERASE")


# multi to singlepart
sg_multi = arcpy.MultipartToSinglepart_management(sg_erase, "seagrass_all_04MULTI")


# eliminate polygon part
sg_elimpart = arcpy.EliminatePolygonPart_management(sg_multi, "seagrass_all_05ELIM", "PERCENT","", 99, "CONTAINED_ONLY")


# erase with comox polygon
# the issue is here that after all these steps I get a super long polygon that extends for all of Denman Island to Comox. It crosses a shallow (but too deep for seagrass) channel. Talking to Coreen, it is unlikely that there is seagrass here. The polygons are also old and look like theyre modeled from bathymetry.
# http://fishing-app.gpsnauticalcharts.com/i-boating-fishing-web-app/fishing-marine-charts-navigation.html?title=Baynes%20Sound%20boating%20app&fbclid=IwAR2vThmPXTBxRTRg80iuF_OVC8MtXl_n3PG1mkYqnnbKBmW--26OQZWq_Zk#13.14/49.6556/-124.8670
# Would be an interesting free dive!
erase_comox = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\spatial\seagrass\seagrass_prep\seagrass_prep.gdb\comox'
sg_erase_comox = arcpy.Erase_analysis(sg_elimpart, erase_comox, "seagrass_all_05ERASEC")

# aggregate polygons
# this will help cut down on the number of polys. There are many areas where the polygons may have been mapped at high resolution and therefore there are many small polygons in an area. These should be combined. The aggregate tool works surprisingly well.
# use coastline as barrier
sg_agg = arcpy.AggregatePolygons_cartography(sg_erase_comox, "seagrass_all_06AGGREGATE", 100, 0, 0,"NON_ORTHOGONAL", barrier_features=coast)


# do a final check by doing multi-to-single and eliminate holes again
sg_multi_check = arcpy.MultipartToSinglepart_management(sg_agg, "seagrass_all_07MULTI")
sg_elimpart_check = arcpy.EliminatePolygonPart_management(sg_multi_check, "seagrass_all_08ELIM", "PERCENT","", 99, "CONTAINED_ONLY")


# deal with Hornby
# the above steps result in some bigger islands with completely fringing seagrass being filled in
# I used seagrass_all_04MULTI and this python code block: https://gis.stackexchange.com/questions/134066/how-to-find-if-polygon-has-a-hole-using-field-calculator-in-arcgis
# This allows me to investigate which ones get filled in.
# Hornby seems to be the only major issue.
# Therefore, I will use just a hornby featureclass to erase, then a split feature to split it up.
# ONE THING TO NOTE: what happens if I seed particles on land? does it matter? Opendrift correctly moves the points to the ocean, which is great. So now I don't need  to worry about manually breaking polygons apart to account for islands.

hornby = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\spatial\seagrass\seagrass_prep\seagrass_prep.gdb\hornby'
hornby_split = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\spatial\seagrass\seagrass_prep\seagrass_prep.gdb\hornby_split'
sg_erase_hornby = arcpy.Erase_analysis(sg_elimpart_check, hornby, "seagrass_all_09ERASE")
arcpy.Split_analysis(sg_erase_hornby, hornby_split, 'name', arcpy.env.workspace)
# split produce individual features classes for each split
# merge them, erase original feature, merge the split features in
ns_merge = arcpy.Merge_management(['north','south'], 'northsouth_merge')
sg_erase = arcpy.Erase_analysis(sg_erase_hornby, ns_merge, 'seagrass_all_10ERASE')
sg_merge = arcpy.Merge_management(([sg_erase, ns_merge]), 'seagrass_all_11MERGE')


# make copy before near analysis
# this isn't totally necessary, but I want it clear that the NEAR fc will be the last one that includes all polys before deletion occurs
arcpy.FeatureClassToFeatureClass_conversion(sg_merge, arcpy.env.workspace, 'seagrass_all_12NEAR')


# Near analysis - find the distance to the closest feature
sg_near = arcpy.Near_analysis('seagrass_all_12NEAR','seagrass_all_12NEAR')



# create field if it doesn't exist
fields = arcpy.ListFields(sg_near)
fnames = []
for field in fields:
    fnames.append(field.name)
if "delete_patch" not in fnames:
    arcpy.AddField_management(sg_near, 'delete_patch', 'TEXT')


######
# This is where deletion decisions begin.
# Ideally I could go through every meadow, but this is way too time consuming and would be hard to make notes on all my decisions. Therefore, I will script it. There may be patches that fall through the cracks, but I need consistency, so this is the best option.
# I am deleting based on area and distance to closest patch.
# Extremely small meadows are not important if they are adjacent to larger ones. However, if they are isolated, then I still want to keep this meadow.
# However, what about if two meadows are small and next to each other but not next to any others? Prime example is Ucluelet and small group in Deer Group near Bamfield.
######


# if less than 2000m2 and near distance less than 1000m, mark for deletion.
# However, if there are two small meadows isolated then it will delete them both (look at 2 small meadows in deer group north of wizard). Therefore, just delete the smaller one for now.
with arcpy.da.UpdateCursor(sg_near, ['NEAR_FID', 'Shape_Area', 'NEAR_DIST', 'delete_patch', 'OBJECTID']) as cursor:
    for row in cursor:
        if (row[1] < 2000 and row[2] < 1000):
            where = """OBJECTID = {}""".format(row[0])
            with arcpy.da.SearchCursor(sg_near, ['OBJECTID', 'Shape_Area', 'delete_patch'], where) as subcursor:
                for subrow in subcursor:
                    if subrow[1] > row[1]: # delete the smaller one
                        row[3] = 'y'
                        cursor.updateRow(row)


# issue now is that if I just run the next iteration, some of the near features are also ones that were already slated for deletion on the previous step. I should remvove these first then run near again.
# copy to new layer withouth 'y' features
delimitedField = arcpy.AddFieldDelimiters(arcpy.env.workspace, 'delete_patch')
expression = delimitedField + " NOT IN ('y')"
sg_delsmall = arcpy.FeatureClassToFeatureClass_conversion(sg_near, arcpy.env.workspace, 'seagrass_all_13DELSMALL', expression)
# delete near fields
arcpy.DeleteField_management(sg_delsmall, ['ORIG_FID','NEAR_FID','NEAR_DIST'])
# run near again
arcpy.Near_analysis(sg_delsmall, sg_delsmall)


# run the above steps again
# if there are a bunch of little patches clustered together then we would have only deleted half of them
# we can run the same thing again and delete a few more
with arcpy.da.UpdateCursor(sg_delsmall, ['NEAR_FID', 'Shape_Area', 'NEAR_DIST', 'delete_patch', 'OBJECTID']) as cursor:
    for row in cursor:
        if (row[1] < 2000 and row[2] < 1000):
            where = """OBJECTID = {}""".format(row[0])
            with arcpy.da.SearchCursor(sg_delsmall, ['OBJECTID', 'Shape_Area', 'delete_patch'], where) as subcursor:
                for subrow in subcursor:
                    if subrow[1] > row[1]: # delete the smaller one
                        row[3] = 'y'
                        cursor.updateRow(row)
delimitedField = arcpy.AddFieldDelimiters(arcpy.env.workspace, 'delete_patch')
expression = delimitedField + " NOT IN ('y')"
sg_delsmall2 = arcpy.FeatureClassToFeatureClass_conversion(sg_delsmall, arcpy.env.workspace, 'seagrass_all_14DELSMALL2', expression)
arcpy.DeleteField_management(sg_delsmall2, ['ORIG_FID','NEAR_FID','NEAR_DIST'])
arcpy.Near_analysis(sg_delsmall2, sg_delsmall2)


# if less than 30,000m2 and near distance less than 1000m, and if near feature is greater than 30,000m2, then mark for deletion
# so basically, even if it is a decent size, if it is also close to another large one, then get rid of it
with arcpy.da.UpdateCursor(sg_delsmall2, ['NEAR_FID','Shape_Area', 'NEAR_DIST', 'delete_patch']) as cursor:
    for row in cursor:
        if (row[1] < 30000 and row[2] < 1000):
            where = """OBJECTID = {}""".format(row[0])
            with arcpy.da.SearchCursor(sg_delsmall2, ['OBJECTID', 'Shape_Area', 'delete_patch'], where) as subcursor:
                for subrow in subcursor:
                    if subrow[1] > 30000:
                        row[3] = 'y'
                        cursor.updateRow(row)
                    # I'll still delete a patch if one of them is greater than 10,000. In that case, delete the smaller one.
                    elif subrow[1] > 10000 or row[1] > 10000:
                        if subrow[1] > row[1]:
                            row[3] = 'y'
                            cursor.updateRow(row)
            

# copy to new layer withouth 'y' features
delimitedField = arcpy.AddFieldDelimiters(arcpy.env.workspace, 'delete_patch')
expression = delimitedField + " NOT IN ('y')"
sg_delmed = arcpy.FeatureClassToFeatureClass_conversion(sg_delsmall2, arcpy.env.workspace, 'seagrass_all_15DELMED', expression)
# delete near fields
arcpy.DeleteField_management(sg_delmed, ['NEAR_FID','NEAR_DIST'])
# run near again
arcpy.Near_analysis(sg_delmed, sg_delmed)


# the past deletetion should also be ran again as there are still meadows that fit this criteria (imagine a few meadows in a row that are all within 1000m of each other)
# I could probably run this until it know longer deletes meadows. I should have put this into a function, oh well.
with arcpy.da.UpdateCursor(sg_delmed, ['NEAR_FID','Shape_Area', 'NEAR_DIST', 'delete_patch']) as cursor:
    for row in cursor:
        if (row[1] < 30000 and row[2] < 1000):
            where = """OBJECTID = {}""".format(row[0])
            with arcpy.da.SearchCursor(sg_delmed, ['OBJECTID', 'Shape_Area', 'delete_patch'], where) as subcursor:
                for subrow in subcursor:
                    if subrow[1] > 30000:
                        row[3] = 'y'
                        cursor.updateRow(row)
                    elif subrow[1] > 10000 or row[1] > 10000:
                        if subrow[1] > row[1]:
                            row[3] = 'y'
                            cursor.updateRow(row)
delimitedField = arcpy.AddFieldDelimiters(arcpy.env.workspace, 'delete_patch')
expression = delimitedField + " NOT IN ('y')"
sg_delmed2 = arcpy.FeatureClassToFeatureClass_conversion(sg_delmed, arcpy.env.workspace, 'seagrass_all_16DELMED2', expression)
arcpy.DeleteField_management(sg_delmed2, ['NEAR_FID','NEAR_DIST'])
arcpy.Near_analysis(sg_delmed2, sg_delmed2)

# Remove meadows where the hydrodynamic model is not resolved
if remove:
    # clip to a model extent
    sg_clip = arcpy.Clip_analysis(sg_delmed2, grid_outline, 'seagrass_all_17CLIP')
    # select meadows that are within a certain distance to grid
    # For SalishSea: I initially chose 440m since that is the value of one grid cell. However, there were many meadows that would be left out that with diffusion could likely drift into areas where the model is resolved. I then chose 880m. It's never going to be perfect, but I think I can justify that.
    sg_select = arcpy.SelectLayerByLocation_management(sg_clip, 'WITHIN_A_DISTANCE', grid, 880, 'NEW_SELECTION')
    sg_select2 = arcpy.CopyFeatures_management(sg_select, 'seagrass_all_18SELECT')

# delete all unecessary fields
# make a copy first so that I will still have a record of the last geoprocessing step
if not remove:
    sg_final = arcpy.FeatureClassToFeatureClass_conversion(sg_delmed2, arcpy.env.workspace, 'seagrass_all_19FINAL')
else:
    sg_final = arcpy.FeatureClassToFeatureClass_conversion(sg_select2, arcpy.env.workspace, 'seagrass_all_19FINAL')


fieldObjList = arcpy.ListFields(sg_final)     
fieldNameList = []
for field in fieldObjList:
    if not field.required:
        fieldNameList.append(field.name)
arcpy.DeleteField_management(sg_final, fieldNameList)


# add in unique ID, starting from 1
arcpy.AddField_management(sg_final, 'uID', 'SHORT')
arcpy.CalculateField_management(sg_final, 'uID', '!OBJECTID!')
# add in area field that will be used in the biology script
arcpy.AddField_management(sg_final, 'area', 'DOUBLE')
arcpy.CalculateField_management(sg_final, 'area', '!SHAPE.area!')


############################################
# create separate buffer feature class
# This version is used in the biology script just for settlement. It is my workaround for dealing with the many slivers between the seagrass polys and the coastline. Otherwise, when things strand right next to a patch, they wouldn't be considered settled there.
############################################

# buffer seagrass polys by 100 meters. Do not dissolve
# this is based on looking at random patches and seeing how far they are from the coast
sgb_buff = arcpy.Buffer_analysis('seagrass_all_19FINAL', 'seagrass_buff_01BUFF', 100, dissolve_option = 'NONE')

# However, some buffers now overlap. I do not want to dissolve these. But if a particle settles in this overlap then it will likely throw an error. To handle this:
# Do intersect on buffer (just need to add fc once). This gives me just the areas of overlap of the buffers.
sgb_intersect = arcpy.Intersect_analysis(sgb_buff, 'seagrass_buff_02INT')

# then erase from first buffer with this intersect layer
sgb_erase = arcpy.Erase_analysis(sgb_buff, sgb_intersect, 'seagrass_buff_03ERASE')

# However, for areas where a buffer extends into the original area of another seagrass poly, this area has now been erased. We need to add this original portion of the seagrass back in.
# But, it will be missing the buffer. I did aggregate polys by 100m, therefore, the closest polys can now be is ~100.1 meters. With a 100m buffer, this doesn't create much of an issue, so I will leave it as is for now, but if I do less than 100m for aggregate, then I may want to consider doing another  buffer/intersect/erase after this one, but with a smaller buffer, then include this in the merge as well.
# merge
sgb_merge = arcpy.Merge_management(['seagrass_all_19FINAL', sgb_erase], 'seagrass_buff_04MERGE')

# dissolve on uID, uncheck create multipart
sgb_diss = arcpy.Dissolve_management(sgb_merge, "seagrass_buff_05DISS", ['uID'], multi_part="SINGLE_PART")

# Check that the feature count matches the original seagrass feature count. There are likely many slivers created OR if a large buffer is used there are a few places where meadows were right next to each but not connected and therefore their buffers cancel most of each other out.
#  Use FIND IDENTICAL tool on uID.
sgb_identical = arcpy.FindIdentical_management(sgb_diss, 'seagrass_buff_06ID', ['uID'], output_record_option = "ONLY_DUPLICATES")

# OLD WAY (may still want to consider this): using this list, go through and delete the identical record (the sliver portion). Sort by area and deleting slivers will get rid of most of them (these will be <1m2).

# These buffers are useless. Delete the buffer portion entirely. Since it is just four meadows, this is acceptable.
# Delete these pieces from sgb_diss
# Get their uID before deleting
# Select out the uID from the last non-buffered dataset
# Merge this back into the buffered dataset
sgb_dupdelete = arcpy.FeatureClassToFeatureClass_conversion(sgb_diss, arcpy.env.workspace, 'seagrass_buff_07DUPDELETE')
uID_list = []
with arcpy.da.SearchCursor(sgb_identical, ['IN_FID']) as cursor:
    for row in cursor:
        where = """OBJECTID = {}""".format(row[0])
        with arcpy.da.UpdateCursor(sgb_dupdelete, ['OBJECTID','uID'], where) as subcursor:
            for subrow in subcursor:
                uID_list.append(subrow[1])
                subcursor.deleteRow()
uID_list = list(set(uID_list)) # get unique values
sql = "{} IN ({})".format("uID", ", ".join([str(n) for n in uID_list]))
sg_seldup = arcpy.FeatureClassToFeatureClass_conversion(sg_final, arcpy.env.workspace, 'seagrass_buff_07SELECTDUP', sql)
sgb_merge = arcpy.Merge_management([sg_seldup, sgb_dupdelete], 'seagrass_buff_08MERGE')


# run Find Identical again to make sure you got them all
#arcpy.FindIdentical_management(sgb_merge, 'seagrass_buff_99', ['uID'], output_record_option = "ONLY_DUPLICATES")

# With the larger buffer holes will have been created.
# Fill these with Eliminate Polygon part tool. It's ok if islands are filled in for the buffer since they will not seed anything in opendrift.
sgb_elimpart = arcpy.EliminatePolygonPart_management(sgb_merge, "seagrass_buff_09ELIM", "PERCENT","", 99, "CONTAINED_ONLY")

# check for multi part
#arcpy.MultipartToSinglepart_management(sgb_elimpart, "seagrass_buff_999")
# check that the fill didn't create new overlap (intersect)
#arcpy.Intersect_analysis(sgb_elimpart, 'seagrass_buff_9999')


# delete unecessary attributes
sgb_final = arcpy.FeatureClassToFeatureClass_conversion(sgb_elimpart, arcpy.env.workspace, 'seagrass_buff_10FINAL')
fieldObjList = arcpy.ListFields(sgb_final)     
fieldNameList = []
for field in fieldObjList:
    if not field.required and not field.name == "uID":
        fieldNameList.append(field.name)
arcpy.DeleteField_management(sgb_final, fieldNameList)

# check that uID is sequential and without gaps (pretty much that the number matches the feature count)
# spot check that it lines up with the unbuffered version

# you now have your buffer settlement fc

################################################




# output both to shapefile
out_folder = sg_folder
arcpy.FeatureClassToShapefile_conversion(sgb_final, out_folder)
arcpy.FeatureClassToShapefile_conversion(sg_final, out_folder)
arcpy.DeleteField_management(os.path.join(out_folder, 'seagrass_all_19FINAL.shp'), ['Shape_Area', 'Shape_Leng'])
arcpy.DeleteField_management(os.path.join(out_folder, 'seagrass_buff_10FINAL.shp'), ['Shape_Area', 'Shape_Leng'])




#################################################
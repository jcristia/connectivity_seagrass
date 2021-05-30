# split seagrass dataset into multiple parts based on area
# based on intial testing... with the number of particles I will need to seed, it will take way too long to run
# therefore, I will break the patches up and run separate simulations and merge the output shapefiles at the end.
# see notes in Evernote for determing particle number

# area per release
area_per_release = 132771337.1134
total_area = 1194942034.0206974
total_particles = 3798754.26376
num_of_releases = 84
total_split = 9


sg_all = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_runs\seagrass\seagrass\seagrass_20200228_SS\seagrass_prep\seagrass.gdb\seagrass_all_19FINAL'
sg_split_gdb = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_runs\seagrass\seagrass\seagrass_20200228_SS\seagrass_prep\seagrass_split\sg_split.gdb'
sg_split_folder = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_runs\seagrass\seagrass\seagrass_20200228_SS\seagrass_prep\seagrass_split'
# don't need to split buffered dataset since it is only used for settlement in the biology script


import arcpy
import os
arcpy.env.workspace = sg_split_gdb
arcpy.env.overwriteOutput = True


# copy main dataset to split folder
sg = arcpy.FeatureClassToFeatureClass_conversion(sg_all, sg_split_gdb, 'seagrass_all_split')

# add in another area field that will be used for particle count
# Opendrift calculates particles to seed based on areas in WGS84 Auxilary Sphere
# I can keep my shapefiles in BC Albers, but I need areas based on WGS84AS
arcpy.AddField_management(sg, 'areaW84A', 'DOUBLE')
crs = arcpy.SpatialReference(3857)
arcpy.CalculateGeometryAttributes_management(sg, [['areaW84A', 'AREA']], area_unit='SQUARE_METERS', coordinate_system=crs)

# add attribute for split number
arcpy.AddField_management(sg, 'split', 'TEXT')

# with an Update Cursor...
# code each split number starting from 1
# maintain a running total of area
# once I hit area_per_release (round up so that I don't have remainders at the end), add one to the split number and reset total area
split = 1
area_cumu = 0.0
with arcpy.da.UpdateCursor(sg, ['areaW84A', 'split']) as cursor:
    for row in cursor:
        area_cumu += row[0]
        if area_cumu < area_per_release:
            row[1] = 'sg' + str(split)
            cursor.updateRow(row)
        elif split < total_split:
            split += 1
            row[1] = 'sg' + str(split)
            area_cumu = row[0]
            cursor.updateRow(row)
        else:
            row[1] = 'sg' + str(total_split)
            cursor.updateRow(row)

# split based on split number
arcpy.SplitByAttributes_analysis(sg, arcpy.env.workspace, ['split'])

# calc how many particles to seed for each section
fcs = arcpy.ListFeatureClasses()
fcs.remove('seagrass_all_split')
for fc in fcs:
    arcpy.AddField_management(fc, 'particles', 'LONG')
    tot_area_section = 0.0
    with arcpy.da.SearchCursor(fc, 'areaW84A') as cursor:
        for row in cursor:
            tot_area_section += row[0]
    particle = int(((tot_area_section * total_particles) / total_area) / num_of_releases)
    arcpy.CalculateField_management(fc, 'particles', particle)

# check total
# it won't be exactly 32 mil because of the int rounding, but it should be close
fcs = arcpy.ListFeatureClasses()
fcs.remove('seagrass_all_split')
particle_check = 0
for fc in fcs:
    with arcpy.da.SearchCursor(fc, 'particles') as cursor:
        for row in cursor:
            particle_check += row[0]
            break
print(particle_check)

# output as shapefiles
# remove unecessary fields
fcs = arcpy.ListFeatureClasses()
fcs.remove('seagrass_all_split')
for fc in fcs:
    arcpy.FeatureClassToShapefile_conversion(fc, sg_split_folder)
    arcpy.DeleteField_management(os.path.join(sg_split_folder, fc + '.shp'), ['Shape_Area', 'Shape_Leng', 'split'])



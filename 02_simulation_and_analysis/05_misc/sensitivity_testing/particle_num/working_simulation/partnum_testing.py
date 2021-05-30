
# refer to particle number notes in Evernote for details on how I am determining particle 
# this script is used to see how many particles get seeded per meadow
# I am changing basemodel.py so that if calculates 0 particles to seed for small patches, it will add 1 particle and take this difference from the largest meadow

#####################################################

import os
import ogr
import osr
import numpy as np

number = 11306 * 16


sg_path = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_dev_scratch\sensitivity_testing\particle_num\working_simulation\seagrass_og'
sg_files = os.listdir(sg_path)
shapefiles = []
for file in sg_files:
    if file.endswith('.shp'):
        shapefiles.append(os.path.join(sg_path, file))

for shp in shapefiles:
    shp = shp
s = ogr.Open(shp)
for layer in s:
  
    targetSRS = osr.SpatialReference()
    targetSRS.ImportFromEPSG(4326)
    coordTrans = osr.CoordinateTransformation(layer.GetSpatialRef(),
                                                            targetSRS)
    
    layer.ResetReading()
    area_srs = osr.SpatialReference()
    area_srs.ImportFromEPSG(3857)
    areaTransform = osr.CoordinateTransformation(layer.GetSpatialRef(), area_srs)

    featurenum = range(1, layer.GetFeatureCount() + 1)

    areas = np.zeros(len(featurenum))
    for i, f in enumerate(featurenum):
        feature = layer.GetFeature(f - 1)  # Note 1-indexing, not 0
        if feature is not None:
            gom = feature.GetGeometryRef().Clone()
            gom.Transform(areaTransform)
            areas[i] = gom.GetArea()

    total_area = np.sum(areas)
    layer.ResetReading()  # Rewind to first layer
    print('Total area of all polygons: %s m2' % total_area)
    # Find number of points per polygon
    numbers = np.round(number*areas/total_area).astype(int)
    numbers = np.where(numbers == 0, 1, numbers)
    numbers[numbers.argmax()] += np.int(number-sum(numbers))

    for i, f in enumerate(featurenum):
        feature = layer.GetFeature(f - 1)
        if feature is None:
            continue
        num_elements = numbers[i]
        geom = feature.GetGeometryRef()
        print('\tSeeding %s elements within polygon number %s' %
                        (num_elements, featurenum[i]))

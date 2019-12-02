#!/usr/bin/python
###############################################################################
# $Id$
#
# Project:  CARBON SINKS
# Purpose:  set of functions to do resampling, called from resampling_alg.py.
# Author:   Marc Pienaar, marc@saeon.ac.za
#
###############################################################################

from osgeo import ogr, gdal, osr, gdalconst,gdal_array
from math import ceil
import time
import subprocess, glob
import raster_chunking as raster_chunk
import numpy as np
import os, sys
from joblib import Parallel, delayed
import multiprocessing
import time
import subprocess, glob
import gdal_merge as gm

start_time = time.time()
gdal.UseExceptions()
gdal.SetConfigOption("GDAL_DATA", "gdal-data/")

def polygonize(original_raster_binary_1000,original_shp_1000):
	if os.path.exists(original_shp_1000):
		""
	else:
		sourceRaster = gdal.Open(original_raster_binary_1000)	
		band = sourceRaster.GetRasterBand(1)
		driver = ogr.GetDriverByName("ESRI Shapefile")
		outDatasource = driver.CreateDataSource(original_shp_1000)            
		# get proj from raster            
		srs = osr.SpatialReference()
		srs.ImportFromWkt( sourceRaster.GetProjectionRef() )
		# create layer with proj
		outLayer = outDatasource.CreateLayer(original_shp_1000,srs)
		# Add class column (1,2...) to shapefile  
		newField = ogr.FieldDefn('Class', ogr.OFTReal)
		outLayer.CreateField(newField)	
		gdal.Polygonize(band, None,outLayer, 0,[],callback=None)
		outDatasource.Destroy()
		sourceRaster=None	
		band=None	
		ioShpFile = ogr.Open(original_shp_1000, update = 1)
		lyr = ioShpFile.GetLayerByIndex(0)	
		field_defn = ogr.FieldDefn( "Value", ogr.OFTReal )
		lyr.CreateField(field_defn)
		lyr.ResetReading()    
		accum=0
		for i in lyr:
			geom = i.GetGeometryRef()
			area = geom.GetArea()
			lyr.SetFeature(i)
			i.SetField( "Value", str(accum) )
			lyr.SetFeature(i)
			accum=accum+1
		ioShpFile.Destroy()				
		
		print("Polygonization translation took --- %s seconds ---" % (time.time() - start_time))
		return original_shp_1000 


def check_intersect(inputraster,inputgrid,tempshape):
	
	returnval = False
	raster_to_shape(inputraster,tempshape)

	vector2 = ogr.Open(inputgrid)
	layer2 = vector2.GetLayer()
	featureCount2= layer2.GetFeatureCount()
	targetprj = layer2.GetSpatialRef()
	feature2 = layer2.GetFeature(0)
	vectorGeometry2 = feature2.GetGeometryRef()
	ring2 = vectorGeometry2.GetGeometryRef(0)
	
	xmin=ring2.GetPoint(0)[0]
	ymin=ring2.GetPoint(0)[1]
	xmax=ring2.GetPoint(0)[0]
	ymax=ring2.GetPoint(0)[1]
	
	for i in range(0,featureCount2):
		feature2 = layer2.GetFeature(i)
		vectorGeometry2 = feature2.GetGeometryRef()
		ring2 = vectorGeometry2.GetGeometryRef(0)
		points = ring2.GetPointCount()
		for p in range(0,points):
			lon, lat, z = ring2.GetPoint(p)
			xmin=min(xmin,lon)
			xmax=max(xmax,lon)
			ymin=min(ymin,lat)
			ymax=max(ymax,lat)
		
	ring2 = ogr.Geometry(ogr.wkbLinearRing)
	rasterGeometry0 = ogr.Geometry(ogr.wkbPolygon)
	ring2.AddPoint(xmin,ymin)
	ring2.AddPoint(xmax,ymin)
	ring2.AddPoint(xmax,ymax)
	ring2.AddPoint(xmin,ymax)
	ring2.AddPoint(xmin,ymin)
	rasterGeometry0.AddGeometry(ring2)
	
	vector1 = ogr.Open(tempshape)
	layer = vector1.GetLayer()
	sourceprj = layer.GetSpatialRef()
	feature = layer.GetFeature(0)
	vectorGeometry = feature.GetGeometryRef()
	ring = vectorGeometry.GetGeometryRef(0)
	points = ring.GetPointCount()
	
	transform = osr.CoordinateTransformation(sourceprj, targetprj)
	
	ring2 = ogr.Geometry(ogr.wkbLinearRing)
	for p in range(0,points):
		lon, lat, z = ring.GetPoint(p)
		t=transform.TransformPoint(lon, lat)
		ring2.AddPoint(t[0], t[1])
	
	rasterGeometry = ogr.Geometry(ogr.wkbPolygon)
	rasterGeometry.AddGeometry(ring2)
	if(rasterGeometry.Intersect(rasterGeometry0)):
		returnval = True
	else:
		returnval = False

	raster = None#close the resource
	return returnval

def create_temp_shapefile(temp_dir, layer, res_in_meters,sourceprj, i):
	#create the resampling shapefile
	driverName = "ESRI Shapefile"
	drv = ogr.GetDriverByName(driverName)
	feature = layer.GetFeature(i)
	geometry = feature.GetGeometryRef()
	(xmin, xmax, ymin, ymax) = geometry.GetEnvelope()
	xmin = float(xmin)
	xmax = float(xmax)
	ymin = float(ymin)
	ymax = float(ymax)
	gridWidth = float(res_in_meters)
	gridHeight = float(res_in_meters)	
	rows = ceil((ymax-ymin)/gridHeight)
	cols = ceil((xmax-xmin)/gridWidth)	
	# start grid cell envelope
	ringXleftOrigin = xmin
	ringXrightOrigin = xmin + gridWidth
	ringYtopOrigin = ymax
	ringYbottomOrigin = ymax-gridHeight	
	#create a temp shapefile
	outDriver = ogr.GetDriverByName(driverName)
	#create the prj file
	temp_shape_prj = str(temp_dir) + ''.join([str(i), ".prj"])
	sourceprj.MorphToESRI()
	file = open(temp_shape_prj, 'w')
	file.write(sourceprj.ExportToWkt())
	file.close()
	#create the shp file
	temp_shape= str(temp_dir)+''.join([str(i),".shp"])
	outDataSource = outDriver.CreateDataSource(temp_shape)
	outLayer = outDataSource.CreateLayer(temp_dir,geom_type=ogr.wkbPolygon )
	field_name = ogr.FieldDefn("temp", ogr.OFTString)
	field_name.SetWidth(75)
	outLayer.CreateField(field_name)
	countcols = 0
	accum=0
	while countcols < cols:
		countcols += 1
		# reset envelope for rows
		ringYtop = ringYtopOrigin
		ringYbottom =ringYbottomOrigin
		countrows = 0
		while countrows < rows:
			accum+=1
			countrows += 1
			ring = ogr.Geometry(ogr.wkbLinearRing)
			ring.AddPoint(ringXleftOrigin, ringYtop)
			ring.AddPoint(ringXrightOrigin, ringYtop)
			ring.AddPoint(ringXrightOrigin, ringYbottom)
			ring.AddPoint(ringXleftOrigin, ringYbottom)
			ring.AddPoint(ringXleftOrigin, ringYtop)
			poly = ogr.Geometry(ogr.wkbPolygon)
			# reproject the geometry
			# get the input geometry
			geom = feature.GetGeometryRef()				
			poly.AddGeometry(ring)
			# add new geom to layer
			featureDefn = outLayer.GetLayerDefn()
			outFeature = ogr.Feature(featureDefn)
			outFeature.SetGeometry(geom)
			outFeature.SetField("label", ''.join([str(i),"_",str(accum)]))
			# Set new geometry
			outFeature.SetGeometry(poly)
			# Add new feature to output Layer
			outLayer.CreateFeature(outFeature)
			# outLayer.SetProjection(sourceprj)
			outFeature.Destroy
			# new envelope for next poly
			ringYtop = ringYtop - gridHeight
			ringYbottom = ringYbottom - gridHeight
		# new envelope for next poly
		ringXleftOrigin = ringXleftOrigin + gridWidth
		ringXrightOrigin = ringXrightOrigin + gridWidth
	outDataSource.Destroy()
	return temp_shape
	
def raster_to_shape(inraster,outputshape):
	raster = gdal.Open(inraster)
	transform = raster.GetGeoTransform()
	pixelWidth = transform[1]
	pixelHeight = transform[5]
	cols = raster.RasterXSize
	rows = raster.RasterYSize
	xLeft = transform[0]
	yTop = transform[3]
	xRight = xLeft+cols*pixelWidth
	yBottom = yTop+rows*pixelHeight
	ring = ogr.Geometry(ogr.wkbLinearRing)
	ring.AddPoint(xLeft, yTop)
	ring.AddPoint(xLeft, yBottom)
	ring.AddPoint(xRight, yBottom)
	ring.AddPoint(xRight, yTop)
	ring.AddPoint(xLeft, yTop)
	rasterGeometry = ogr.Geometry(ogr.wkbPolygon)
	rasterGeometry.AddGeometry(ring)
	
	sourceprj = osr.SpatialReference(wkt = raster.GetProjection())
	(xmin, xmax, ymin, ymax) = rasterGeometry.GetEnvelope()
	xmin = float(xmin)
	xmax = float(xmax)
	ymin = float(ymin)
	ymax = float(ymax)
	#create the prj file
	temp_shape_prj = str(outputshape) + ''.join([".prj"])#outputGridfn3
	sourceprj.MorphToESRI()
	file = open(temp_shape_prj, 'w')
	file.write(sourceprj.ExportToWkt())
	file.close()
	
	outDriver = ogr.GetDriverByName("Esri Shapefile")
	outDataSource = outDriver.CreateDataSource(outputshape)
	outLayer = outDataSource.CreateLayer('', sourceprj, ogr.wkbPolygon)
	field_name = ogr.FieldDefn("temp", ogr.OFTString)
	field_name.SetWidth(75)
	outLayer.CreateField(field_name)
	poly = ogr.Geometry(ogr.wkbPolygon)
	poly.AddGeometry(ring)
	featureDefn = outLayer.GetLayerDefn()
	outFeature = ogr.Feature(featureDefn)
	outFeature.SetGeometry(rasterGeometry)
	outFeature.SetGeometry(poly)
	outLayer.CreateFeature(outFeature)
	outFeature.Destroy
	outDataSource.Destroy()
	raster=None
	
		
def hirescrop(inputgrid,inputraster,outputraster,temp_raster,temp_raster2,buffer2):	
	#read in the shapfile
	driverName = "ESRI Shapefile"
	drv = ogr.GetDriverByName(driverName)
	dataSource_in = drv.Open(inputgrid, 0) # 0 means read-only. 1 means writeable.
	layer = dataSource_in.GetLayer(0)
	#get target projection
	targetprj = layer.GetSpatialRef()
	minXl, maxXl, minYl, maxYl = layer.GetExtent()
	dataSource_in.Destroy();  # destroy ref to source
	#read in raster and get projection
	ds=gdal.Open(inputraster)
	geotransform = ds.GetGeoTransform()
	pixelWidth = abs(geotransform[1])
	pixelHeight = abs(geotransform[5])
	temp_buffer=max(geotransform[5],abs(geotransform[5]))
	prj=ds.GetProjection()
	sourceprj=osr.SpatialReference(wkt=prj)	
	ds = None
	transform = osr.CoordinateTransformation(targetprj,sourceprj)
	#get extent of target layer in source projection 
	minX,maxY = transform.TransformPoint(minXl, maxYl)[0:2]
	maxX,minY = transform.TransformPoint(maxXl, minYl)[0:2]
	#use gdal.translate to clip the raster + add a small buffer to ensure that we return overlapping areas
	out_ds =gdal.Translate(temp_raster,
			       inputraster,
			       projWin=[minX-temp_buffer-(buffer2), maxY+temp_buffer+(buffer2), maxX+temp_buffer+(buffer2), minY-temp_buffer-(buffer2)])
	out_ds=None
	dataset = None # release resources
	#now apply the new projection
	out_ds=gdal.Warp(temp_raster2,temp_raster,dstSRS=targetprj,xRes = 1.0, yRes = 1.0,multithread =True)
	out_ds=None
	##reclip
	out_ds =gdal.Translate(outputraster,
			       temp_raster2,
			       projWin=[minXl, maxYl, maxXl, minYl])  # geotransform[1], abs(geotransform[5])])
	out_ds=None

def resample(i,inputraster,outputraster,inputgrid,Temp_polygon,mean_dir,sum_dir,temp_dir,path_for_outputs_mean,path_for_outputs_sum,outputGridfn,start,total):
	driverName = "ESRI Shapefile"
	drv = ogr.GetDriverByName(driverName)
	dataSourcein = drv.Open(inputgrid, 0)  # 0 means read-only. 1 means writeable.
	layer = dataSourcein.GetLayer(0)
	featureCount = layer.GetFeatureCount()
	sourceprj2 = layer.GetSpatialRef()
	factor=100/(total-start)
	current_val=(i-start)+1
	percent=current_val*factor
	for j in range(0,layer.GetFeatureCount()):
		print("{} {} {} {} {} {} {} {}".format(i, "of",total,percent,"%", j, " of ", layer.GetFeatureCount()))
		feature2 = layer.GetFeature(j)
		geometry2 = feature2.GetGeometryRef()
		(minX, maxX, minY, maxY) = geometry2.GetEnvelope()
		dataSource_out2 = drv.CreateDataSource(Temp_polygon)
		layer2 = dataSource_out2.CreateLayer(Temp_polygon, sourceprj2)
		layer2.CreateFeature(feature2)
		feature2.Destroy()
		dataSource_out2.Destroy()
		numList = [mean_dir, str(j), ".tif"]
		numList2 = [sum_dir, str(j), ".tif"]
		numList3 = [temp_dir, str(j), ".tif"]		
		tempfilemean = ''.join(numList)
		tempfilesum = ''.join(numList2)
		tempfile = ''.join(numList3)
		if os.path.exists(tempfilemean):
			os.remove(tempfilemean)
		else:
			""
		if os.path.exists(tempfilesum):
			os.remove(tempfilesum)
		else:
			""
		if os.path.exists(tempfile):
			os.remove(tempfile)
		else:
			""
		#clip
		out_ds =gdal.Translate(tempfile,outputraster,projWin=[minX, maxY, maxX, minY]) 
		geotransform = out_ds.GetGeoTransform()
		Projection=out_ds.GetProjection()
		out_ds=None
		rasterArray = gdal_array.LoadFile(tempfile)
		b=rasterArray.flatten()
		ds=gdal.Open(inputraster)
		geotransform = ds.GetGeoTransform()
		pixelWidth = abs(geotransform[1])
		pixelHeight = abs(geotransform[5])
		bandlist = ds.GetRasterBand(1)
		nodata = bandlist.GetNoDataValue()
		ds=None
		#disaggregate
		m=nodata
		try:
			if(len(b[b !=nodata])==0):
				m=nodata
				m2 =nodata
			else:
				m1=(b[b !=nodata]/(pixelWidth*pixelHeight))
				m2 = b[b !=nodata].sum()/len(b[b !=nodata])
				m = m1.sum()
		except:
			m=nodata
			m2=nodata
		source_ds = ogr.Open(Temp_polygon)
		mb_l = source_ds.GetLayer()
		minXl, maxXl, minYl, maxYl = mb_l.GetExtent()
		target_ds = gdal.GetDriverByName('GTiff').Create(tempfilesum, 1, 1, 1, gdal.GDT_Float32) #sum
		target_ds2 = gdal.GetDriverByName('GTiff').Create(tempfilemean, 1, 1, 1, gdal.GDT_Float32) #sum
		target_ds.SetGeoTransform((minX,abs(maxXl-minXl),0,maxY,0,-abs(maxYl-minYl))) 
		target_ds2.SetGeoTransform((minX,abs(maxXl-minXl),0,maxY,0,-abs(maxYl-minYl))) 
		target_ds.SetProjection(Projection)
		target_ds2.SetProjection(Projection)
		bandlist1 = target_ds.GetRasterBand(1)
		bandlist1.SetNoDataValue(nodata)
		bandlist2 = target_ds2.GetRasterBand(1)
		bandlist2.SetNoDataValue(nodata)
		gdal.RasterizeLayer(target_ds, [1], mb_l, None, None,burn_values=[m])
		gdal.RasterizeLayer(target_ds2, [1], mb_l, None, None,burn_values=[m2])
		target_ds = None
		target_ds2 = None
		source_ds = None
	dataSourcein.Destroy();  # source
	mean_merge = [mean_dir, "*.tif"]
	sum_merge = [sum_dir, "*.tif"]		
	meanmerge = ''.join(mean_merge)
	summerge = ''.join(sum_merge)
	files_to_mosaic = glob.glob(meanmerge)
	files_string = ' '.join(files_to_mosaic)
	command = "python gdal_merge.py -a_nodata " + str(nodata)+ " -o "+path_for_outputs_mean +"merge_" +str(i)+ ".tif -of gtiff " + files_string
	os.system(command)
	files_to_mosaic = glob.glob(summerge)
	files_string = ' '.join(files_to_mosaic)
	command = "python gdal_merge.py -a_nodata " + str(nodata)+ " -o "+path_for_outputs_sum +"merge_" +str(i)+ ".tif -of gtiff " + files_string
	os.system(command)
	#delete the temp resources for the next run
	mean_delete = [mean_dir, "*"]
	sum_delete = [sum_dir, "*"]	
	temp_delete = [temp_dir, "*"]
	rest_delete = [outputGridfn, "*"]		
	meanmdelete = ''.join(mean_merge)
	sumdelete= ''.join(sum_merge)
	tempdelete= ''.join(sum_merge)
	restdelete= ''.join(rest_delete)
	
	files = glob.glob(meanmdelete)
	for f in files:
		os.remove(f)
	files = glob.glob(sumdelete)
	for f in files:
		os.remove(f)
	files = glob.glob(tempdelete)
	for f in files:
		os.remove(f)
	files = glob.glob(restdelete)
	for f in files:
		os.remove(f)
	
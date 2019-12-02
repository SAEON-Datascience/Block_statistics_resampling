#!/usr/bin/python
###############################################################################
# $Id$
#
# Project:  CARBON SINKS
# Purpose:  Resampling code example.
# Author:   Marc Pienaar, marc@saeon.ac.za
#
###############################################################################
import os, sys
from osgeo import ogr, gdal, osr,gdal_array
import gdal_merge as gm
from math import ceil
import time
import subprocess, glob
import functions as func

start_time = time.time()
gdal.UseExceptions()
#set the relative path to the gdal-data directory 
gdal.SetConfigOption("GDAL_DATA", "gdal-data/")
#get input data
roi="Data/RSA_BSU_20000.shp"#tiles (20km x 20km resolution) ROI covering South africa
inputraster = "Data/p170r081_TC_2015.tif"# e.g. of the input raster to resample - user to define their own
path_for_outputs = 'output/'#the relative path for the outputs
#both a weighted mean and weighted sum will be created in the output path
path_for_outputs_mean = 'output/mean/'
path_for_outputs_sum = 'output/sum/'
#temp files
temp_dir='outputtemp/'#give a unique name for each chunk
outputraster="outputtemp/temp_output.tif"#give a unique name for each chunk
path_for_outputs_mean_temp='output/mean_temp1/'#give a unique name for each chunk
path_for_outputs_sum_temp='output/sum_temp1/'#give a unique name for each chunk
path_for_outputs_temp='output/temp1/'#give a unique name for each chunk
#range
i1=0
i2=3780#total = 3780 for Sout Africa using the roi file
#resampling size to specify
res_in_meters=1000 #i.e 1km = 1000m output resampling, whereas 500 = 500m resampling nested into the RSA BSU (basic spatial unit)
#the resampling loop - i1 and 12 are the tiles (in 20km x 20km resolution) in the roi - there are 3780 tiles in total, allowoing the data to be processed in chunks
def run(i1,i2,outputraster,temp_dir,path_for_outputs_mean_temp,path_for_outputs_sum_temp,path_for_outputs_temp):
	#create all the directories is they don't exist
	if not os.path.exists(temp_dir):
		os.makedirs(temp_dir)
	if not os.path.exists(path_for_outputs):
		os.makedirs(path_for_outputs)
	if not os.path.exists(path_for_outputs_mean):
		os.makedirs(path_for_outputs_mean)
	if not os.path.exists(path_for_outputs_sum):
		os.makedirs(path_for_outputs_sum)
	if not os.path.exists(path_for_outputs_temp):
		os.makedirs(path_for_outputs_temp)
	if not os.path.exists(path_for_outputs_sum_temp):
		os.makedirs(path_for_outputs_sum_temp)
	if not os.path.exists(path_for_outputs_mean_temp):
		os.makedirs(path_for_outputs_mean_temp)
	#create the resampling shapefile
	driverName = "ESRI Shapefile"
	drv = ogr.GetDriverByName(driverName)
	dataSource_in = drv.Open(roi,0)  # 0 means read-only. 1 means writeable.
	layer = dataSource_in.GetLayer(0)
	featureCount = layer.GetFeatureCount()
	sourceprj = layer.GetSpatialRef()#
	transform = osr.CoordinateTransformation(sourceprj, sourceprj)
	factor=100/(i2-i1)	
	for i in range(i1,i2):
		current_val=(i-i1)+1
		percent=current_val*factor
		print("{} {} {} {} {}".format(i, " of ", i2, percent, "%"))
		#remove contents of folder during each iteration of the loop
		files = glob.glob(temp_dir + '*')
		for f in files:
			os.remove(f)
		#create the temp resmapling shape
		inputgrid = func.create_temp_shapefile(temp_dir, layer, res_in_meters,sourceprj, i)
		#create additional temporary files
		numList1 = [temp_dir, "temp_raster.tif"]
		numList2 = [temp_dir, "temp_raster2.tif"]
		numList3 = [temp_dir, "shpclip.shp"]
		numList4 = [temp_dir, "shp_reproj.shp"]
		temp_raster=''.join(numList1)
		temp_raster2=''.join(numList2)
		Temp_polygon=''.join(numList3)
		tempshape=''.join(numList4)
		#check in input raster intersects the roi, if so perform resmapling
		#comment out the if statement if you want to brute force every chunk in the roi
		if(func.check_intersect(inputraster,inputgrid,tempshape)):
#			print("{} {} {} {} {}".format(i, " of ", layer.GetFeatureCount(), percent, "%"))
			buffer2=res_in_meters*2#for a geographic projection will have to use some other value like 0.001 degrees
			#do a high resolution crop and reprojection ~1m
			func.hirescrop(inputgrid,inputraster,outputraster,temp_raster,temp_raster2,buffer2)
			total =(i2-i1)
			func.resample(i,inputraster,outputraster,inputgrid,Temp_polygon,path_for_outputs_mean_temp,path_for_outputs_sum_temp,path_for_outputs_temp,path_for_outputs_mean,path_for_outputs_sum,temp_dir,i1,i2)
#			return
		else:
			pass	
	dataSource_in.Destroy()  # close the source shapefile
		
run(i1,i2,outputraster,temp_dir,path_for_outputs_mean_temp,path_for_outputs_sum_temp,path_for_outputs_temp)

print("--- %s seconds ---" % (time.time() - start_time))
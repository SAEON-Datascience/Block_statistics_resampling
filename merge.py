#!/usr/bin/python
###############################################################################
# $Id$
#
# Project:  CARBON SINKS
# Purpose:  Merge code to join the outputs from resampling_alg.py
# Author:   Marc Pienaar, marc@saeon.ac.za
#
###############################################################################
import natsort
from osgeo import ogr, gdal, osr, gdalconst, gdal_array

from math import ceil
import time
import subprocess, glob
import numpy as np
import os, sys
import gdal_merge as gm
import functions as func
import raster_chunking as raster_chunk

start_time = time.time()
gdal.UseExceptions()
gdal.SetConfigOption("GDAL_DATA", "gdal-data/")
driver = gdal.GetDriverByName('GTiff')

def merge():
	outmerge="merge/"
	if not os.path.exists(outmerge):
		os.makedirs(outmerge)

#	remove all old files if the exist
	files = glob.glob('merge/*')
	for f in files:
		os.remove(f)
	
	#do mean
	p="_mean"
	linss=['output/mean/*.tif']
	files_to_process = sorted(glob.glob(''.join(linss)))
	files_to_process=natsort.natsorted(files_to_process)
	#get info from a file
	ds = gdal.Open(files_to_process[0])
	x_res = ds.RasterXSize
	y_res = ds.RasterYSize	
	geotransform = ds.GetGeoTransform()
	band = ds.GetRasterBand(1)
	nodata = band.GetNoDataValue()	
	ds=None
	
	linss=['merge_class',str(p),'_']	
	for i in range(0,len(files_to_process),5000):
		files_string = ' '.join(files_to_process[i:i+5000])
		command = "python gdal_merge.py -a_nodata " + str(nodata)+ " -o " + outmerge + ''.join(linss) + str(i)+".tif -of gtiff " + files_string
		os.system(command)
	
	#do sum
	p="_sum"
	linss=['output/mean/*.tif']
	files_to_process = sorted(glob.glob(''.join(linss)))
	files_to_process=natsort.natsorted(files_to_process)	
	linss=['merge_class',str(p),'_']
	for i in range(0,len(files_to_process),5000):
		files_string = ' '.join(files_to_process[i:i+5000])
		command = "python gdal_merge.py -a_nodata " + str(nodata)+ " -o " + outmerge + ''.join(linss) + str(i)+".tif -of gtiff " + files_string
		os.system(command)
merge()


print("--- %s seconds ---" % (time.time() - start_time))
#!/usr/bin/python

###############################################################################
#
# Purpose:  Break larger images into smaller chunks for easier processing
# Code has been adoptede from https://gis.stackexchange.com/questions/92207/split-a-large-geotiff-into-smaller-regions-with-python-and-gdal/92215
#
###############################################################################



import os, sys
import time
import argparse
try:
	from osgeo import ogr, gdal, osr,gdal_array
except:
	import gdal
	import ogr
	import osr
	import gdal_array

global output_path


def get_extent(dataset):
	cols = dataset.RasterXSize
	rows = dataset.RasterYSize
	transform = dataset.GetGeoTransform()
	minx = transform[0]
	maxx = transform[0] + cols * transform[1] + rows * transform[2]
	miny = transform[3] + cols * transform[4] + rows * transform[5]
	maxy = transform[3]

	return {
			"minX": str(minx), "maxX": str(maxx),
			"minY": str(miny), "maxY": str(maxy),
			"cols": str(cols), "rows": str(rows)
			}

def create_tiles(minx, miny, maxx, maxy, n_h,n_w):
	width = maxx - minx
	height = maxy - miny
	matrix = []
	for j in range(n_h, 0, -1):
		for i in range(0, n_w):
			ulx = minx + (width/n_w) * i # 10/5 * 1
			uly = miny + (height/n_h) * j # 10/5 * 1
			lrx = minx + (width/n_w) * (i + 1)
			lry = miny + (height/n_h) * (j - 1)
			matrix.append([[ulx, uly], [lrx, lry]])
	return matrix

def print_factors(x):
	# This function takes a number and prints the factors
	print("The factors of",x,"are:")
	for i in range(1, x + 1):
		 if x % i == 0:
			  print(i)
			
def split(file_name, row_split,col_split,random):	
	
	raw_file_name = os.path.splitext(os.path.basename(file_name))[0].replace("_downsample", "")
	driver = gdal.GetDriverByName('GTiff')
	dataset = gdal.Open(file_name)
	band = dataset.GetRasterBand(1)
	nodata = band.GetNoDataValue()
	transform = dataset.GetGeoTransform()
	extent = get_extent(dataset)
	cols = int(extent["cols"])
	rows = int(extent["rows"])
	print("Columns: ", cols)
	print("Rows: ", rows)
	minx = float(extent["minX"])
	maxx = float(extent["maxX"])
	miny = float(extent["minY"])
	maxy = float(extent["maxY"])

	width = maxx - minx
	height = maxy - miny
	output_path = os.path.join("data", raw_file_name)
	if not os.path.exists(output_path):
		os.makedirs(output_path)
	print("Relative output path: ", output_path)
	#get transform etc
	transform = dataset.GetGeoTransform()
	xOrigin = transform[0]
	yOrigin = transform[3]
	pixelWidth = transform[1]
	pixelHeight = -transform[5]

	n_h=int(rows/row_split)
	n_w=int(cols/col_split)
	#create tile indices
	tiles = create_tiles(minx, miny, maxx, maxy, n_h,n_w)
	tile_num = 0
	
	for tile in tiles:
		minx = tile[0][0]
		maxx = tile[1][0]
		miny = tile[1][1]
		maxy = tile[0][1]

		p1 = (minx, maxy)
		p2 = (maxx, miny)

		i1 = int((p1[0] - xOrigin) / pixelWidth)
		j1 = int((yOrigin - p1[1])  / pixelHeight)
		i2 = int((p2[0] - xOrigin) / pixelWidth)
		j2 = int((yOrigin - p2[1]) / pixelHeight)
		new_cols = i2-i1
		new_rows = j2-j1
		if random:
			import numpy as np
			data = np.random.rand(new_rows,new_cols)
		else:
			data = band.ReadAsArray(i1, j1, new_cols, new_rows)
		new_x = xOrigin + i1*pixelWidth
		new_y = yOrigin - j1*pixelHeight
		print(str(round(i1/pixelWidth)),str(round(j1/pixelHeight)))
		new_transform = (new_x, transform[1], transform[2], new_y, transform[4], transform[5])
		#save tif file as name_row_col
		output_file_base = raw_file_name + "_" + str(round(j1/pixelHeight)) + "_" + str(round(i1/pixelWidth))+ ".tif"
		output_file = os.path.join("data", raw_file_name, output_file_base)
		dst_ds = driver.Create(output_file,
					       new_cols,
					       new_rows,
					       1,
					       gdal.GDT_Float32)
		#writting output raster
		dst_ds.GetRasterBand(1).WriteArray( data )
		tif_metadata = {
			"minX": str(minx), "maxX": str(maxx),
			"minY": str(miny), "maxY": str(maxy)
		}
		dst_ds.SetMetadata(tif_metadata)
		dst_ds.GetRasterBand(1).SetNoDataValue( nodata ) #assuming your raster has 1 band.

		#setting extension of output raster
		# top left x, w-e pixel resolution, rotation, top left y, rotation, n-s pixel resolution
		dst_ds.SetGeoTransform(new_transform)
		wkt = dataset.GetProjection()
		# setting spatial reference of output raster
		srs = osr.SpatialReference()
		srs.ImportFromWkt(wkt)
		dst_ds.SetProjection( srs.ExportToWkt() )
		#Close output raster dataset
		dst_ds = None
		tile_num += 1 
	dataset = None

def getoutput_dir():
	return output_path

#test case
#start_time = time.time()
#gdal.UseExceptions()
#gdal.SetConfigOption("GDAL_DATA", "/Users/privateprivate/SAEON_data/DATA/gdal-data/")
#input_raster = "/Users/privateprivate/Documents/Jean/data/env 2/bathy.tif"
#split(input_raster,10,10)#split into approx 10 by 10 pixels 
#print("--- %s seconds ---" % (time.time() - start_time))

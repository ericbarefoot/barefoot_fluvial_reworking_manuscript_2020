#!/usr/bin/env Rscript
# Aggregate and collapse model-based mapping data into relevant measurements.
# Eric Barefoot
# Sep 2019

# This set of routines takes data as mapped in Agisoft metashape and exported as 
# .shp files, then puts it into a table format as expected for the rest of the 
# analysis for this project (i.e. a tidy data frame with each entry being a measurement)
# this data frame is later appended to the field-based dataset

# first, require libraries

# require()

# now introduce and define necessary functions for collapsing data

# this function calculates the distance along the mapped linear feature.

barfaceDistFun = function(X, Y, Z) 
{
	x = cbind(X,Y,Z)
	
	dists = dist(x)
	
	disc = (1 - 4 * -2 * length(dists))

	if (disc > 0) 
	{
		size = (-1 + sqrt(disc)) / 2
	} 
	else 
	{
		error('no solutions!')
	}

	face = numeric(size)

	j = 0

	for (i in 1:size) 
	{
		face[i] = i + j
		j = j + size - i
	}

	barface_length = sum(dists[face])
	
	return(barface_length)
}

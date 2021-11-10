#EXPORT SpatialPolygon data  TO WINBUGS

# install.packages("maptools")
# install.packages("spdep")
# install.packages("sf")
library(sf)
library(spdep)
library(maptools)



##Export SpatialPolygons object as S-Plus map for WinBUGS.

#sp2WB:The function exports an sp SpatialPolygons object into a S-Plus map format to be import by WinBUGS.
#sp2WB(map, filename, Xscale = 1, Yscale = Xscale, plotorder = FALSE)

setwd("/home/ajo/gitRepos/spatial/lattice/Chapter5")

load("/home/ajo/gitRepos/spatial/lattice/Chapter1/data_nc_sids.Rdata")

nc <-sf::st_read(system.file("shape/nc.shp", package="sf"))
class(nc)
spd <- sf::as_Spatial(st_geometry(nc), IDs = as.character(1:nrow(nc)))

class(spd)
# [1] "SpatialPolygons"
# attr(,"package")
# [1] "sp"

## grab the data from the sf object
df <- nc
df$geometry <- NULL
df <- as.data.frame(df)

## create the SpatialPolygonsDataFrame
nc.sids <- sp::SpatialPolygonsDataFrame(spd, data = df)

#SP2WB: The function exports an sp SpatialPolygons object into a S-Plus map format to be import by WinBUGS.

sp2WB(nc.sids, filename="North_carolina.txt") # Makes a txt-file and saves it in the working directory.

#Results

#The Splus import file is in three parts:

#The first line contains the key word 'map' (lower case) followed by a colon and an integer, N,
#where N is the number of distinct areas in the map (note that one area can consist of more than one polygon).
#The 2nd and 3rd lines are optional, and can be used to specify the units for the map scale.
#By default, GeoBUGS assumes that the polygon coordinates are measured in metres.
#If the coordinates are measured in kilometres, say, then specify Xscale and Yscale to be 1000.
#GeoBUGS will then multiply all polygon co-ordinates by Xscale and Yscale as appropriate before storing the map file.
#If Xscale and Yscale are not specified, then the default units (metres) are assumed.


#The next part of the import file is a 2 column list giving:

#(column 1) the numeric ID of the area - this must be a unique integer between 1 and N;
#the areas should be labelled in the same order as the corresponding data for that area appears in the model.
#(column 2) the area label - this must start with a character, and can be a maximum of 79 alphanumeric characters (no spaces allowed)

#The final part of the import file is a 3 column list giving the co-ordinates of the polygons. The format is:

#(col 1) the label of the area to which the polygon belongs
#(col 2) x-coordinate
#(col 3) y-xoordinate

#The polygon coordinates can be listed either clockwise or anticlockwise. Polygons should be separated by a row of NA's
#The import file should end with the key word:   END

##### Obrir el fitxer exportat en el winbugs i tot seguit Map- Import Splus.
##Open the file txt in WinBUGS and in the menu Map- Import Splus  import de map into WinBUGS

# map:100
# Xscale:1
# Yscale:1
#
# 1 area0
# 2 area1
# 3 area2
# 4 area3
#  ...............
#
# area0 -79.2462 35.8682
# area0 -79.2380 35.8372
# area0 -79.5410 35.8370
# area0 -79.5378 35.8910
# area0 -79.5306 36.2362
# area0 -79.5305 36.2461
# area0 -79.2585 36.2357
# area0 -79.2598 36.0479
# area0 -79.2708 35.9046
# area0 -79.2462 35.8682
# NA NA NA
# area1 -81.1089 35.7719
# area1 -81.1273 35.7890
# area1 -81.1414 35.8233
# area1 -81.3281 35.7951
# ..........................
#
# END


#Import WinBUGS map into R

#ReadSplus Read exported WinBUGS maps:
#The function permits an exported WinBUGS map to be read into an sp package class SpatialPolygons object.

##First, export Winbugs map as Splus (file with txt extension)
##Second: read file in R using ReadSplus
scotland<-readSplus("scotland.txt")
plot(scotland)


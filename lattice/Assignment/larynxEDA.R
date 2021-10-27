setwd("/home/ajo/gitRepos/spatial/lattice/Assignment")

library(sf)
library(spdep)
library(class)

df <- read.delim("larynx_data.txt") # Import Larynx-data. 
head(df)
summary(df)
dim(df)
df$SMR <- df$O/df$E # Add Standardized Morbidity Rate to data frame. 
head(df)

nc <-sf::st_read("NWEngland.shp") # Import geometric data about New England. 
class(nc)
#convert a sf to spatial*Dataframe
spd <- sf::as_Spatial(nc)

# 1) Build a spatial polygon data frame for the England larynx cancer data.
nc.sids <- sp::SpatialPolygonsDataFrame(spd, data = df)
nc.sids <- cbind(nc.sids, nc)
summary(nc.sids)
plot(nc.sids)

# 2) Calculate the adjacency matrix, considering neighbors, those regions that share
# geographic limits and the spatial weights, the standardised to sum unity row.
# a. Identify the regions with higher and lower number or neighbours

xxnb <- poly2nb(nc.sids)
class(xxnb)
plot(nc.sids)#, border="grey")   
plot(xxnb, coordinates(nc.sids), add=TRUE, col="blue")

summary.nb(xxnb) # The summary shows the 10 least connected regions. 
# It also shows that the most connected region is 88 with 16 links. 

# Plot of SMR. 
brks <- round(quantile(nc.sids$SMR, probs=seq(0,1,0.2)), digits=2)   
colours <- c("yellow", "orange2", "red3", "brown", "black")

plot(nc.sids, col=colours[findInterval(nc.sids$SMR, brks,all.inside=TRUE)])
legend(x=c(-84, -80), y=c(33, 34.5), legend=leglabs(brks),fill=colours, bty="n",cex=0.5)
title(main="Standardized Morbidity Rate (SMR) in New England")
text(coordinates(nc.sids),label=as.factor(nc.sids$region_id),cex=0.5)

# Different color Palette

colours<-gray.colors(5,0.95,0.2)
colours<- terrain.colors(5)
colours<- rev(terrain.colors(5))

plot(nc.sids, col=colours[findInterval(nc.sids$SMR, brks,all.inside=TRUE)])
legend(x=c(-84, -80), y=c(33, 34.5), legend=leglabs(brks),fill=colours, bty="n",cex=0.5)
invisible(title(main="Standardized Morbidity Rate (SMR) in New England"))
text(coordinates(nc.sids),label=as.factor(nc.sids$region_id),cex=0.5)


# --------- Choose a neighbourhood criterion. 

#Poly2nb:Construct neighbours list from polygon list based on contiguous boundaries


#Show in the map the regions
plot(nc.sids)#, border="grey")   
label.colors<-c("red",rep("black",99))
label.colors[xxnb[[1]]]<-"green"
text(coordinates(nc.sids),label=as.factor(nc.sids$poly_id),cex=0.5,col=label.colors)




# 3) Calculate the Moran’s I and the Geary’s C indicators and check the independence of the
# process.

# 4) Calculate the local Moran indicators and plot the results. Identify clusters of regions and
# hotspots.

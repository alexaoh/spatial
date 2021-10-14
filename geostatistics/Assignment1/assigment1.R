setwd("/home/ajo/gitRepos/spatial/geostatistics/Assignment1")
library(geoR)
library(sm)
library(gstat)

coordinates <- read.table("poly84.txt", head = TRUE, sep = "\t", dec = ".")
elevs <- read.table("elevationsIslet.txt", head = TRUE, sep = "\t", dec = ".")
head(elevs)

# 1. Exploration of the large scale variability of the elevation.
hist(elevs$data) # Looks ok, close to normal distribution. 
geoelevs <- as.geodata(elevs)
summary(geoelevs)  
plot(geoelevs) # Looks like most of the large values are to the left in the leftmost plot. 
# Also, looks like there is a trend of a second order polynomial in the longitude (X).
# Might also be some sort of rising trend in latitude (Y), from negative to positive values of Y. 

# Remember that: 
###   (circles) : 1st quantile
###   (triangles) : 2nd
###   (plus)  : 3rd
###   (crosses)  : 4th

sm.regression(cbind(elevs$x,elevs$y),elevs$data,display="image",xlab="x",ylab="y")
sm.regression(cbind(elevs$x,elevs$y),elevs$data,display="slice",col=4,add=TRUE)
# It can also be seen here that the large values are at small X-values and large Y-values. 

with(geoelevs, plot(coords[, 1], data, xlab = "W-E",ylab = "Elevation data", pch = 20, cex = 0.7))
lines(with(geoelevs, lowess(data ~ coords[, 1])))
# Here we see that we perhaps could fit a quadratic polynomial to W-E, to remove the trend. 

with(geoelevs , plot(coords[, 2], data, xlab = "S-N",ylab = "Elevation data", pch = 20, cex = 0.7))
lines(with(geoelevs, lowess(data ~ coords[, 2])))
# Here we see that we perhaps could fit a regular linear model to N-S, to remove the trend.

# Remove the trend with a linear regression model. 
lm.fit <- lm(data ~ y + poly(x,2), data = elevs)
summary(lm.fit)

# Add residuals to a new data frame. 
elevs2 <- data.frame(elevs, residuals = lm.fit$residuals)
geoelevs2 <- as.geodata(elevs2, data.col = 4)
plot(geoelevs2)

with(geoelevs2, plot(coords[, 1], data, xlab = "W-E",ylab = "Elevation data", pch = 20, cex = 0.7))
lines(with(geoelevs2, lowess(data ~ coords[, 1])))

with(geoelevs2, plot(coords[, 2], data, xlab = "S-N",ylab = "Elevation data", pch = 20, cex = 0.7))
lines(with(geoelevs2, lowess(data ~ coords[, 2])))
# Looks like some of the trend has been removed. 

# 2. Exploration of the small scale variability of the elevation.

# Calculate and plot variogram cloud and empirical variaogram for the residuals. 
# The robust estimator is used (modulus).
maxdist<-max(dist(cbind(elevs$x,elevs$y))) # Unsure if we should use this to set an upper limit of distances used when calculating variograms. 
cloud <- variog(geoelevs2, option = "cloud", estimator.type = "modulus")
bp <- variog(geoelevs2, option = "bin", bin.cloud = T,
               pairs.min=30, max.dist=maxdist/2,
               estimator.type = "modulus") # Not sure about all these arguments. 
bin <- variog(geoelevs2, option = "bin", pairs.min=25, max.dist=maxdist/2,
                estimator.type = "modulus")

par(mfrow=c(1,3))
plot(cloud,main="CLOUD", cex.main=1, cex.lab=1, cex=2)
plot(bp, bin.cloud=T, cex.lab=1)
title(main = "BINNED BOXPLOTS", cex.main=1)
plot(bin, main="EMPIRICAL VARIOGRAMS (MODULUS)\nBINNED",cex.main=1, cex.lab=1, cex=2, pch=16)

# Will calculate directional variograms to study isotropic/anistropic properties also. 
par(mfrow=c(1, 1))
variod <- variog4(geoelevs2,max.dist=maxdist/2, pairs.min=30,estimator.type = "modulus")
plot(variod,lyt=2,legend=FALSE)
legend(x="bottomright", inset=0.01, lty=c(1,2,3,4), col=c("black", "red", "green","blue"),
       legend=c("0º", "45º", "90º","135º"), cex=0.5)
# I think the data looks relatively isotropic, as all 4 directions checked look relatively similar to each other. 

# 3. Exploration of the spatial independence. Fixed seed set.seed(1000).
set.seed(1000)
par(mfrow=c(1,1))
indep.env<-variog.mc.env(geoelevs2,coords=geoelevs2$coords, data=geoelevs2$data,obj.variog=bp,nsim=200)
plot(bp, envelope = indep.env, main="CONFIDENCE BANDS FOR INDEPENDENT MODEL", lwd=2, pch=16)
# Not sure how to interpret these confidence bands right now. Perhaps I have some notes somewhere explaining it. 

# Add a RoseDiagram also?

# 4. Four theoretical variograms and estimations. 

# 5. Predict elevations along the area of study using the two variaograms selected in 4.

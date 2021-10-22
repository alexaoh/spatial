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

# These plots are perhaps not that informative, could decide later if we want to use them or not in our answer. 
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
# Looks like some om the trend in X has been removed at least. 
# Some of the "rising" trend might also have been slightly removed in Y, even though this was never as clear as X. 

with(geoelevs2, plot(coords[, 1], data, xlab = "W-E",ylab = "Elevation data", pch = 20, cex = 0.7))
lines(with(geoelevs2, lowess(data ~ coords[, 1])))

with(geoelevs2, plot(coords[, 2], data, xlab = "S-N",ylab = "Elevation data", pch = 20, cex = 0.7))
lines(with(geoelevs2, lowess(data ~ coords[, 2])))
# Looks like some of the trend has been removed. 


# 2. Exploration of the small scale variability of the elevation.

# Calculate and plot variogram cloud and empirical variaogram for the residuals. 
# The robust estimator is used (modulus).
maxdist<-max(dist(cbind(elevs$x,elevs$y))) # Practical rules: lags only up to half of maxdist. That is why this is saved here. 
cloud <- variog(geoelevs2, option = "cloud", estimator.type = "modulus")
bp <- variog(geoelevs2, option = "bin", bin.cloud = T,
               pairs.min=30, max.dist=maxdist/2,
               estimator.type = "modulus") # Not sure about all these arguments. 
bin <- variog(geoelevs2, option = "bin", pairs.min=30, max.dist=maxdist/2,
                estimator.type = "modulus")
# Also: pairs,min = 30 since the sample variogram should only be considered for lags that have more than 30 pairs. 

par(mfrow=c(1,3))
plot(cloud,main="CLOUD", cex.main=1, cex.lab=1, cex=2)
plot(bp, bin.cloud=T, cex.lab=1)
title(main = "BINNED BOXPLOTS", cex.main=1)
plot(bin, main="EMPIRICAL VARIOGRAMS (MODULUS)\nBINNED",cex.main=1, cex.lab=1, cex=2, pch=16)
# From these different variaogram plots it looks like the range of a variogram is approx. 30 (distance), 
# the sill around 150 (semivariance) and the nugget might be around 30. 
# Thus, one might say that there is spatial correlation, at least in a range of 30. 

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
# The confidence bands show a range (min to max values) of estimated variaograms when the measurements are randomly permuted in the spatial points. 
# In order for the hypothesis that the process is independent, i.e. that the apparent increase in the estimated variogram 
# (from the original data) might be attributed to chance, the entire variogram should fall inside these confidence bands. 
# It is not clear what to conclude based on this plot only, but the highest point is outside the confidence bands, which might indicate
# that complete spatial randomness in the underlying process may be unplausible! 

# Add a Rose Diagram to study the anisotropy also. 
source("../Chapter3/RoseDiagram.R")
rose.diagram(data.var=geoelevs2$data,data.cds=geoelevs2$coord,max.dist=maxdist/2,numcases=10,numdirec=4,poly.tnd="cte",crit.val=5)

# The idea behind this plot is that one can see if the data looks anisotropic. If this is the case, one can try to apply a anisotropic 
# coordinate correction (rotational transformation) and then make the same plot as above again 
# (with the estimated variogram in four different directions), 
# to see if the data looks more isotropic. This is done in "exploreSmallScaleDep.R", at the bottom of the file. Done in "study_isotropy.R" also. 

# From the Rose Diagram plotted above, it looks like the Anisotry angle is 0 (angle between y-axis and the direction with the maximum range).
# Moreover, it looks like the Anisotropy ration is 2 (ratio between maximum and minimum ranges). 

# Thus we can try the rotational transformation below and observe if it looks more isotropic (even though it looked good to begin with IMO).
angle <- 0 # No rotation in this case, only "squishing".
ratio <- 2   
elevs2.ani <- cbind(coords.aniso(geoelevs2$coords, c(angle,ratio),reverse = FALSE),geoelevs2$data)
geoelevs2.ani<-as.geodata(elevs2.ani)
plot(geoelevs2.ani)
par(mfrow=c(1,1))
variod.ani <- variog4(geoelevs2.ani, max.dist=maxdist/2, pairs.min=30,estimator.type = "modulus");
plot(variod.ani,lyt=2,legend=FALSE,main="Anisotropy angle=0,ratio=2")
legend(x="bottomright", inset=0.01, lty=c(1,2,3,4), col=c("black", "red", "green","blue"),
       legend=c("0\u00B0", "45\u00B0", "90\u00B0","135\u00B0"), cex=0.5)
# These variograms looks about the same as earlier IMO. 
# This tranformation-process only makes sense for geometric anisotropies?

maxdist2<-max(dist(cbind(elevs2.ani[,1],elevs2.ani[,2]))) # Practical rules: lags only up to half of maxdist. That is why this is saved here. 
cloud2 <- variog(geoelevs2.ani, option = "cloud", estimator.type = "modulus")
bp2 <- variog(geoelevs2.ani, option = "bin", bin.cloud = T,
             pairs.min=30, max.dist=maxdist2/2,
             estimator.type = "modulus") # Not sure about all these arguments. 
bin2 <- variog(geoelevs2.ani, option = "bin", pairs.min=30, max.dist=maxdist2/2,
              estimator.type = "modulus")
# Also: pairs.min = 30 since the sample variogram should only be considered for lags that have more than 30 pairs. 

par(mfrow=c(1,3))
plot(cloud2,main="CLOUD", cex.main=1, cex.lab=1, cex=2)
plot(bp2, bin.cloud=T, cex.lab=1)
title(main = "BINNED BOXPLOTS", cex.main=1)
plot(bin2, main="EMPIRICAL VARIOGRAMS (MODULUS)\nBINNED",cex.main=1, cex.lab=1, cex=2, pch=16)

# Plot this sample variogram (after transformation) with the old one.
par(mfrow=c(1,2))
plot(bin, main="EMPIRICAL VARIOGRAMS (MODULUS)\nBINNED Before",cex.main=1, cex.lab=1, cex=2, pch=16)
plot(bin2, main="EMPIRICAL VARIOGRAMS (MODULUS)\nBINNED Transf",cex.main=1, cex.lab=1, cex=2, pch=16)
# Looks like the sills are different. Do we have to conclude on whether this is a geometric, zonal or combined anisotropy?

# 4. Four theoretical variograms and estimations. 
# Will propose four theoretical variograms and estimate them via restricted maximum likelihood (REML). 
# Later will select the two variograms that best fit the data and explain the parameters of the chosen variogram. 

# The ones I propose (at least for now) are the Exponential, Gaussian, Spherical and Matern. 
# This is done with 'bin' for now (which is before the rotation), but not sure if we should use the rotated data instead. 
par(mfrow=c(1, 1))
eyefit(bin, silent = FALSE) # Can be used to find initial values NOT DONE, SHOULD BE DONE BEFORE RUNNING FINAL OPTIMIZATIONS.
# Hence, the initial values used below are just copied from code in chapter 4, should probably be changed.


lk1 <- likfit(geoelevs2, cov.model = "exponential", ini =c(3,0.3), 
              fix.nugget = F, nugget =0.4 ,lik.method = "REML") 

lk2 <- likfit(geoelevs2, cov.model = "gaussian", ini = c(3,0.3),
              fix.nugget = F, nugget =0.4 ,lik.method = "REML")

lk3 <- likfit(geoelevs2, cov.model = "spherical", ini = c(3,0.3),
              fix.nugget = F, nugget =0.4,lik.method = "REML")

lk4 <- likfit(geoelevs2, cov.model = "matern", ini = c(3,0.3),
              fix.nugget = F, nugget =0.4 ,fix.kappa = FALSE, kappa=1,lik.method = "REML") # As we can see, can also estimate kappa (parameter of Matern).

# Increase number of bins to estimate the variograms (as done in 'variogramaEstimationScallops3.R')?
# This estimates a new sample variogram I think. Should this be done earlier, in problem 2/3 already perhaps?

# 5. Predict elevations along the area of study using the two variaograms selected in 4.

# Use kriging to estimate, with the two best variograms selected in problem 4. 

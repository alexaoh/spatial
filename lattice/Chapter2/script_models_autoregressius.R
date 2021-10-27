
# SAR and  CAR models 

##Exemple New York leukemia data

#Data file:  nydata.dbf
#Contains:
# Areaname: name of census tract
# Areakey: unique FIPS code for each tract
# X: x-coordinate of tract centroid (in km)
# UY: y-coordinate of tract centroid (in km)
# Pop8: population size (1980 U.S. Census)
# Tractcas: number of cases 1978-1982 (see text for thorough description).
# Propcas: proportion of cases per tract (Tractcas/Pop8)
# Pctownhome: percentage of people in each tract owning their own home (X3 in regressions).
# Pctage65p: percentage of people in each tract aged 65 or more (X2 in regressions).
# Z: transformed proportions (log(1000*(Tractcas+1)/Pop8)).
# Avgdist: average distance between centroid and  Trichloroethene  sites (TCE).
# Exposure:  "exposure potential" (X1 in regressions).  This is the inverse distance between each census tract centroid and the nearest TCE site, IDIST, transformed via log(100*IDIST) (see pp. 351-353 in text).

#References:
#Waller, L. A. and C. A. Gotway (2004). Applied Spatial Statistics for Public Health Data. John Wiley & Sons, Hoboken, New Jersey. http://web1.sph.emory.edu/users/lwaller/WGindex.htm

library(spdep)

#Definition sptial weights matrix:

library(foreign)
nydata <- read.dbf(system.file("etc/misc/nydata.dbf", package="spdep")[1])
coordinates(nydata) <- c("X", "Y")
plot(nydata)

#spatial weights matrix
nyadjmat <- as.matrix(read.dbf(system.file("etc/misc/nyadjwts.dbf",package="spdep")[1])[-1])
ID <- as.character(names(read.dbf(system.file("etc/misc/nyadjwts.dbf",package="spdep")[1]))[-1])

# Convert a square spatial weights matrix to a weights list object
nyadjlw <- mat2listw(nyadjmat, as.character(nydata$AREAKEY))

#Spatial weights for neighbours lists
listB_NY <- nb2listw(nyadjlw$neighbours, style="B")   #B:Binary
listW_NY <- nb2listw(nyadjlw$neighbours, style="W")   #W:row standardized

#Moran's I test for spatial autocorrelation
ny.moran.B<-moran.test(nydata$Z,listB_NY)
ny.moran.B
# Since the p-values are low, the data are not independent. Thus the null-hypothesis of spatial independence is "forkastet".
# Same conclusion with both types of matrices (weights).

ny.moran.W<-moran.test(nydata$Z,listW_NY)
ny.moran.W
# Since the p-values are low, the data are not independent. Thus the null-hypothesis of spatial independence is "forkastet".
# Same conclusion with both types of matrices (weights).

#Relation between transformed incidence proportions and PEXPOSURE,PCTAGE65P and PCTOWNHOME
nylm<-lm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME,data=nydata)
summary(nylm)
# This model does not take the correlation between the data into account. 


#Moran's I test for residual spatial autocorrelation
lm.morantest(nylm, listB_NY)
lm.morantest(nylm, listW_NY)
# From these tests we can see that the residuals have spatial correlation to a 5% significance level, 
# since the null-hypothesis is "forkastet" under 5% significance (since the p-value is small enough to this significance. 


#Model SARerror or CARerror:
library(spatialreg)

# spautolm {spatialreg}:Spatial conditional and simultaneous autoregression model estimation
#SAR assumptions: Binary   \sigma^2*(I-rho*B)^(-1)(I-rho*B)^(-T)
nysar.B<- spautolm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, listw= listB_NY ,family="SAR",data=nydata)
summary(nysar.B)
summary(nylm)$coefficients

#SAR assumptions: W:row standardized   \sigma^2*(I-rho*W)^(-1)(I-rho*W)^(-T)
nysar.W<- spautolm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, listw= listW_NY ,family="SAR",data=nydata)
summary(nysar.W)


#CAR assumptions Binary   \sigma^2*(I-rho*B)
nycar.B<- spautolm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, listw= listB_NY ,family="CAR",data=nydata)
summary(nycar.B)

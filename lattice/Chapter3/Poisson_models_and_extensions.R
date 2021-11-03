library(sf)
library(foreign)
library(spdep)
library(sp)
library(maptools)
library(MASS)

#EXEMPLE:POISSON REGRESSION

#Lung cancer example from Introductory Statistics with R (ISwR) Book
#E.B. Andersen (1977), Multiplicative Poisson models with unequal cell rates, Scandinavian Journal of Statistics, 4:153-158.
#Data preparation
## Load ISwR package
library(ISwR)
## Load data
data(eba1977)

# Lung cancer incidence in four Danish cities 1968-1971
# 
# Description:
#   This data set contains counts of incident lung cancer cases and
# population size in four neighbouring Danish cities by age group.
# 
# Format:
#   A data frame with 24 observations on the following 4 variables:
#   city a factor with levels Fredericia, Horsens, Kolding, and Vejle.
#   age a factor with levels 40-54, 55-59, 60-64, 65-69, 70-74, and 75+.
#   pop a numeric vector, number of inhabitants.
#   cases a numeric vector, number of lung cancer cases.

head(eba1977)

#Descriptive analysis
eba1977$rate<-eba1977$cases/eba1977$pop
by(eba1977$rate,eba1977$city,summary)
by(eba1977$rate,eba1977$age,summary)


## Fit Poisson model
model.1 <- glm(cases ~ city + age, offset = log(pop), family = poisson(link = "log"), data = eba1977)
summary(model.1)
deviance(model.1)/df.residual(model.1) # Dispersion parameter - check if the Poisson model can be used (or if we have overdispersion).

# The dispersion parameter is 1.5, and is not so large. The Poisson model is probably acceptable. 
# However, the significance of the main effect of living in Kolding city is qualitatively different 
# (significant in Poisson, and non-significant quasi-Poisson).
# Goodness of fit test If the residual deviance is close enough to the residual degrees of freedom, it is a good fit. 
# It can be tested by Chi-squared test.

list(residual.deviance           = deviance(model.1),
     residual.degrees.of.freedom = df.residual(model.1),
     chisq.p.value               = pchisq(deviance(model.1), df.residual(model.1), lower = F)
) # This p-value shows that there is no overdispersion, since the test of dispersion parameter vs degrees of freedom is not significant.
# What hypothesis test is done? Google it or something. The dispersion parameter does not significantly differ (larger) than the 
# degrees of freedom --> this is what the p-value here tells us. 

## Check dispersion parameter with quasi-Poisson regression
model.1q <- glm(cases ~ city + age, offset = log(pop), family = quasipoisson(link = "log"), data = eba1977)
summary(model.1q)
# Can see that there is no overdispersion from the quasi-Poisson regression as well. 

## Results from regular Poisson
summary(model.1)
anova(model.1,test="Chisq")
#Confidence interval for coefficients
confint(model.1)
#Confidence interval for RR
data.frame(RR=round(exp(model.1$coefficients),3), Lower=round(exp(confint(model.1))[,1],3),Upper=round(exp(confint(model.1))[,2],3))

######################################################################
#New York leukemia data taken from the data sets supporting Waller and Gotway 2004

# A data frame with 281 observations on the following 12 variables, and the binary coded spatial weights used in the source.
# AREANAME: name of census tract
# AREAKEY:unique FIPS code for each tract
# X:x-coordinate of tract centroid (in km)
# Y:y-coordinate of tract centroid (in km)
# POP8:population size (1980 U.S. Census)
# TRACTCAS:number of cases 1978-1982
# PROPCAS:proportion of cases per tract
# PCTOWNHOME:percentage of people in each tract owning their own home
# PCTAGE65P:percentage of people in each tract aged 65 or more
# Z:transformed propoprtions
# AVGIDIST:average distance between centroid and TCE sites (TCE:trichloroethylene is a halocarbon commonly used as an industrial solvent.)
# PEXPOSURE:"exposure potential": inverse distance between each census tract centroid and the nearest TCE site, IDIST, transformed via log(100*IDIST)


#############################################################################################################

library(foreign)
library(spdep)
library(sp)
library(maptools)
library(MASS)


#Read dataset 
nydata <- read.dbf(system.file("etc/misc/nydata.dbf", package="spdep")[1])

# Fit Poisson model
head(nydata)
nydata$TRACTCAS<-round(nydata$TRACTCAS,0)

model.poisson.NY <- glm(TRACTCAS ~ PEXPOSURE +PCTAGE65P+ PCTOWNHOME , offset = log(POP8), family = poisson(link = "log"), data = nydata)
summary(model.poisson.NY)
#The dispersion parameter
model.poisson.NY$deviance/model.poisson.NY$df.residual

#Goodness of fit test If the residual deviance is close enough to the residual degrees of freedom, it is a good fit. 
#It can be tested by Chi-squared test.

list(residual.deviance           = deviance(model.poisson.NY),
     residual.degrees.of.freedom = df.residual(model.poisson.NY),
     dispersion.parameter        = deviance(model.poisson.NY)/df.residual(model.poisson.NY),
     chisq.p.value               = pchisq(deviance(model.poisson.NY), df.residual(model.poisson.NY), lower = F)
)

## Quasi-Poisson regression with an estimated dispersion parameter
model.q.poisson.NY<- glm(TRACTCAS ~ PEXPOSURE +PCTAGE65P+ PCTOWNHOME , offset = log(POP8), family = quasipoisson(link = "log"), data = nydata)
summary(model.q.poisson.NY)
summary(model.q.poisson.NY)$coefficients
summary(model.poisson.NY)$coefficients


#Fit Negative binomial regression
library(MASS)
model.nb.NY <- glm.nb(TRACTCAS ~ PEXPOSURE +PCTAGE65P+ PCTOWNHOME+offset(log(POP8)), data = nydata)
summary(model.nb.NY)
model.nb.NY$aic
model.poisson.NY$aic

# Not better fit with negative binomial, since the source of the overdispersion is the spatial correlation in the data
# (which we have seen that this data has earlier). Fitting a negative binomial does not solve this type of overdispersion. 
# Below we check if the overdispersion is due to presence of spatial correlation (as if we did not know this about this specific data from before).

####################################################
#Zero inflated poisson.
library(pscl)

result.zero.inf.NY<-zeroinfl(TRACTCAS ~ PEXPOSURE +PCTAGE65P+ PCTOWNHOME+offset(log(POP8))|1, data = nydata,dist = "poisson", link = "logit")
summary(result.zero.inf.NY)   

#Compute AIC
-2*(result.zero.inf.NY$loglik)+2*5
model.nb.NY$aic
model.poisson.NY$aic


#Check overdispersion is due to presence of spatial correlation

head(nydata)
coordinates(nydata) <- c("X", "Y")

#Read neighours
nyadjmat <- as.matrix(read.dbf(system.file("etc/misc/nyadjwts.dbf",package="spdep")[1])[-1])
ID <- as.character(names(read.dbf(system.file("etc/misc/nyadjwts.dbf",package="spdep")[1]))[-1])
identical(substring(ID, 2, 10), substring(as.character(nydata$AREAKEY), 2, 10))

# Convert a square spatial weights matrix to a weights list object
nyadjlw <- mat2listw(nyadjmat, as.character(nydata$AREAKEY))

#Spatial weights for neighbours lists
#listB_NY <- nb2listw(nyadjlw$neighbours, style="B") #B: Binary
listW_NY <- nb2listw(nyadjlw$neighbours, style="W") #W:row standardized

#Test I Moran for Poisson regression residuals.
moran.mc(model.poisson.NY$residuals,listW_NY,nsim=100)

##########################################################.
#DATA: 100 COUNTIES OF NORTH CAROLINA INCLUDES COUNTS OF NUMBERS OF LIVER BIRTHS
#(ALSO NON-WHITE LIVE BIRTHS) AND NUMBERS OF SUDDEN INFANT DEATHS, FOR 1974-1978 AND 1979-1984.
#THE COUNTY SEAT LOCATION COORDINATES ARE GIVEN IN MILES.


#The nc.sids data frame has 100 rows and 21 columns. It contains data given in Cressie (1991,
#pp. 386-9), Cressie and Read (1985) and Cressie and Chan (1989) on sudden infant deaths in North
#Carolina for 1974-78 and 1979-84. The data set also contains the neighbour list given by Cressie
#and Chan (1989) omitting self-neighbours (ncCC89.nb), and the neighbour list given by Cressie
#and Read (1985) for contiguities (ncCR85.nb). The data are ordered by county ID number, not
#alphabetically as in the source tables sidspolys is a "polylist" object of polygon boundaries, and
#sidscents is a matrix of their centroids.
#Usage

#Format
#This data frame contains the following columns:
#SP_ID SpatialPolygons ID
#CNTY_ID county ID
# east eastings, county seat, miles, local projection
#north northings, county seat, miles, local projection
#L_id Cressie and Read (1985) L index
#M_id Cressie and Read (1985) M index
#names County names
#AREA County polygon areas in degree units
#PERIMETER County polygon perimeters in degree units
#CNTY_ Internal county ID
#NAME County names
#FIPS County ID
#FIPSNO County ID
#CRESS_ID Cressie papers ID
#BIR74 births, 1974-78
#SID74 SID deaths, 1974-78
#NWBIR74 non-white births, 1974-78
#BIR79 births, 1979-84
#SID79 SID deaths, 1979-84
#NWBIR79 non-white births, 1979-84

#Read shapefile into Map object; the file should be given including its ".shp" extension, and the
#function will reconstruct the names of the database (dbf) file and the index (shx) file from these.
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

#Rate
rate<-(1000)*nc.sids$SID74/nc.sids$BIR74
hist(rate)


###Expected cases
r<-sum(nc.sids$SID74)/sum(nc.sids$BIR74)
expected<-r*nc.sids$BIR74
SMR<-nc.sids$SID74/expected
hist(SMR)
nc.sids$expected<-expected
nc.sids$SMR<-SMR

#Fit a Poisson regression
result.pois.SIDS<-glm(SID74~1+offset(log(expected)),data=nc.sids,family="poisson")
summary(result.pois.SIDS)
nc.sids$residuals<-result.pois.SIDS$residuals


list(residual.deviance           = deviance(result.pois.SIDS),
     residual.degrees.of.freedom = df.residual(result.pois.SIDS),
     dispersion.parameter        = deviance(result.pois.SIDS)/df.residual(result.pois.SIDS),
     chisq.p.value               = pchisq(deviance(result.pois.SIDS), df.residual(result.pois.SIDS), lower = F)
)



#The variance is slighty more than twice the mean. 
#The quasibinomial and quasipoisson families differ from the binomial and poisson families only in that the dispersion parameter is not fixed at one, so they can model over-dispersion. 

result.qpois.SIDS <- glm(SID74~1+offset(log(expected)), family = quasipoisson(link = log),data = nc.sids)
summary(result.qpois.SIDS)
summary.glm(result.qpois.SIDS)$dispersion

###############################################
#Fit Binomial Negative.

result.BN.SIDS <- glm.nb(SID74~1+offset(log(expected)), data = nc.sids)
summary(result.BN.SIDS)
result.BN.SIDS$aic
result.pois.SIDS$aic


####################################################
#Fit Zero inflated poisson.
library(pscl)

result.zero.inf.SIDS<-zeroinfl(SID74~1+offset(log(expected)), data = nc.sids,dist = "poisson", link = "logit")
summary(result.zero.inf.SIDS)  

-2*result.zero.inf.SIDS$loglik+2*2
result.pois.SIDS$aic
result.BN.SIDS$aic

##############################.
#Calculate test de moran from residuals of glm model. model.frame()
#########################################3
#Poly2nb:Construct neighbours list from polygon list based on contiguous boundaries

xxnb <- poly2nb(nc.sids)    
plot(nc.sids)#, border="grey")   
plot(xxnb, coordinates(nc.sids), add=TRUE, col="blue")

#nb2listw:The function supplements a neighbours list with spatial weights 
#for the chosen coding scheme.

#B:Binary (1: neighbour, 0:otherwise)
#W: standardised to sum unity row.

w.sids<-nb2listw(xxnb, glist=NULL, style="W",  zero.policy=TRUE)   #Spatial weights for neighbours lists

#Test I Moran for Poisson regression residuals.
moran.mc(result.pois.SIDS$residuals, listw=w.sids, nsim=100)

# This p-value is large --> there is no spatial correlation. She thinks something might be wrong here, see will check it!?

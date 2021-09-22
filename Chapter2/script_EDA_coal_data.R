###########################################
###R script: EXPLORATORY DATA ANALYSIS####
##########################################


### a. Description geostatistical dada. Spaced grid of locations
### b. Description geostatistical dada. Non spaced grid of locations

library(geoR)
library(gstat)

# setwd("E:\\Master d'Estadistica i IO\\Geoestadistica")
# getwd()

## a.Description geostatistical dada. Spaced grid of locations

#The coal seam on the Rodena Mine Property in Greence County, Pennsylvania. 
#The data frame contains 208 coal 
#core samples collectes on a grid given by x and y planar coordinates. 

data(coalash) # in the gstat library
head(coalash)
summary(coalash)

dim(coalash)
plot(coalash[,1],coalash[,2],type="n",xlab="x",ylab="y")
text(coalash[,1],coalash[,2],coalash[,3])

#Converts an Object to the Class "geodata"#
#as.geodata(obj, coords.col = 1:2, data.col = 3,....)

geocoal <- as.geodata(coalash)

names(geocoal)
head(geocoal$coords)
head(geocoal$data)
summary(geocoal)

#Plot of the geodata object 
#Plot of Z_i versus each marginal coordinate (latitude, longitude). 
#Histogram of Z_i

plot(geocoal) 

#Returned:
#2-D Scatterplot:
##Symbols, theirs sizes (and colors) separates data from different quartiles as follows
###   (circles) : 1st quantile
###   (triangles) : 2nd
###   (plus)  : 3rd
###   (crosses)  : 4th


#Alternatives 
#points(geodata) (points.geodata function) generates a plot with the data positions 
#(and by default with the size of the points proportional to the value)

points(geocoal, xlab = "Coord X", ylab = "Coord Y", pt.divide = "rank.proportional")

#pt.divide=c("data.proportional","rank.proportional","quintiles", "quartiles", "deciles", "equal"),

#"data.proportional":sizes proportional to the data values.
#"rank.proportional" :sizes proportional to the rank of the data.
#"quintiles":five different sizes according to the quintiles of the data.
#"quartiles"four different sizes according to the quartiles of the data.
#"deciles"ten different sizes according to the deciles of the data.
#"equal"all points with the same size.

par(mfrow = c(2, 2))
points(geocoal, xlab = "Coord X", ylab = "Coord Y")
points(geocoal, xlab = "Coord X", ylab = "Coord Y", pt.divide = "rank.proportional")
points(geocoal, xlab = "Coord X", ylab = "Coord Y",col = "gray",pt.divide = "equal")
points(geocoal, pt.divide = "quintile", xlab = "Coord X",ylab = "Coord Y")
points(geocoal, pt.divide = "quintile", col = "gray",xlab = "Coord X",ylab = "Coord Y")

#Fix the size of the point.
par(mfrow = c(1, 2))
par(mfrow = c(1, 1))
points(geocoal, cex.min = 2, cex.max = 2, pt.div = "quint")
points(geocoal, cex.min = 2, cex.max = 2, col = "gray")

# Box-plot by columns and row

boxplot(coalash~x , data=coalash, main="Column summaries")
points(sort(unique(coalash[,1])),tapply(coalash[,3],coalash[,1],mean),pch=19,col=2)

boxplot(coalash~y,data=coalash ,main="Row summaries")
points(sort(unique(coalash[,2])),tapply(coalash[,3],coalash[,2],mean),pch=19,col=2)


# we detected an outlier.
# Search outlier:
i.outlier <- which.max(geocoal$data)
i.outlier 
coalash[i.outlier,]
geocoal2 <- as.geodata(coalash[-i.outlier,])
plot(geocoal2)
coalash[i.outlier,]

#Remove outlier
par(mfrow=c(2,2))
points(geocoal2, cex.min = 2, cex.max = 2, col = "gray")
boxplot(coalash~y , data=coalash, subset=-50,horizontal=TRUE, main="Row summaries")
points(tapply(coalash[-50,3],coalash[-50,2],mean),sort(unique(coalash[-50,2])),pch=19,col=2)
boxplot(coalash~x , data=coalash, subset=-50, main="Column summaries")
points(sort(unique(coalash[-50,1])),tapply(coalash[-50,3],coalash[-50,1],mean),pch=19,col=2)

#regres. no-parametrica

library(sm)
par(mfrow=c(1,1))
sm.regression(coalash[-50,1:2],coalash[-50,3],display="image",xlab="x",ylab="y")
sm.regression(coalash[-50,1:2],coalash[-50,3],display="slice",add=TRUE,col=2)


######Linear regression

lm1<-lm(coalash~x+y,data=coalash[-50,])

summary(lm1)

residus<-round(residuals(lm1),digits=3)
coalash_res<-cbind(coalash[-50,],residus) 
head(coalash_res)

geocoal_res <- as.geodata(coalash_res,data.col=4)
plot(geocoal_res)

lm2<-lm(coalash~poly(x,2),data=coalash)
summary(lm2)


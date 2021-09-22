#Exercise 

###Data meuse analyse

#The meuse data set is a data set comprising of four heavy metals measured in
#the top soil in a flood plain along the river Meuse. The governing process seems
#that polluted sediment is carried by the river, and mostly deposited close to the
#river bank.

#In relation to zinc, 
#- Explore the large scale variation.  
#Do you think the mean of the process is constant? Or is it related to the geographical coordinates?  


library(sp) # For the data.
library(geoR) # For as.geodata.
library(sm) # For the contour-plot regression.

data(meuse)
head(meuse)

hist(meuse$zinc)
hist(log(meuse$zinc)) # Bedre med log! Først brukte jeg nederst, men kan gjøres her også.
meuse$logzinc <- log(meuse$zinc) # Legger til log av zinc.


# Col 1:2 er standard, så trenger egentlig ikke å være med!
df <- as.geodata(meuse,data.col=15) # Data is zinc assumed. Change to geodata to be able to make nice plots later. 
summary(df)
head(df)
plot(df)


sm.regression(cbind(meuse$x, meuse$y),meuse$zinc,display="image",xlab="x",ylab="y")
sm.regression(cbind(meuse$x, meuse$y),meuse$zinc,display="slice",col=4,add=TRUE)

# Both the plots seam to show that the level of zinc changes with the coordinates. 
# Looks like the levels of zinc increase in the North-West direction.

# Dette gir samme resultat som plottingen nedenfor!
#plot(df$coords[, 1], df$data, pch = 20, cex = 0.7)
with(df , plot(coords[, 1], data, xlab = "W-E",ylab = "Zinc data", pch = 20, cex = 0.7))
lines(lowess(df$data ~ df$coords[, 1]))
with(df , plot(coords[, 2], data, xlab = "S-N",ylab = "Zinc data", pch = 20, cex = 0.7))
lines(lowess(df$data ~ df$coords[, 2]))

# I will try to rotate the axes, as shown from the example on Scallops. 
rotate.axis <- function(xy,theta=0){
  # xy is a nx2 matrix of coordinates
  # theta is the angle of rotation desireed in degrees
  # theta can be negative
  # XY is the new coordinates
  pimult <- 2*pi*theta/360
  newx <- c(cos(pimult),sin(pimult))
  newy <- c(-sin(pimult),cos(pimult))
  XY <- as.matrix(xy) %*% cbind(newx,newy)
  as.data.frame(XY)
}

xy <- meuse[c("x","y")]
meuse.rot <- cbind(rotate.axis(xy,52),lgzinc=meuse$logzinc)
summary(as.geodata(meuse.rot))
geomeuse.rot<-as.geodata(meuse.rot)
plot(geomeuse.rot)

with(geomeuse.rot , plot(coords[, 1], data, xlab = "newx",ylab = "Zinc data", pch = 20, cex = 0.7))
lines(lowess(geomeuse.rot$data ~ geomeuse.rot$coords[, 1]))
with(geomeuse.rot , plot(coords[, 2], data, xlab = "newy",ylab = "Zinc data", pch = 20, cex = 0.7))
lines(lowess(geomeuse.rot$data ~ geomeuse.rot$coords[, 2]))

# Try to fit a linear regression of the input data. 
lm.meuse.rot <- lm(logzinc ~ poly(newy, 2), data=meuse.rot) # Den beste modellen er ifølge profesora den uten x.
summary(lm.meuse.rot)

#Calculate residuals and analyse
meuse.rot<-data.frame(meuse.rot,residuals=lm.meuse.rot$residuals)
geomeuse.rot.res<-as.geodata(meuse.rot,data.col=4)
plot(geomeuse.rot.res)


###########################################
###R script: EXPLORATORY DATA ANALYSIS####
##########################################
#b. Description geostatistical dada. Non spaced grid of locations


# The data contains the following fields:
#   
# strata: a factor indicating the National Marine Fisheries Service (NMFS) 4-digit strata designator in which the sample was taken.
# sample: sample number per year ranging from 1 to approximately 450.
# lat: latitude location of each sample in the Atlantic Ocean.
# long: longitude location of each sample in the Atlantic Ocean.
# tcatch: total number of scallops caught at the ith sample location. This is prerec + recruits.
# prerec: number of scallops whose shell length is smaller than 70 millimeters.
# recruits: number of scallops whose shell length is 70 millimeters or larger.

setwd("/home/ajo/gitRepos/spatial/Chapter2")
library(geoR)
library(sm)


scallops <- read.table("Scallops_R.txt",head=TRUE,sep=" ",dec=".")
head(scallops)

scallops$long <- -scallops$long

hist(scallops$tcatch)
scallops$lgcatch <- log(1 + scallops$tcatch)
hist(scallops$lgcatch )
geoscalg <- as.geodata(scallops[,-c(1,2)],coords.col = 2:1,data.col=6)
summary(geoscalg)

plot(geoscalg)

par(mfrow=c(1,1))

library(maps)
map("usa")
points(scallops$long,scallops$lat,cex=.5)

par(mfrow = c(1, 1))
map("usa", xlim=c(-75.2,-71),ylim=c(38.5,41.5))
points(scallops$long,scallops$lat,cex=.75)

plot(geoscalg)

sm.regression(cbind(scallops$long,scallops$lat),scallops$lgcatch,display="image",xlab="x",ylab="y")
sm.regression(cbind(scallops$long,scallops$lat),scallops$lgcatch,display="slice",col=4,add=TRUE)

with(geoscalg , plot(coords[, 1], data, xlab = "W-E",ylab = "Scallops data", pch = 20, cex = 0.7))
lines(lowess(geoscalg$data ~ geoscalg$coords[, 1]))
with(geoscalg , plot(coords[, 2], data, xlab = "S-N",ylab = "Scallops data", pch = 20, cex = 0.7))
lines(with(geoscalg, lowess(data ~ coords[, 2])))
theta=52

# rotation axes
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

xy <- scallops[c("long","lat")]
scall.rot <- cbind(rotate.axis(xy,52),lgcatch=scallops$lgcatch)
summary(as.geodata(scall.rot))
geoscall.rot<-as.geodata(scall.rot)
plot(as.geodata(scall.rot))

with(geoscall.rot , plot(coords[, 1], data, xlab = "W-E",ylab = "Scallops data", pch = 20, cex = 0.7))
lines(lowess(geoscall.rot$data ~ geoscall.rot$coords[, 1]))
with(geoscall.rot , plot(coords[, 2], data, xlab = "S-N",ylab = "Scallops data", pch = 20, cex = 0.7))
lines(with(geoscall.rot, lowess(data ~ coords[, 2])))


lm.scp.rot <- lm(lgcatch ~ newx+poly(newy,3), data=scall.rot)
summary(lm.scp.rot )

#Calculate residuals and analyse
scall.rot<-data.frame(scall.rot,residuals=lm.scp.rot$residuals)
geoscall.rot.res<-as.geodata(scall.rot,data.col=4)
plot(geoscall.rot.res)


save.image("results_P1.Rdata")
setwd("/home/ajo/gitRepos/spatial/geostatistics/Chapter3/")

library(geoR)

source("RoseDiagram.R")


#############################################################
# Datos: Scallops.
##############################################################

scallops<- read.table("../Chapter2/Scallops_R.txt",head=TRUE,sep=" ",dec=".")
scallops$lgcatch <- log(scallops$tcatch + 1)
scallops$long <- -scallops$long;
geoscalg <- as.geodata(scallops[,-c(1,2)],coords.col = 2:1,data.col=6)

# rotacion de ejes
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
geoscal.rot<-as.geodata(scall.rot)
plot(geoscal.rot)
lm.scp.rot<-glm(lgcatch ~poly(newx,1)+poly(newy,3),data=scall.rot)
summary(lm.scp.rot)

residus.scallops<-round(residuals(lm.scp.rot),digits=3)
scall.rot.residus <- cbind(scall.rot,residus.scallops)
geoscal.rot.residus<-as.geodata(scall.rot.residus,data.col=4)
plot(geoscal.rot.residus)

#####Variogram empirico
maxdist <- variog(geoscal.rot.residus, option = "cloud")$max.dist;

cloud.2 <- variog(geoscal.rot.residus, option = "cloud", estimator.type = "modulus");
bp.2 <- variog(geoscal.rot.residus, option = "bin", bin.cloud = T,
               pairs.min=30, max.dist=maxdist/2,
               estimator.type = "modulus");
bin.2 <- variog(geoscal.rot.residus, option = "bin", pairs.min=25, max.dist=maxdist/2,
                estimator.type = "modulus");

par(mfrow=c(1,3))

plot(cloud.2,main="CLOUD", cex.main=1, cex.lab=1, cex=2);
plot(bp.2, bin.cloud=T, cex.lab=1);
title(main = "BINNED BOXPLOTS", cex.main=1);
plot(bin.2,main="EMPIRICAL VARIOGRAMS (MODULUS)\nBINNED",cex.main=1, cex.lab=1, cex=2, pch=16);


##############################################################################
#  DETECTION ANISOTROPIA

par(mfrow=c(1,1))
variod <- variog4(geoscal.rot.residus, max.dist=maxdist/2, pairs.min=30,estimator.type = "classical");
plot(variod,lyt=2,legend=FALSE)
legend(x="bottomright", inset=0.01, lty=c(1,2,3,4), col=c("black", "red", "green","blue"),
       legend=c("0\u00B0", "45\u00B0", "90\u00B0","135\u00B0"), cex=0.75);

0.8/0.4
0.8/0.2

rose.diagram(data.var=geoscal.rot.residus$data,data.cds=geoscal.rot.residus$coord,max.dist=maxdist/2,numcases=10,numdirec=6,poly.tnd="cte",crit.val=3)
anisotropy.plot


#Anisotropy angle and anisotrop ratio
#Anisotropy angle
#defined here as the azimuth angle of the direction with greater spatial continuity, 
#i.e. the angle between the y-axis and the direction with the maximum range.

#Anisotropy ratio
#defined here as the ratio between the ranges of the directions with greater and smaller continuity,
#i.e. the ratio between maximum and minimum ranges. 
#Therefore, its value is always greater or equal to one.

angle<-pi/2
ratio<-2   #0.8/0.4
ratio<-3

coords.aniso(geoscal.rot.residus$coords, c(angle,ratio), reverse = FALSE)
scall.ani<-cbind(coords.aniso(geoscal.rot.residus$coords, c(angle,ratio),reverse = FALSE),geoscal.rot.residus$data)
geoscal.ani<-as.geodata(scall.ani)
plot(geoscal.ani)
par(mfrow=c(1,1))
variod.ani <- variog4(geoscal.ani, max.dist=maxdist/2, pairs.min=30,estimator.type = "classical");
plot(variod.ani,lyt=2,legend=FALSE,main="Anisotropy angle=pi/2,ratio=2.9")
legend(x="bottomright", inset=0.01, lty=c(1,2,3,4), col=c("black", "red", "green","blue"),
       legend=c("0\u00B0", "45\u00B0", "90\u00B0","135\u00B0"), cex=0.5)

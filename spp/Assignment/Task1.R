# Nest data "nucli23" and window data "poly23".
library(spatstat)

setwd("/home/ajo/gitRepos/spatial/spp/Task1")
nucli.23=read.table(file="nucli23.txt",header=T,sep="\t")
head(nucli.23)
summary(nucli.23)

poly.23=read.table(file=paste("poly23.txt"),header=T,sep="\t")
summary(poly.23)

# Change reference of X and Y (to 0). 
min.X=min(nucli.23$X)
min.Y=min(nucli.23$Y)
nucli.23$X=nucli.23$X-min.X
nucli.23$Y=nucli.23$Y-min.Y
summary(nucli.23)

poly.23$X=poly.23$X-min.X
poly.23$Y=poly.23$Y-min.Y
summary(poly.23)

# 1) Build a ppp object using the “poly23” data as a window. 

pol.illa<-list(x=poly.23$X,y=poly.23$Y)
n23=ppp(nucli.23$X,nucli.23$Y,poly=pol.illa)
plot(n23,main="nucli23")
axis(1,at=c(seq(0,75,by=25)),pos=c(-10,0))
axis(2,at=c(seq(0,100,by=25)),pos=c(-10,-0))

# 2) Briefly describe the point pattern process. 
# Some EDA follows. 
summary(n23)

n23.den<-density(n23,dimyx=c(256,256))
plot(n23.den)
plot(density(n23, dimyx=c(256,256), sigma=1))
plot(density(n23, dimyx=c(256,256), sigma=5))
plot(density(n23, dimyx=c(256,256), sigma=7.5))
axis(1,at=c(seq(0,75,by=25)),pos=c(-10,0))
axis(2,at=c(seq(0,100,by=25)),pos=c(-10,-0))
plot(density(n23, dimyx=c(256,256), sigma=10))
plot(density(n23, dimyx=c(256,256), sigma=20))
plot(density(n23, dimyx=c(256,256), sigma=50))

# We can see that sigma = 20 and sigma = 50 gives to much smoothing. 
# Similarly, we can see that sigma = 1 gives to little smoothing. 
# The density with sigma between 5 and 10 looks more informative. 

# We can see that the process is not homogeneous. There are some regions that clearly have a larger intensity of points.
# The regions around (x,y) = (60, 25) and (x,y) = (30,60) are two examples, where the first one has the largest intensity. 

# 3) Rebuild the ppp object using the “time to nest” (data_pos variable) as marks.
(n23T=ppp(nucli.23$X,nucli.23$Y,poly=pol.illa,marks=nucli.23$data_pos))

# 4) Briefly describe the marked point process. 
summary(n23T)
plot(n23T,main="nucli23",markscale=0.2,leg.side="right")
axis(1,at=c(seq(0,75,by=25)),pos=c(-10,0))
axis(2,at=c(seq(0,100,by=25)),pos=c(-10,-0))
# Hard to interpret, categorize time into categories for every 4 days instead. 

DPOScat=cut(nucli.23$data_pos,breaks=c(9, 12, 14,16, 18, 20, 23),labels=c("10-12","12-14","14-16","16-18", "18-20", "20-22"))
table(DPOScat) # Are no points 12-14. Thus, merge 10-12 and 12-14 into 10-14.
DPOScat=cut(nucli.23$data_pos,breaks=c(9, 14,16, 18, 20, 22),labels=c("10-14","14-16","16-18", "18-20", "20-22"))
table(DPOScat) # Are no points 12-14. Thus, merge 10-12 and 12-14 into 10-14.
(n23Tcat=ppp(nucli.23$X,nucli.23$Y,poly=pol.illa,marks=DPOScat))
summary(n23Tcat)
plot(split(n23Tcat))

plot(n23Tcat,main="nucli23",cex=0.75,cols=c("purple", "green", "blue","red", "black"),
     chars=rep(16,4),leg.side="right")
axis(1,at=c(seq(0,75,by=25)),pos=c(-10,0))
axis(2,at=c(seq(0,100,by=25)),pos=c(-10,-0))

plot(density(split(n23Tcat), sigma = 7.5), ribbon = F)

nucli.23.df=data.frame(nucli.23,DPOScat)
n23new=split(nucli.23.df,DPOScat)
for (i in 1:5){
  temp.ppp=ppp(n23new[[i]]$X,n23new[[i]]$Y,poly=pol.illa)
  plot(temp.ppp,pch=i,main=c("Group",i))
  axis(1,at=c(seq(0,75,by=25)),pos=c(-10,0))
  axis(2,at=c(seq(0,100,by=25)),pos=c(-10,-0))
  par(ask=T)
}
par(ask=F)

# It is hard to see a clear pattern in how the point process evolves with time, but we can make some remarks: 
# It looks like the area we noted as high intensity at first ((x,y) = (60, 25)) stays high intensity through all times
# that were recorded. One can think that it never gets saturated with nests, i.e. it is a good area for nesting, with room, 
# throughout all recorded time. The other area we noted as relatively high intensity ((x,y) = (30,60)) is most "popular" 
# for nesting between weeks 16 and 20. There is also an area around (x,y) = (10,90) (in the north-west) 
# that has some nests in all weeks. 

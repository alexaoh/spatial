library(spatstat)

# Empty space function
data(cells)
?cells
plot(cells)
class(cells)

Fc <- Fest(cells)
Fc
par(pty = "s")
plot(Fest(cells,correction="best"),legend=F)

# NN function
Gc <- Gest(cells)
Gc
par(pty = "s")
plot(Gest(cells,correction="best"),legend=F)

# Pairwise distances
K <- Kest(cells)
K
par(pty = "s")
plot(Kest(cells,correction="best"),legend=F)


# Envelope

E <- envelope(cells, Kest, nsim = 39,correction="best")
E
plot(E, main = "pointwise envelopes",legend=F)

E <- envelope(cells, Kest, nsim = 19, rank = 1, global = TRUE,correction="best")
E
plot(E, main = "global envelopes",legend=F)

E <- envelope(cells, Lest, nsim = 19, rank = 1, global = TRUE,correction="best")
plot(E, main = "global envelopes of L(r)", legend=F)



#### Fitting models

# Fit cluster point process
# Thomas
data(redwood)
?redwood
plot(redwood)

E <- envelope(redwood, Lest, nsim = 19, rank = 1, global = TRUE,correction="best")
plot(E, main = "global envelopes of L(r)",legend=F)
fit <- kppm(redwood, ~1, "Thomas")
fit



plot(Kest(redwood,correction="best"),legend=F)
plot(fit,legend=F,what="statistic")

#Simulate data
Xsim=rThomas(fit$par[1], sqrt(fit$par[2]), fit$mu,win=redwood$window)
plot(Xsim)

sqrt(fit$par[2])*2

#Matern
fitM <- kppm(redwood, ~1, "MatClust")
fitM

plot(fitM,legend=F,what="statistic")

Xsim=rMatClust(fitM$par[1],fitM$par[2],fitM$mu,win=redwood$window)
plot(Xsim)

qunif(c(0.95),0,0.0865)

# Log-Gaussian Cox process
library(RandomFields)
fitLG <- kppm(redwood, ~1, "LGCP")
fitLG

plot(fitLG,legend=F,what="statistic")

fitLG$par
sLG<-fitLG$par[1]
aLG<-fitLG$par[2]

# Covariance at distance 0.1
sLG*exp(-(0.1/aLG))

#Autocorrelation function
r<-seq(0,1,0.01)
plot(r,exp(-(r/aLG)),type="l")

# Range with autocorrelation = 0.2
aLG*log(0.2)*(-1)


Xsim=rLGCP("exp", mu = fitLG$modelpar[3], var=fitLG$modelpar[1],scale=fitLG$modelpar[2],win=redwood$window)
plot(Xsim)


# Gaussian. It took 32s with i7-8565U CPU @ 1.80GHz 8th Gen  + 16Gb RAM
date()
fitLG2 <- kppm(redwood, ~1, "LGCP",covmodel=list(model="gauss"))
date()
fitLG2

plot(fitLG2,legend=F,what="statistic")
plot(fitLG,add=T,lty=2)

fitLG2$par
sLG2<-fitLG2$par[1]
aLG2<-fitLG2$par[2]

# Covariance at distance 0.1
sLG2*exp(-(0.1/aLG2))

#Autocorrelation function
r<-seq(0,1,0.01)
plot(r,exp(-(r/aLG2)^2),type="l")
lines(r,exp(-(r/aLG)),type="l",lty=2)


# Faster
fitLG2b <- kppm(redwood, ~1, "LGCP",covmodel=list(model="gauss"),method="palm")

plot(fitLG2b,add=T,lty=3)

fitLG2b$par
sLG2b<-fitLG2b$par[1]
aLG2b<-fitLG2b$par[2]

# Covariance at distance 0.1
sLG2b*exp(-(0.1/aLG2b))


#Autocorrelation function
r<-seq(0,1,0.01)
plot(r,exp(-(r/aLG2)),type="l")
lines(r,exp(-(r/aLG)),type="l",lty=2)
lines(r,exp(-(r/aLG2b)),type="l",lty=3)




# The four models together (only available for minimum contrast method)
plot(fit,what="statistic",legend=F)
lines(fitM$Fit$Stat$r,fitM$Fit$mcfit$fit$fit,col="red")
lines(fitLG$Fit$Stat$r,fitLG$Fit$mcfit$fit$fit,col="green")
lines(fitLG2$Fit$Stat$r,fitLG2$Fit$mcfit$fit$fit,col="blue")

# Compare exp vs gauss
fitLG3 <- kppm(redwood, ~1, "LGCP",method="palm")
fitLG3
AIC(fitLG2b)
AIC(fitLG3)


# Tolerance bands
E <- envelope(fit, Lest, nsim = 19, nrank = 1, 
              global = TRUE,correction="best")
plot(E, main = "global envelopes of L",legend=F)

E <- envelope(fitLG, Lest, nsim = 19, nrank = 1, 
              global = TRUE,correction="best")
plot(E, main = "global envelopes of L",legend=F)


##################
# Gibbs process ##
##################

E <- envelope(redwood, Lest, nsim = 19, rank = 1, global = TRUE,correction="best")
plot(E, main = "global envelopes of L(r)", legend=F)


## Cluster process

# Area Interaction
# It takes some time
df=data.frame(r=seq(0.025,0.15,by=0.01))
pfit=profilepl(df,AreaInter,redwood,~1)

pfit
plot(pfit)

fitA<-ppm(redwood, ~1, AreaInter(r=0.085))
summary(fitA)
plot(fitA)

E <- envelope(fitA, Lest, nsim = 19, rank = 1, global = TRUE,correction="best")
plot(E, main = "global envelopes of L(r)",legend=F)


# Geyer

df=data.frame(r=seq(0.05,0.15,by=0.01),sat=2:12)
pfit=profilepl(df,Geyer,redwood,~1)
pfit
plot(pfit)

# Set r=0.1
df=data.frame(r=rep(0.1,11),sat=2:12)
pfit=profilepl(df,Geyer,redwood,~1)
pfit
plot(pfit)

fitG=ppm(redwood, ~1, Geyer(r = 0.1, sat = 7))
fitG
plot(fitG)

fitG2=ppm(redwood, ~1, Geyer(r = 0.1, sat = 4))
fitG2
plot(fitG)
AIC(fitG)
AIC(fitG2)

E <- envelope(fitG2, Lest, nsim = 19, rank = 1, global = TRUE,correction="best")
plot(E, main = "global envelopes of L(r)",legend=F)

AIC(fitA)
AIC(fitG2)


## Regular process

## Hardcore

fitH<-ppm(cells, ~1, Hardcore)
summary(fitH)
plot(fitH)
E <- envelope(fitH, Lest, nsim = 19, rank = 1, global = TRUE,correction="best")
plot(E, main = "global envelopes of L(r)")

## Strauss

df=data.frame(r=seq(0.05,0.15,by=0.01))

pfit=profilepl(df,Strauss,cells,~1)
pfit
plot(pfit)
fitSt<-ppm(cells, ~1, Strauss(r=0.1))
fitSt
plot(fitSt)
E <- envelope(fitSt, Lest, nsim = 19, rank = 1, global = TRUE,correction="best")
plot(E, main = "global envelopes of L(r)")

AIC(fitH)
AIC(fitSt)



########################
# Inhomogeneous Data ##
#######################

nucli.84=read.table(file="nucli84.txt",header=T,sep="\t")
poly.84=read.table(file="poly84.txt",header=T,sep="\t")

#Create ppp object
summary(nucli.84)

min.X=min(nucli.84$X)
min.Y=min(nucli.84$Y)

nucli.84$X=nucli.84$X-min.X
nucli.84$Y=nucli.84$Y-min.Y
summary(nucli.84)

#Create ppp object
# Poligon
poligon=poly.84
summary(poly.84)
poligon$X=poligon$X-min.X
poligon$Y=poligon$Y-min.Y
summary(poligon)


pol.illa<-list(x=poligon$X,y=poligon$Y)
plot(owin(poly=pol.illa),main="nucli84")
n84p=ppp(nucli.84$X,nucli.84$Y,poly=pol.illa)
plot(n84p,main="nucli84")
axis(1,at=c(seq(-25,100,by=25)),pos=c(-75,0))
axis(2,at=c(seq(-75,100,by=25)),pos=c(-25,0))



# Covariables: height and vegetation

# Load data
grid<-read.table("grid.txt",header=T)
grid_veg<-read.table("grid_veg.txt",header=T)

mat<-as.matrix(read.table("height.txt"))
height<-im(mat,grid$x,grid$y)
plot(height,axis=T)
plot(n84p, add=T, cex=0.5)


mat<-as.matrix(read.table("veg.txt"))
veg<-im(mat,grid_veg$x,grid_veg$y)
plot(veg,axis=T)
plot(n84p, add=T, cex=0.5)


# Inhomogeneous PP
fit <- ppm(n84p, ~height + veg)
lam <- predict(fit, locations = n84p)
Ki <- Kinhom(n84p, lam,correction="none")
plot(Ki, main = "Inhomogeneous K function",legend=F)

# Envelope
E <- envelope(fit, Lest, nsim = 39, rank = 1,correction="none")
plot(E, main = "pointwise envelopes")

E <- envelope(fit, Lest, nsim = 19, rank = 1, global = TRUE,correction="none")
plot(E, main = "global envelopes of L(r)")


# Models

#Thomas
fit <- kppm(n84p, ~height + veg, "Thomas")
fit

plot(fit, what="statistic")
plot(fit)

plot(predict(fit))
plot(n84p, add = TRUE, pch = "+")


# Envelope
E <- envelope(fit, Lest, nsim = 39, rank = 1,correction="none")
plot(E, main = "pointwise envelopes")

E <- envelope(fit, Lest, nsim = 19, rank = 1, global = TRUE,correction="none")
plot(E, main = "global envelopes of L(r)")


#Simulate data
Xsim=rThomas(fit$par[1], sqrt(fit$par[2]), fit$mu, win=n84p$window)
plot(Xsim)


# Matern
fitM <- kppm(n84p, ~height + veg, "MatClust")
fitM
plot(fitM,what="statistic")


#Simulate data
Xsim=rMatClust(fitM$par[1],fitM$par[2],fitM$mu,win=n84p$window)
plot(Xsim)

E <- envelope(fitM, Lest, nsim = 19, rank = 1, global = TRUE,correction="none")
plot(E, main = "global envelopes of L(r)")


plot(predict(fitM))
plot(n84p, add = TRUE, pch = "+")


# Log Gaussian Cox
fitLG <- kppm(n84p, ~height+veg, "LGCP")
fitLG

plot(fitLG,what="statistic")
plot(fitLG)

sLG<-fitLG$par[1]
aLG<-fitLG$par[2]

# Covariance at distance 10
sLG*exp(-(6/aLG))

#Autocorrelation function
r<-seq(0,20,1)
plot(r,exp(-(r/aLG)),type="l")

# Range with autocorrelation = 0.2
aLG*log(0.2)*(-1)

# simulation

# Trend function

trend=as.im(log(predict(fitLG)),W=n84p$window)
Xsim=rLGCP("exp", mu =trend, var=fitLG$modelpar[1],scale=fitLG$modelpar[2],win=n84p$window)
plot(Xsim)

E <- envelope(fitLG, Lest, nsim = 19, rank = 1, global = TRUE,correction="none",
              simulate=expression(Xsim=rLGCP("exp", mu = trend, var=fitLG$modelpar[1],scale=fitLG$modelpar[2],win=n84p$window)
))

plot(E, main = "global envelopes of L(r)", legend=F)



# Gibbs process

# Area Interaction
# It takes some time 80s
df=data.frame(r=seq(3,6,by=0.5))
date()
pfit=profilepl(df,AreaInter,n84p,~height+veg,
               correction="none")
date()
pfit
plot(pfit)

fitA<-ppm(n84p, ~height+veg, AreaInter(r=3.5))
summary(fitA)
plot(fitA)

E <- envelope(fitA, Lest, nsim = 19, rank = 1, global = TRUE,correction="none")
plot(E, main = "global envelopes of L(r)")


# Geyer
df=data.frame(r=seq(4,15,by=0.5),sat=1:23)
pfit=profilepl(df,Geyer,n84p,~height+veg,correction="none")
pfit
plot(pfit)

sat<-2:7
r<-5:10

library(tidyr)

df=as.data.frame(crossing(r,sat))
pfit=profilepl(df,Geyer,n84p,~height+veg,correction="none")
pfit
plot(pfit)


fitG=ppm(n84p, ~height+veg, Geyer(r = 6, sat = 7))
fitG
summary(fitG)
plot(fitG)


E <- envelope(fitG, Lest, nsim = 19, rank = 1, global = TRUE,correction="none",
              simulate=expression(Xsim=rLGCP("exp", mu = trend, var=fitLG$modelpar[1],scale=fitLG$modelpar[2],win=n84p$window)
              ))
plot(E, main = "global envelopes of L(r)")


AIC(fitA)
AIC(fitG)


fitG2=ppm(n84p, ~height, Geyer(r = 6, sat = 7))
AIC(fitG2)
summary(fitG2)
plot(fitG2)


#Smoothed residuals
#Smoothed residuals
diagnose.ppm(fitG2, which = "smooth",type="pearson",rbord=0)

# Null standard deviation = 1/(2*bandwidth*sqrt(pi))
0.0158*2


# Lurking variable plot
lurking(fitG2, expression(x), type = "raw",envelope=T)
lurking(fitG2, expression(y), type = "raw",envelope=T)
lurking(fitG2, height, type = "raw",envelope=T)



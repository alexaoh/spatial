library(spatstat)

setwd("/home/ajo/gitRepos/spatial/spp")

# Empty space function
data(cells)
?cells
plot(cells)
class(cells)

# Empty space function. 
Fc <- Fest(cells)
Fc
par(pty = "s")
plot(Fest(cells,correction="best"),legend=T)

# NN (Nearest Neighbor) function
Gc <- Gest(cells)
Gc
par(pty = "s")
plot(Gest(cells,correction="best"),legend=T)

# Pairwise distances
K <- Kest(cells)
K
par(pty = "s")
plot(Kest(cells,correction="best"),legend=T)


# Envelope

E <- envelope(cells, Kest, nsim = 39,correction="best")
E
plot(E, main = "pointwise envelopes",legend=T)

E <- envelope(cells, Kest, nsim = 19, rank = 1, global = TRUE,correction="best")
E
plot(E, main = "global envelopes",legend=T)

# This is the test with more power, as shown in slides (Lest is used instead of Kest).
E <- envelope(cells, Lest, nsim = 19, rank = 1, global = TRUE,correction="best")
plot(E, main = "global envelopes of L(r)", legend=T)



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

# Log-Gaussian Cox process
library(RandomFields)
fitLG <- kppm(redwood, ~1, "LGCP")

plot(fitLG,legend=F,what="statistic")


r<-0.1
1.048*exp(-(r/0.1))

r<-seq(0,1,0.01)
plot(r,exp(-(r/0.1)),type="l")



Xsim=rLGCP("exp", mu = fitLG$modelpar[3], var=fitLG$modelpar[1],scale=fitLG$modelpar[2],win=redwood$window)
plot(Xsim)


# Gaussian. It takes too much time
date()
fitLG2 <- kppm(redwood, ~1, "LGCP",covmodel=list(model="gauss"))
date()
fitLG2

r<-0.1
0.929*exp(-(r/0.108))

r<-seq(0,1,0.01)
plot(r,exp(-(r/0.108)^2),type="l")

# Faster
fitLG2b <- kppm(redwood, ~1, "LGCP",covmodel=list(model="gauss"),method="palm")
fitLG2b



# The three models together (only available for minimum contrast method)
plot(fit,what="statistic",legend=F)
lines(fitM$Fit$Stat$r,fitM$Fit$mcfit$fit$fit)
lines(fitLG$Fit$Stat$r,fitLG$Fit$mcfit$fit$fit)
lines(fitLG2$Fit$Stat$r,fitLG$Fit$mcfit$fit$fit)

# Compare exp vs gauss
fitLG3 <- kppm(redwood, ~1, "LGCP",method="palm")
fitLG3
AIC(fitLG2b)
AIC(fitLG3)


# Tolerance bands
E <- envelope(fit, Lest, nsim = 19, nrank = 1, global = TRUE,correction="best")
plot(E, main = "global envelopes of L",legend=F)

E <- envelope(fitLG2, Lest, nsim = 19, nrank = 1, global = TRUE,correction="best")
plot(E, main = "global envelopes of L")


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


fitG=ppm(redwood, ~1, Geyer(r = 0.1, sat = 7))
fitG
plot(fitG)

E <- envelope(fitG, Lest, nsim = 19, rank = 1, global = TRUE,correction="best")
plot(E, main = "global envelopes of L(r)",legend=F)

AIC(fitA)
AIC(fitG)


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

# Envelope
E <- envelope(n84p, Lest, nsim = 39, rank = 1,correction="none")
plot(E, main = "pointwise envelopes")

E <- envelope(n84p, Lest, nsim = 19, rank = 1, global = TRUE,correction="none")
plot(E, main = "global envelopes of L(r)")

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
Ki
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
summary(fit)
plot(fit, what="statistic")

# Envelope
E <- envelope(fit, Lest, nsim = 39, rank = 1,correction="none")
plot(E, main = "pointwise envelopes")

E <- envelope(fit, Lest, nsim = 19, rank = 1, global = TRUE,correction="none")
plot(E, main = "global envelopes of L(r)")


#Simulate data
Xsim=rThomas(fit$par[1], sqrt(fit$par[2]), fit$mu, win=n84p$window)
plot(Xsim)

plot(predict(fit))
plot(n84p, add = TRUE, pch = "+")



# Matern
fitM <- kppm(n84p, ~height + veg, "MatClust")
fitM
summary(fitM)
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
summary(fitLG)

plot(fitLG,what="statistic")

plot(fitLG)

# simulation

# Trend function

trend=as.im(log(predict(fitLG)),W=n84p$window)
Xsim=rLGCP("exp", mu =trend, var=fitLG$modelpar[1],scale=fitLG$modelpar[2],win=n84p$window)
plot(Xsim)




E <- envelope(fitLG, Lest, nsim = 19, rank = 1, global = TRUE,correction="none",
              simulate=expression(Xsim=rLGCP("exp", mu = trend, var=fitLG$modelpar[1],scale=fitLG$modelpar[2],win=n84p$window)
))
plot(E, main = "global envelopes of L(r)")



# Gibbs process


# Area Interaction
# It takes some time
df=data.frame(r=seq(3,6,by=0.5))
pfit=profilepl(df,AreaInter,n84p,~height+veg,correction="none")
pfit
plot(pfit)

fitA<-ppm(n84p, ~height+veg, AreaInter(r=3.5))
summary(fitA)
plot(fitA)

E <- envelope(fitA, Lest, nsim = 19, rank = 1, global = TRUE,correction="none")
plot(E, main = "global envelopes of L(r)")



# Geyer
df=data.frame(r=seq(0.5,10,by=0.5),sat=1:20)
pfit=profilepl(df,Geyer,n84p,~height+veg,correction="none")
pfit
plot(pfit)
fitG=ppm(n84p, ~height+veg, Geyer(r = 5.5, sat = 11))
fitG
summary(fitG)
plot(fitG)


E <- envelope(fitG, Lest, nsim = 19, rank = 1, global = TRUE,correction="none",
              simulate=expression(Xsim=rLGCP("exp", mu = trend, var=fitLG$modelpar[1],scale=fitLG$modelpar[2],win=n84p$window)
              ))
plot(E, main = "global envelopes of L(r)")


AIC(fitA)
AIC(fitG)

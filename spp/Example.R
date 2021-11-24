library(spatstat)

## Read data
setwd("/home/ajo/gitRepos/spatial/spp")
nucli.84=read.table(file="nucli84.txt",header=T,sep="\t")
summary(nucli.84)
min.X=min(nucli.84$X)
min.Y=min(nucli.84$Y)

nucli.84$X=nucli.84$X-min.X
nucli.84$Y=nucli.84$Y-min.Y
summary(nucli.84)

#Create ppp object
n84=ppp(nucli.84$X,nucli.84$Y,c(-10,150),c(-10,150))
plot(n84)
axis(1,at=c(seq(0,150,by=25)))
axis(2,at=c(seq(0,150,by=25)),pos=c(-10,0))


# Window
# Poligon
poligon=read.table(file=paste("poly84.txt"),header=T,sep="\t")
summary(poligon)
poligon$X=poligon$X-min.X
poligon$Y=poligon$Y-min.Y
summary(poligon)


pol.illa<-list(x=poligon$X,y=poligon$Y)
plot(owin(poly=pol.illa),main="nucli84")
n84p=ppp(nucli.84$X,nucli.84$Y,poly=pol.illa)
plot(n84p,main="nucli84")
axis(1,at=c(seq(0,150,by=25)),pos=c(-100,0))
axis(2,at=c(seq(-100,150,by=25)),pos=c(-10,0))

# point process descriptive analysis

summary(n84p)
n84p.den<-density(n84p)
n84p.den
plot(n84p.den)

n84p.den2<-density(n84p,dimyx=c(256,256))
n84p.den2
plot(n84p.den2)

n84p.den3<-density(n84p,dimyx=c(512,512))
n84p.den3
plot(n84p.den3)

plot(density(n84p, dimyx=c(256,256), sigma=10))

plot(density(n84p, dimyx=c(256,256), sigma=20))

plot(density(n84p, dimyx=c(256,256), sigma=50))

plot(density(n84p, dimyx=c(256,256), sigma=5))

plot(density(n84p, dimyx=c(256,256), sigma=1))
# Want to find a trade-off between smooting and information by changing this sigma!
# Sigma is the standard deviation of isotropic smoothing kernel.
# https://www.rdocumentation.org/packages/spatstat/versions/1.64-1/topics/density.ppp

# marked point process descriptive analysis
# Create marked point pattern
# Time

n84T=ppp(nucli.84$X,nucli.84$Y,poly=pol.illa,marks=nucli.84$DATAPOS)
n84T
summary(n84T)
plot(n84T,main="nucli84")
plot(n84T,main="nucli84",maxsize=10)
plot(n84T,main="nucli84",markscale=0.25)
plot(n84T,main="nucli84",markscale=0.25,leg.side="right")
axis(1,at=c(seq(0,150,by=25)),pos=c(-100,0))
axis(2,at=c(seq(-100,150,by=25)),pos=c(-10,0))

# Time categorized

DPOScat=cut(nucli.84$DATAPOS,breaks=c(0,7,14,21,28),labels=c("W1","W2","W3","W4"))
table(DPOScat)
n84Tcat=ppp(nucli.84$X,nucli.84$Y,poly=pol.illa,marks=DPOScat)
n84Tcat

summary(n84Tcat)

plot(split(n84Tcat))


plot(n84Tcat,main="nucli84",cex=0.75,cols=c("purple", "green", "blue","red"),
     chars=rep(16,4),leg.side="right")
axis(1,at=c(seq(0,150,by=25)),pos=c(-100,0))
axis(2,at=c(seq(-100,150,by=25)),pos=c(-10,0))

plot(n84Tcat,main="nucli84",cex=0.75,cols=c("purple", "green", "blue","red"),
     chars=rep(16,4),legend=F)
axis(1,at=c(seq(0,150,by=25)),pos=c(-100,0))
axis(2,at=c(seq(-100,150,by=25)),pos=c(-10,0))
legend(110,0,levels(DPOScat),col=c("purple", "green","blue","red"),pch=16)


# Alternative to watch the plot
nucli.84=data.frame(nucli.84,DPOScat)
n84new=split(nucli.84,DPOScat)
for (i in 1:4){
  temp.ppp=ppp(n84new[[i]]$X,n84new[[i]]$Y,poly=pol.illa)
  plot(temp.ppp,pch=i,main=c("Week",i))
  axis(1,at=c(seq(0,150,by=25)),pos=c(-100,0))
  axis(2,at=c(seq(-100,150,by=25)),pos=c(-10,0))
  par(ask=T)
}

par(ask=F)



plot(density(split(n84Tcat)), ribbon = F)


plot(density(split(n84Tcat),sigma=10), ribbon = FALSE)

# Relative proportions of intensity

Y <- density(split(n84Tcat),sigma=20)
attach(Y)


pW1 <- eval.im(W1/(W1 + W2 + W3 + W4))
plot(pW1)
pW2 <- eval.im(W2/(W1 + W2 + W3 + W4))
plot(pW2)
pW3 <- eval.im(W3/(W1 + W2 + W3 + W4))
plot(pW3)
pW4 <- eval.im(W4/(W1 + W2 + W3 + W4))
plot(pW4)
detach(Y)


data(lansing)
plot(lansing)

plot(split(lansing))


Y <- density(split(lansing))
attach(Y)

W1<-Y$blackoak
W2<-Y$hickory
W3<-Y$maple
W4<-Y$misc
W5<-Y$redoak
W6<-Y$whiteoak

pW1 <- eval.im(W1/(W1 + W2 + W3 + W4+W5+W6))
plot(pW1)
pW2 <- eval.im(W2/(W1 + W2 + W3 + W4+W5+W6))
plot(pW2)
pW3 <- eval.im(W3/(W1 + W2 + W3 + W4+W5+W6))
plot(pW3)
pW4 <- eval.im(W4/(W1 + W2 + W3 + W4+W5+W6))
plot(pW4)
pW5 <- eval.im(W5/(W1 + W2 + W3 + W4+W5+W6))
plot(pW5)
pW6 <- eval.im(W6/(W1 + W2 + W3 + W4+W5+W6))
plot(pW6)
detach(Y)

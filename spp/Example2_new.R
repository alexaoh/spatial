library(spatstat)
data(bei)
plot(bei)


#Homogeneous
summary(bei)
lamb=summary(bei)$intensity
lamb

quadratcount(bei, nx = 4, ny = 2)
Q <- quadratcount(bei, nx = 4, ny = 2)
plot(bei, cex = 0.5, pch = "+")

plot(Q, add = TRUE, cex = 2)


#nonparametric estimation of density

# Searching bandwidth
b1<-bw.diggle(bei)

b2<-bw.ppl(bei)

b3<-bw.scott(bei)

b4<-bw.CvL(bei)



b1
b2
b3
b4
plot(b1)
plot(b2)
plot(b4)


# Gaussian kernel

den <- density(bei, sigma = b1)
plot(den)
plot(bei, add = TRUE, cex = 0.5)

den <- density(bei, sigma = b2)
plot(den)
plot(bei, add = TRUE, cex = 0.5)

den <- density(bei, sigma = b3[1])
plot(den)
plot(bei, add = TRUE, cex = 0.5)

den <- density(bei, sigma = b3[2])
plot(den)
plot(bei, add = TRUE, cex = 0.5)

den <- density(bei, sigma = b4)
plot(den)
plot(bei, add = TRUE, cex = 0.5)


den <- density(bei, sigma = 50)
plot(den)
plot(bei, add = TRUE, cex = 0.5)


# Exploratory analysis
# quadrats determined by a covariate
# Gradient
Z <- bei.extra$grad
b <- quantile(Z, probs = (0:4)/4)
b
Zcut <- cut(Z, breaks = b, labels = 1:4)
V <- tess(image = Zcut)
plot(V,valuesAreColours=FALSE)
plot(bei, add = TRUE, pch = 16)

qb <- quadratcount(bei, tess = V)

plot(qb,valuesAreColours=FALSE)

# elevation
Z <- bei.extra$elev
b <- quantile(Z, probs = (0:4)/4)
b
Zcut <- cut(Z, breaks = b, labels = 1:4)
V <- tess(image = Zcut)
plot(V,valuesAreColours=FALSE)
plot(bei, add = TRUE, pch = "+")
qb <- quadratcount(bei, tess = V)

plot(qb,valuesAreColours=FALSE)


# Nucli 84 data

nucli.84=read.table(file="nucli84.txt",header=T,sep="\t")
min.X=min(nucli.84$X)
min.Y=min(nucli.84$Y)
nucli.84$X=nucli.84$X-min.X
nucli.84$Y=nucli.84$Y-min.Y


poligon=read.table(file=paste("poly84.txt"),header=T,sep="\t")
poligon$X=poligon$X-min.X
poligon$Y=poligon$Y-min.Y

pol.illa<-list(x=poligon$X,y=poligon$Y)
n84p=ppp(nucli.84$X,nucli.84$Y,poly=pol.illa)
plot(n84p,main="nucli84")
axis(1,at=c(seq(0,150,by=25)),pos=c(-100,0))
axis(2,at=c(seq(-100,150,by=25)),pos=c(-10,0))


summary(n84p)

quadratcount(n84p, nx = 2, ny = 4)

Q <- quadratcount(n84p, nx = 2, ny = 4)
plot(n84p, cex = 0.5, pch = "+")

plot(Q, add = TRUE, cex = 2)
axis(1,at=c(seq(0,150,by=25)),pos=c(-100,0))
axis(2,at=c(seq(-100,150,by=25)),pos=c(-10,0))

# An other option
quadratcount(n84p, nx = 2, ybreaks = c(-100,0,25,50,100))

Q<-quadratcount(n84p, nx = 2, ybreaks = c(-100,0,25,50,100))
plot(n84p, cex = 0.5, pch = "+")

plot(Q, add = TRUE, cex = 2)


#nonparametric estimation of density

# Searching bandwidth
b1<-bw.diggle(n84p,correction="none")
b2<-bw.ppl(n84p,edge=F)
b3<-bw.scott(n84p)
b4<-bw.CvL(n84p)
b1
b2
b3
b4
plot(b1)
plot(b2)
plot(b4)


den <- density(n84p, sigma = b1,edge=F)
plot(den)
plot(n84p, add = TRUE, cex = 0.5)

den <- density(n84p, sigma = b2,edge=F)
plot(den)
plot(n84p, add = TRUE, cex = 0.5)

den <- density(n84p, sigma = b3[1],edge=F)
plot(den)
plot(n84p, add = TRUE, cex = 0.5)

den <- density(n84p, sigma = b3[2],edge=F)
plot(den)
plot(n84p, add = TRUE, cex = 0.5)

den <- density(n84p, sigma = b4,edge=F)
plot(den)
plot(n84p, add = TRUE, cex = 0.5)

den <- density(n84p, sigma = 10,edge=F)
plot(den)
plot(n84p, add = TRUE, cex = 0.5, pch=16)



##############
#test of CSR #
##############


#Chi-square test

### Bei

M <- quadrat.test(bei, nx = 4, ny = 2)
M
plot(bei)
plot(M, add = TRUE)

### n84p

M <- quadrat.test(n84p, nx = 2, ny = 4)
M
plot(n84p)
plot(M, add = TRUE)
axis(1,at=c(seq(0,150,by=25)),pos=c(-100,0))
axis(2,at=c(seq(-100,150,by=25)),pos=c(-10,0))

# An other option
M<-quadrat.test(n84p, nx = 2, ybreaks = c(-100,0,25,50,100))
M
plot(n84p, cex = 0.5, pch = "+")
plot(M, add = TRUE, cex = 2)
axis(1,at=c(seq(0,150,by=25)),pos=c(-100,0))
axis(2,at=c(seq(-100,150,by=25)),pos=c(-10,0))


# In case expected<5- P-value Montecarlo simulation
M<-quadrat.test(n84p, nx = 2, ybreaks = c(-100,-50,-25,0,25,50,100))
M
plot(n84p, cex = 0.5, pch = "+")
plot(M, add = TRUE, cex = 1.5)
axis(1,at=c(seq(0,150,by=25)),pos=c(-100,0))
axis(2,at=c(seq(-100,150,by=25)),pos=c(-10,0))


obs<-M$observed

# Compute the area of each subregion
Z<-quadratcount(n84p, nx = 2, ybreaks = c(-100,-50,-25,0,25,50,100))
quadratsM <- tiles(as.tess(Z))
quadrat.areas <- unlist(lapply(quadratsM, area.owin))

# Expected proportions
p.exp<-quadrat.areas/(sum(quadrat.areas))

chisq.test(M$observed,p=p.exp)$expected
chisq.test(M$observed,p=p.exp,simulate.p.value = T)



## NZ trees
data(nztrees)
plot(nztrees)
M <- quadrat.test(nztrees, nx = 3, ny = 2)
M
plot(nztrees)
plot(M, add = TRUE)


#Kolmogorov-Smirnov test

### Bei
KS=cdf.test(bei,covariate="x")
KS
plot(KS)

KS=cdf.test(bei,covariate="y")
KS
plot(KS)


### nucli 84
KS=cdf.test(n84p,covariate="x")
KS
plot(KS)

KS=cdf.test(n84p,covariate="y")
KS
plot(KS)


### NZ trees

KS=cdf.test(nztrees,covariate="x")
KS
plot(KS)

KS=cdf.test(nztrees,covariate="y")
KS
plot(KS)


################################################
# Exploratory analysis using covariate data ####
################################################

### Bei

# Gradient
data(bei)
Z <- bei.extra$grad
b <- quantile(Z, probs = (0:4)/4)
Zcut <- cut(Z, breaks = b, labels = 1:4)
V <- tess(image = Zcut)
M=quadrat.test(bei, tess = V)
M
plot(M,valuesAreColours=FALSE)
KS <- cdf.test(bei, Z)
plot(KS,style="QQ")


# Elevation
Z <- bei.extra$elev
b <- quantile(Z, probs = (0:4)/4)
Zcut <- cut(Z, breaks = b, labels = 1:4)
V <- tess(image = Zcut)
M=quadrat.test(bei, tess = V)
M
plot(M,valuesAreColours=FALSE)
KS <- cdf.test(bei, Z)
plot(KS,style="QQ")
KS




### nucli84

# Covariables: height and vegetation

# Load data
grid<-read.table("grid.txt",header=T)
grid_veg<-read.table("grid_veg.txt",header=T)

mat<-as.matrix(read.table("height.txt"))
height<-im(mat,grid$x,grid$y)
plot(height,axis=T)
plot(n84p, add=T, cex=0.5, pch=16)


mat<-as.matrix(read.table("veg.txt"))
veg<-im(mat,grid_veg$x,grid_veg$y)
plot(veg,axis=T)
plot(n84p, add=T, cex=0.5)



# Height

Z<-height
b<-quantile(Z, probs = (0:4/4),na.rm = T)
b<-quantile(Z[n84p$window], probs = (0:4/4),na.rm = T)
b

Zcut <- cut(Z, breaks = b, labels = 1:4)
V <- tess(image = Zcut)
plot(V,valuesAreColours=FALSE)
plot(n84p, add = TRUE, pch = 16)

plot(height)
plot(n84p, add = TRUE, pch = 16)

M=quadrat.test(n84p, tess = V)
M
plot(M,valuesAreColours=FALSE)

KS <- cdf.test(n84p, Z)
plot(KS,style="QQ")


# Vegetation
Z<-veg
b<-quantile(Z, probs = (0:4/4),na.rm = T)
b<-quantile(Z[n84p$window], probs = (0:4/4),na.rm = T)
b

Zcut <- cut(Z, breaks = b, labels = 1:4)
V <- tess(image = Zcut)
plot(V,valuesAreColours=FALSE)
plot(n84p, add = TRUE, pch = 16)

plot(veg)
plot(n84p, add = TRUE, pch = 16)

M=quadrat.test(n84p, tess = V)
M

plot(M,valuesAreColours=FALSE)
KS <- cdf.test(n84p, Z)
plot(KS,style="QQ")






#######################################################
# Maximum likelihood estimation for Poisson process ###
#######################################################


# fit an homogeneous poisson model
model1=ppm(bei, ~1)
model1
log(lamb)
plot(model1)

# inhomogeneous Poisson model
model2=ppm(bei, ~x + y)
model2
plot(model2)
anova(model1,model2,test="Chi")

# Model with spatial covariate
grad <- bei.extra$grad
elev=bei.extra$elev

model3=ppm(bei, ~x+y+grad)
model3
plot(model3,se=F)
anova(model3,model2,test="LRT")

# Increase when gradient increases with 1. 
exp(model3$coef[4])

model4=ppm(bei, ~x+y+grad+elev)
model4
plot(model4,se=F)
anova(model3,model4,test="LRT")


# Other visualizations
plot(model4,trend=F,se=T,superimpose=F)
plot(model4,trend=T,se=F,how="persp",superimpose=F)
plot(model4,trend=T,se=F,how="image",superimpose=F,col=heat.colors(112))
plot(model4,trend=T,se=F,how="contour",superimpose=F,nlevels=20)


#Predictions
pred=predict(model4, type = "trend")
plot(pred)


# Model selection. AIC
AIC(model1)
AIC(model2)
AIC(model3)
AIC(model4)
# Lowest AIC is the better model. 


# Checking the fitted Poisson model
# Goodness-of-fit

M <- quadrat.test(model4, nx = 4, ny = 2)
M
plot(bei, pch = ".")
plot(M, add = TRUE, cex = 1.5, col = "red")

KS=cdf.test(model4,"y")
KS
plot(KS)

KS=cdf.test(model4,"x")
KS
plot(KS)


#Validation with residuals

plot(predict(model4))
plot(bei, add = TRUE, pch = "+")

sum(eem(model4))/area(bei)


#Smoothed residuals
xx<-diagnose.ppm(model4, which = "smooth",type="pearson",sigma=50)

# Null standard deviation = 1/(2*bandwidth*sqrt(pi))
3/(2*50*sqrt(pi))


# Lurking variable plot
lurking(model4, expression(x), type = "raw")
lurking(model4, expression(y), type = "raw")
lurking(model4, grad, type = "raw")
lurking(model4, elev, type = "raw")

lurking(model4, grad, type = "raw", cumulative = FALSE)
lurking(model4, elev, type = "raw", cumulative = FALSE)

# Categorize elevation

elev_cut <- cut(elev, breaks = c(0,135,145,155,200))
V <- tess(image = elev_cut)
plot(V,valuesAreColours=FALSE)
plot(bei, add = TRUE, pch = 16)

model5=ppm(bei, ~x+y+grad+elev_cut)
model5
plot(model5,se=F)
anova(model5,test="LRT")
AIC(model4)
AIC(model5)


#Validation with residuals

plot(predict(model5))
plot(bei, add = TRUE, pch = "+")


#Smoothed residuals

par(mfrow=c(1,2))
diagnose.ppm(model4, which = "smooth",type="pearson",sigma=30)
diagnose.ppm(model5, which = "smooth",type="pearson",sigma=30)
par(mfrow=c(1,1))

# Lurking variable plot
lurking(model5, expression(x), type = "raw")
lurking(model5, expression(y), type = "raw")
lurking(model5, grad, type = "raw")
lurking(model5, elev, type = "raw")

lurking(model5, grad, type = "raw", cumulative = FALSE)
lurking(model5, elev, type = "raw", cumulative = FALSE)


model6=ppm(bei, ~x+y+grad+elev_cut+grad*elev_cut)
model6
plot(model6,se=F)
plot(bei,add=T)
anova(model6,model5,test="LRT")
AIC(model5)
AIC(model6)


#Validation with residuals

plot(predict(model6))
plot(bei, add = TRUE, pch = "+")


#Smoothed residuals

par(mfrow=c(1,2))
diagnose.ppm(model5, which = "smooth",type="pearson",sigma=30)
diagnose.ppm(model6, which = "smooth",type="pearson",sigma=30)
par(mfrow=c(1,1))

3/(2*30*sqrt(pi))

sum(eem(model4))/area(bei)
sum(eem(model5))/area(bei)
sum(eem(model6))/area(bei)



# Lurking variable plot
lurking(model6, expression(x), type = "raw")
lurking(model6, expression(y), type = "raw")
lurking(model6, grad, type = "raw")

par(mfrow=c(1,2))
lurking(model6, elev, type = "raw")
lurking(model4, elev, type = "raw")

par(mfrow=c(1,1))

lurking(model6, grad, type = "raw", cumulative = FALSE)
lurking(model6, elev, type = "raw", cumulative = FALSE)





### Nucli 84


# fit an homogeneous poisson model
model1=ppm(n84p, ~1)
model1
plot(model1)

# inhomogeneous Poisson model
model2=ppm(n84p, ~x + y)
model2
plot(model2, se=F)
anova(model1,model2,test="LRT")

# Model with spatial covariate

plot(height)

plot(n84p, add=T)

plot(veg)
plot(n84p, add=T)

model3=ppm(n84p, ~x+y+height)
model3
plot(model3,se=F,ngrid=c(256,256))
anova(model3,model2,test="LRT")



model4=ppm(n84p, ~height)
model4
plot(model4,se=F,ngrid=c(256,256))
anova(model3,model4,test="LRT")

model5=ppm(n84p, ~height+veg)
model5
plot(model5,se=F,ngrid=c(256,256))
anova(model5,model4,test="LRT")


#Predictions
pred=predict(model5, type = "trend",ngrid=c(512,512))
plot(pred)


# Model selection. AIC
AIC(model1)
AIC(model2)
AIC(model3)
AIC(model4)
AIC(model5)



# Checking the fitted Poisson model
# Goodness-of-fit

M <- quadrat.test(model5, nx = 2, ybreaks = c(-100,0,25,50,100))
M

plot(n84p, pch = ".")
plot(M, add = TRUE, cex = 1.5, col = "red")

# Compute the area of each subregion
Z<-quadratcount(n84p, nx = 2, ybreaks = c(-100,0,25,50,100))
quadratsM <- tiles(as.tess(Z))
quadrat.areas <- unlist(lapply(quadratsM, area.owin))

# Expected proportions
p.exp<-M$expected/sum(M$expected)

chisq.test(M$observed,p=p.exp)$expected
chisq.test(M$observed,p=p.exp,simulate.p.value = T)



KS=cdf.test(model5,"y")
KS
plot(KS)

KS=cdf.test(model5,"x")
KS
plot(KS)

theight=cdf.test(model5,height)
theight
plot(theight,style="QQ")


tveg=cdf.test(model5,veg)
tveg
plot(tveg,style="QQ")


#Validation with residuals

plot(predict(model5,ngrid=c(512,512)) )
plot(n84p, add = TRUE, pch = "+")


#Smoothed residuals

diagnose.ppm(model5, which = "smooth",type="pearson",sigma=5)
3/(2*5*sqrt(pi))


# Lurking variable plot
lurking(model5, height, type = "raw")
lurking(model5, veg, type = "raw")

lurking(model5, expression(height+veg), type = "raw")


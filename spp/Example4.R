#############################
### Case control studies ####
#############################

library(spatstat)

# Exploratory analysis
data(chorley)
chorley.extra$plotit()

plot(chorley$window)
points(split(chorley)$lung, pch = 16, col = 1)
points(split(chorley)$larynx, pch = 16, col = 2)
points(chorley.extra$incin, pch = 1, cex = 2, col = 3)

summary(chorley)
plot(split(chorley))
plot(density(split(chorley), sigma=1.5), ribbon = FALSE)



#########################################
# Spatial variation of the relative risk
#########################################

library(sparr)

cases <- split(chorley)$larynx
controls <- split(chorley)$lung

# Intensities
chp <- risk(cases,controls, log=F,intensity=T)
# Cases
plot(chp$f, ribbon=F)

# Controls
plot(chp$g, ribbon=F)

# Relative risk
plot(chp$rr)

# Reference value
cases$n/controls$n


# Log Intensities
chp <- risk(cases,controls, intensity=T)
# Cases
plot(chp$f, ribbon=F)

# Controls
plot(chp$g, ribbon=F)

# Relative risk
plot(chp$rr)

# Reference value
log(cases$n/controls$n)

# Sealevel --> yellow
max(chp$rr)
min(chp$rr)
plot(chp$rr,col=beachcolours(c(-8,-1), sealevel = log(cases$n/controls$n)))


# Log Densities

chp <- risk(cases,controls)

# Cases
plot(chp$f, ribbon=F)
# Controls
plot(chp$g, ribbon=F)

# Relative risk
plot(chp$rr)

max(chp$rr)
min(chp$rr)
plot(chp$rr,col=beachcolours(c(-6,2), sealevel = 0))




########
# Test #
########

rho0<-0 # For Log-densities
chp <- risk(cases,controls)

cellsize<-chp$rr$xstep*chp$rr$ystep
ratiorho <- cellsize*sum((chp$rr$v-rho0)^2,na.rm=T)

# Permutation function

perm_rr<-function(){
  new_ch<-rlabel(chorley)
  new_cases <- split(new_ch)$larynx
  new_controls <- split(new_ch)$lung
  new_chp <- risk(new_cases,new_controls)
  cellsize<-new_chp$rr$xstep*new_chp$rr$ystep
  ratio_perm <- cellsize*sum((new_chp$rr$v-rho0)^2,na.rm=T)
  ratio_perm
}

nsim<-99
set.seed(2021)
rperm<-sapply(1:99,function(i) perm_rr())

# P-value
(sum(rperm > ratiorho)+1)/(nsim+1)

plot(density(rperm))
ratiorho


# P-values Contour Plot
chp <- risk(cases,controls,tolerate = T,adapt=T)
plot(chp)
plot(chp$P)
chp_p<-cut(chp$P,breaks=c(0,0.05,0.1,0.9,0.95,1))
V <- tess(image = chp_p)
plot(V,valuesAreColours=FALSE)


?risk


# Permutation function

perm_rr<-function(){
  new_ch<-rlabel(chorley)
  cases <- split(new_ch)$larynx
  controls <- split(new_ch)$lung
  new_chp <- risk(cases,controls)
  cellsize<-new_chp$rr$xstep*new_chp$rr$ystep
  ratio_perm <- cellsize*sum((new_chp$rr$v-rho0)^2,na.rm=T)
  ratio_perm
}

nsim<-99
set.seed(2021)
rperm<-sapply(1:99,function(i) perm_rr())

# P-value
(sum(rperm > ratiorho)+1)/(nsim+1)



# Binary regression using gam
# Covariate: Distance to focus
dis=sqrt(((chorley$x-chorley.extra$incin$x)^2)+((chorley$y-chorley.extra$incin$y)^2))
ch.data=data.frame(x1=chorley$x,x2=chorley$y,dis)
ch.data$y=abs(as.integer(chorley$marks)-2)
head(ch.data)
library(mgcv)
ch.gam<-gam(y~1+s(x1,x2), data=ch.data, family=binomial)
summary(ch.gam)

ch.gam2<-gam(y~1+dis+s(x1,x2), data=ch.data, family=binomial)
summary(ch.gam2)


# Point source pollution
library(splancs)

rho0<-cases$n/controls$n

D2_mat <- as.matrix(dis^2)
expsource2<-tribble(ccflag=ch.data$y, vars=D2_mat, rho=rho0, alphas=1, betas=1)
print(expsource2)


# Assessment of general spatial clustering


s=seq(0,5,by=1)
khcases<-Kest(cases, r=s, correction="iso")
khcontrols<-Kest(controls, r=s, correction="iso")

plot(khcases, legend=F)
lines(khcontrols$r,khcontrols$iso,lty=2)

# Difference of k functions with envelope
Xppp=chorley
Kdif<-function(Xppp, r, cr="iso")
{
  k1<-Kest(Xppp[marks(Xppp)=="larynx"], r=r, correction=cr)
  k2<-Kest(Xppp[marks(Xppp)=="lung"], r=r, correction=cr)
  D=k1[[cr]]-k2[[cr]]
  res<-data.frame(r=r, D=D)
  return(fv(res, valu="D", fname="D"))
}

nsim<-39
envKdif<-envelope(chorley, Kdif, r=s ,nsim=nsim, nrank=1,
                  savefuns=TRUE,simulate=expression(rlabel(chorley)))


plot(envKdif,legend=F)

# Test

# Extract the simulated functions
simfuns<-as.data.frame(attr(envKdif, "simfuns"))[,-1]
#Compute diagonal of var-cov matrix of dif. K-functions
khcovdiag<-apply(simfuns, 1, var)


#Test
T0<-sum( ((khcases$iso-khcontrols$iso)/sqrt(khcovdiag))[-1])
T_pm<-apply(simfuns, 2, function(X){
  sum((X/sqrt(khcovdiag))[-1])
})

plot(density(T_pm))
abline(v=T0)

pvalue<-2*(sum(abs(T_pm)>abs(T0))+1)/(nsim+1)
pvalue


# Alternative option for limits
plot(envKdif, legend=F)
lines(s, -1.96*sqrt(khcovdiag), lty=2)
lines(s, +1.96*sqrt(khcovdiag), lty=2)



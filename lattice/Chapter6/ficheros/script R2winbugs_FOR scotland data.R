library(R2WinBUGS)

setwd("/home/ajo/gitRepos/spatial/lattice/Chapter6")
###########Scottish data analysis########################################
#Data is in winbugs needs to export into Rpakacge 

#The  rates of lip cancer in 56 counties in Scotland have been analysed 
#by Clayton and Kaldor (1987) and Breslow and Clayton (1993).  
#The form of the data includes the observed and expected cases (expected 
#numbers based on the population and its age and sex distribution in the county), 
#a covariate measuring the percentage of the population engaged in agriculture, 
#fishing, or forestry. 

N <-56
O <- c(9,39,11,9,15,8,26,7,6,20, 13,5,3,8,17,9,2,7,9,7,16,31,11,7,19,15,7,10,16,11, 5,3,7,8,11, 9,11, 8,6,4,
       10, 8,2,6,19,3,2,3,28,6,1,1,1,1,0,0)

E <- c( 1.4, 8.7, 3.0, 2.5, 4.3, 2.4, 8.1,2.3, 2.0, 6.6,4.4, 1.8, 1.1, 3.3,7.8, 4.6, 1.1,4.2,5.5, 4.4,
        10.5,22.7,8.8,5.6,15.5,12.5,6.0,9.0,14.4,10.2,4.8,2.9,7.0, 8.5,12.3,10.1,12.7,9.4,7.2,5.3,
        18.8,15.8,4.3,14.6,50.7,8.2,5.6,9.3,88.7,19.6,3.4,3.6,5.7,7.0,4.2,1.8)

X <- c(16,16,10,24,10,24,10, 7, 7,16, 7,16,10,24, 7,16,10, 7, 7,10, 7,16,10, 7, 1, 1, 7, 7,10,10,7,24,10, 7, 7, 0,10, 1,16, 0,1,16,16, 0, 1, 7, 1, 1, 0, 1,1, 0, 1, 1,16,10)


#Creation dataset.

data.lips<-list(N=N,O=O,E=E,X=X/10)

#File with model
escocia.Poisson<-paste(getwd(),"ficheros/escocia-Poisson-model.odc",sep="")
escocia.inits<-list(
                list(alpha0=0.00000E+00, alpha1=0.00000E+00),
                list(alpha0=-5.00000E-01, alpha1=-2.00000E+00),
                list(alpha0=rnorm(1),alpha1=rnorm(1))
            )


#Define file
# TRIED TO GET IT TO WORK FOR MY OS, BUT GADD IKKE MER.
path.winbugs1<-getwd()#paste("~/Documents/winbugs14_full_patched/WinBUGS14")
path.winbugs2<-paste("~/Documents/winbugs14_full_patched/WinBUGS14") #path PC of UPC university

mybugsdir<-"D:\\demo_WinBugs\\"

results.Poisson<-bugs(data=data.lips,inits=escocia.inits,
                  parameters.to.save = c("RR","alpha0","alpha1"),
                  model.file=escocia.Poisson,n.chains=3,n.iter=2000,n.burnin = 1000,
                  bugs.directory=path.winbugs2,debug=FALSE,codaPkg=FALSE,working.directory=mybugsdir)

print(results.Poisson)

###################Heterogenity model:unstructured overdispersion

escocia.heterogenity<-paste(getwd(),"//escocia-heterogeneity-model.odc",sep="")
escocia.inits.heterogenity<-list(
                          list(alpha0=0.00000E+00, alpha1=0.00000E+00,h=rnorm(56),tau.h=exp(rnorm(1))),
                          list(alpha0=-5.00000E-01, alpha1=-2.00000E+00,h=rnorm(56),tau.h=exp(rnorm(1))),
                          list(alpha0=rnorm(1),alpha1=rnorm(1),h=rnorm(56),tau.h=exp(rnorm(1))))

results.heterogenity<-bugs(data=data.lips,inits=escocia.inits.heterogenity,
                           parameters.to.save = c("alpha0","alpha1","tau.h","sigma.h"),
                           model.file=escocia.heterogenity,n.chains=3,n.iter=15000,n.burnin = 5000,n.thin = 1,
                           bugs.directory=path.winbugs2,debug=TRUE,codaPkg=FALSE,working.directory=mybugsdir)


results.heterogenity<-bugs(data=data.lips,inits=escocia.inits.heterogenity,
                      parameters.to.save = c("alpha0","alpha1","tau.h","sigma.h"),
                      model.file=escocia.heterogenity,n.chains=3,n.iter=15000,n.burnin = 5000,n.thin = 10,
                      bugs.directory=path.winbugs2,debug=TRUE,codaPkg=FALSE,working.directory=mybugsdir)


print(results.heterogenity)

###################Spatial model:structured overdispersion

escocia.spatial<-paste(getwd(),"//escocia-spatial-model.odc",sep="")

#you need spatial weights.
library(maptools)
library(spdep)
scotland<-readSplus("scotland.txt")
plot(scotland)
#Assigno names a les files del data.frame .
data.lips<-data.frame(O=O,E=E,X=X/10)
row.names(data.lips) <- sapply(slot(scotland, "polygons"), slot, "ID")
nc_scotland <- SpatialPolygonsDataFrame(scotland,data=as(data.lips, "data.frame"))
slot(nc_scotland,"data")



#Option 1
#Spatial weights for neighbours lists obtained from the map exported from WINBUGS

#Poly2nb:Construct neighbours list from polygon list based on contiguous boundaries
xxnb <- poly2nb(nc_scotland)
#nb2listw:The function supplements a neighbours list with spatial weights 
w.scotland<-nb2listw(xxnb, style="W", zero.policy=TRUE)
card(xxnb)
class(w.scotland)
sids.moran.mc<- moran.mc(nc_scotland$residuals, listw=w.scotland, nsim=100,zero.policy=TRUE)
attributes(w.scotland)

###Be carefully with the "adjuntos" vector
adj<-unlist(w.scotland$neighbours)
adj<-adj[adj!=0]
data.lips.spatial<-list(N=N,O=O,E=E,X=X/10,adj=adj,num=card(xxnb),sumNumNeigh=sum(card(xxnb)) )

###Option 2: From WinBUGS

# num = c(3, 2, 1, 3, 3, 0, 5, 0, 5, 4,
#               0, 2, 3, 3, 2, 6, 6, 6, 5, 3,
#               3, 2, 4, 8, 3, 3, 4, 4, 11, 6,
#               7, 3, 4, 9, 4, 2, 4, 6, 3, 4,
#               5, 5, 4, 5, 4, 6, 6, 4, 9, 2,
#               4, 4, 4, 5, 6, 5)
# adj = c(19, 9, 5, 10, 7, 12, 28, 20, 18, 19, 12, 1, 17, 16, 13, 10, 2,  29, 23, 19, 17, 1,
#   22, 16, 7, 2, 5, 3, 19, 17, 7,  35, 32, 31, 29, 25, 29, 22, 21, 17, 10, 7,  29, 19, 16, 13, 9, 7, 56, 55, 33, 28, 20, 4,  17, 13, 9, 5, 1,
#   56, 18, 4,  50, 29, 16, 16, 10, 39, 34, 29, 9, 56, 55, 48, 47, 44, 31, 30, 27,   29, 26, 15, 43, 29, 25,  56, 32, 31, 24, 45, 33, 18, 4,
#   50, 43, 34, 26, 25, 23, 21, 17, 16, 15, 9,  55, 45, 44, 42, 38, 24, 47, 46, 35, 32, 27, 24, 14, 31, 27, 14, 55, 45, 28, 18, 54, 52, 51, 43, 42, 40, 39, 29, 23,
#   46, 37, 31, 14, 41, 37, 46, 41, 36, 35, 54, 51, 49, 44, 42, 30, 40, 34, 23, 52, 49, 39, 34, 53, 49, 46, 37, 36, 51, 43, 38, 34, 30,
#   42, 34, 29, 26, 49, 48, 38, 30, 24, 55, 33, 30, 28, 53, 47, 41, 37, 35, 31, 53, 49, 48, 46, 31, 24, 49, 47, 44, 24, 54, 53, 52, 48, 47, 44, 41, 40, 38,
#   29, 21, 54, 42, 38, 34, 54, 49, 40, 34, 49, 47, 46, 41, 52, 51, 49, 38, 34, 56, 45, 33, 30, 24, 18, 55, 27, 24, 20, 18)
# sumNumNeigh = 234
# data.lips.spatial<-list(N=N,O=O,E=E,X=X/10,adj=adj,num=num,sumNumNeigh=sumNumNeigh) 


escocia.inits.spatial<-list(
  list(alpha0=0.00000E+00, alpha1=0.00000E+00,s=rnorm(56),tau.s=exp(rnorm(1))),
  list(alpha0=-5.00000E-01, alpha1=-2.00000E+00,s=rnorm(56),tau.s=exp(rnorm(1))),
  list(alpha0=rnorm(1),alpha1=rnorm(1),s=rnorm(56),tau.s=exp(rnorm(1))))

results.spatial<-bugs(data=data.lips.spatial,inits=escocia.inits.spatial,
                           parameters.to.save = c("RR","alpha0","alpha1","tau.s","sigma.s"),
                           model.file=escocia.spatial,n.chains=3,n.iter=15000,n.burnin = 5000,n.thin=10,
                           bugs.directory=path.winbugs2,debug=TRUE,codaPkg=FALSE,working.directory=mybugsdir)


print(results.spatial)
results.spatial$sims.list$alpha0
mean(results.spatial$sims.list$alpha0)

###################Convolution model:

escocia.convolution<-paste(getwd(),"//escocia-convolution-model.odc",sep="")
escocia.inits.convo<-list(
  list(alpha0=0.00000E+00, alpha1=0.00000E+00,h=rnorm(56),tau.h=exp(rnorm(1)), s=rnorm(56),tau.s=exp(rnorm(1))),
  list(alpha0=-5.00000E-01, alpha1=-2.00000E+00,h=rnorm(56),tau.h=exp(rnorm(1)),s=rnorm(56),tau.s=exp(rnorm(1))),
  list(alpha0=rnorm(1),alpha1=rnorm(1),h=rnorm(56),tau.h=exp(rnorm(1)),s=rnorm(56),tau.s=exp(rnorm(1))))


#for check convergence
results.convo<-bugs(data=data.lips.spatial,inits=escocia.inits.convo,
                    parameters.to.save = c("alpha0","alpha1","tau.s","sigma.s","tau.h","sigma.h"),
                    model.file=escocia.convolution,n.chains=3,n.iter=2000,n.burnin = 500,
                    bugs.directory=path.winbugs2,debug=TRUE,codaPkg=FALSE,working.directory=mybugsdir)



results.convo<-bugs(data=data.lips.spatial,inits=escocia.inits.convo,
                      parameters.to.save = c("RR","alpha0","alpha1","tau.s","sigma.s","tau.h","sigma.h","frac.spatial","s2.marginal", "sigma2.h"),
                      model.file=escocia.convolution,n.chains=3,n.iter=150000,n.thin=10,n.burnin = 5000,
                      bugs.directory=path.winbugs1,debug=TRUE,codaPkg=FALSE,working.directory=mybugsdir)




print(results.convo)
print(results.spatial)
print(results.heterogenity)

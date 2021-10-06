# Encontrado en la página web de Montserrat Fuentes:
# http://www.stat.ncsu.edu/people/fuentes/st790m/lectures/lectures.html
# http://www.stat.ncsu.edu/people/fuentes/st790m/lectures/Rose.txt
# Citado como
# "Handout with the R code to produce a Rose Diagram by M.S. Park."
#
###########################################################################################
#
# Rose diagram gives us useful information about anisotropy (i.e., rotation and stretch).
# In the output figure, length of a line segment is the relative distance to touch your 
# preassigned semivariogram value (crit.val).
#
# 1. data.var : a vector of the variable
# 2. data.cds : an n x 2 matrix containing coordinates
# 3. max.dist : a numerical value defining the maximum distance for the variogram
# 4. numcases : the number of bins for the variogram
# 5. numdirec : the number of directional angles 
# 6. poly.tnd : specifies the mean structure of the model ("cte", "1st" or "2nd")
# 7. crit.val : should be less than min(max.semi)
# 	   otherwise, crit.val is changed to min(max.semi)*0.99 automatically.
#
###########################################################################################

rose.diagram <- function(data.var,data.cds,max.dist,numcases,numdirec,poly.tnd,crit.val){
	semivarg <- matrix(0,numcases*numdirec,1)
	distance <- matrix(0,numcases*numdirec,1)
	directio <- matrix(0,numcases*numdirec,1)
	sample.n <- matrix(0,numcases*numdirec,1)
	max.semi <- matrix(0,numdirec,1)
	for (i in 1:numdirec){
		temp <- variog(data=data.var,coords=data.cds,bin.cloud=T,uvec=seq(0,max.dist,l=numcases),direction=i*pi/numdirec,tol=pi/numdirec,trend=poly.tnd)
		max.semi[i,] <- max(as.matrix(temp$v))
		semivarg[((i-1)*numcases+1):(i*numcases),] <- as.matrix(temp$v)
		distance[((i-1)*numcases+1):(i*numcases),] <- as.matrix(temp$u)
		sample.n[((i-1)*numcases+1):(i*numcases),] <- as.matrix(temp$n)
		directio[((i-1)*numcases+1):(i*numcases),] <- rep(temp$direction,numcases)
	}
	if (crit.val > min(max.semi)){
		cat("Warning: crit.val should be less than",min(max.semi),"\n")
		crit.val <- min(max.semi)*0.99
	}
	rose.len <- rep(0,numdirec)
	rose.dir <- rep(0,numdirec)
	j <- 1
	for (j in 1:numdirec){
		for (i in ((j-1)*numcases+1):(j*numcases-1)){
			prod <- (semivarg[i] - crit.val)*(semivarg[i+1] - crit.val)
			if (directio[i] != directio[i+1]) prod <- 0
			if (j > 1){
				if (rose.dir[j-1] == directio[i]) prod <- 0
			}
			if (prod < 0){
				rose.len[j] <- distance[i]+(crit.val-semivarg[i])*(distance[i+1]-distance[i])/(semivarg[i+1]-semivarg[i])
				rose.dir[j] <- directio[i]
				j <- j+1
			}
		}
	}
	max.length <- max(rose.len)
	rose.len <- rose.len/max.length
	
	par(mar=c(0.5,1.5,1.5,0.5))
	plot(0,0,type="n",axes=F,xlim=c(-1,1),ylim=c(-1,1),xlab="",ylab="")
	for (j in 1:(length(rose.len))){
		segments(0,0,(rose.len[j]*sin(rose.dir[j])),(rose.len[j]*cos(rose.dir[j])))
		segments(0,0,-(rose.len[j]*sin(rose.dir[j])),-(rose.len[j]*cos(rose.dir[j])))
	}
	mtext("N-S Direction",side=2)
	mtext("E-W Direction",side=3)
}

# instrucción en el fichero original:
#rose.diagram(data.var=pm.dat$pm,data.cds=pm.dat[,4:5],max.dist=800,numcases=9,numdirec=8,poly.tnd="1st",crit.val=50)

# # Usando un conjunto de datos diferente:
# scallops <- read.table("Scallops_R.dat",head=TRUE,sep=" ",dec=".")
# scallops$long <- -scallops$long
# scallops$lgcatch <- log(1 + scallops$tcatch)
# geosca <- as.geodata(scallops[,-c(1,2)],coords.col = 2:1, data.col=3:6)
# geoscalg <- as.geodata(scallops[,-c(1,2)],coords.col = 2:1,data.col=6)
# summary(geoscalg)
# rose.diagram(data.var=scallops$lgcatch,data.cds=scallops[,4:3],max.dist=1,numcases=9,numdirec=8,poly.tnd="1st",crit.val=50)


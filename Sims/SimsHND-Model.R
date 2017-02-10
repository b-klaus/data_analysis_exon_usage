
## ----, echo=FALSE--------------------------------------------------------
print(date())


## ----setup---------------------------------------------------------------
library(knitr)
library(ConcaveFDR)
library(ggplot2)
library(grid)
library(locfdr)
library(fdrtool)
library(mixfdr)
library(BiocParallel)
library(plyr)


opts_chunk$set(fig.width=14, fig.height=10, cache=T, 
error=FALSE, cache.lazy=FALSE)
options(digits = 2, scipen = 10)

multicoreParam <- MulticoreParam(workers = round(parallel::detectCores() * 0.75))
multicoreParam
register(multicoreParam)



## ----setting parameters--------------------------------------------------

### number of reps
B = 1000 
### number of stats
d = 1000
### prop of null
eta0 = 0.8
### sigma
sigma = 1

### grid for true fdrs
s = seq(0.01,10,by=0.01) 

## FDR functions for the HND model
Fdr.HNDscore = function(x)  ConcaveFDR:::HND(abs(x),  ConcaveFDR:::eta02k(eta0)) 
fdr.HNDscore = function(x)  ConcaveFDR:::hnd(abs(x),  ConcaveFDR:::eta02k(eta0))

## true FDRs
truefdr =  fdr.HNDscore(s) 
trueFDR = Fdr.HNDscore(s)


### test it
score <- get.random.HNDscore(1000, eta0)
qplot(score)



## ----simulate HND scores-------------------------------------------------
##  d z-scores in every column 
  HNDMat  <- t(laply(bplapply(seq_len(B), function(i){ 
    get.random.HNDscore(d)
	}), identity))

#save(HNDMat, file = "HNDMat.rda")

hist(get.random.HNDscore(1000), breaks = 100)

hist(HNDMat[1:1000,1], breaks = 100)
z <- HNDMat[1:1000,1]


## ----helper functions----------------------------------------------------


mix.fdr.fun = function(z,P,J,empnull=FALSE){
  theonull = !empnull
	a = NA
	if(theonull) a = 1
	e = mixFdr(z[1:d],J=J, P = P, theonull = theonull, noiseSD = a, calibrate = 
	FALSE, plots = FALSE, nocheck = TRUE)
	ind = abs(e$mu - e$mu[1])<=.1
	tmp = (fdrMixModel(s,e,ind))
	sigma = e$sigma[1]
	eta_0 = e$pi[1]
	return ( c (tmp , eta_0, sigma ))
	#return(tailFDRMixModel(s, e, ind)$tailFDR)
}

#test = mix.fdr.fun(z,P=50,J=3,empnull = TRUE)

loc.fdr.fun = function(z){

	loc = locfdr(z[1:d], plot = 0)
	eta_0 = loc$fp0["mlest","p0"]
  	sigma = loc$fp0["mlest","sigma"]
	
	res = loc$fdr	
	tmp = approx(x = z[1:d], y= res, xout = s)$y
	tmp[which(is.na(tmp))] = 0.0 #replace NAs by 0.0

	return ( c (tmp , eta_0, sigma ))
}

#test = loc.fdr.fun(z)
#plot(s,test[1:d])

#test = loc.fdr.fun(z)


fdrtool.fdr.fun = function(z, met = "fndr"){
	tool = fdrtool(z[1:d], "normal",plot=FALSE, verbose =FALSE, cutoff.method = 
	met)	
	res = tool$lfdr		
	
	#sp = splinefun(z, logit(pmax(pmin(res,0.999),.001)), method = "natural")	
	#pmin(1,pmax(0,logistic(sp(s))))	
	tmp = approx(x = z[1:d], y= res, xout = s)$y
	tmp[which(is.na(tmp))] = 0.0 #replace NAs by 0.0

	return ( c (tmp , tool$param[3], tool$param[5] ))
}

#test = fdrtool.fdr.fun(z)
#plot(s,test[1:d], type = "l")


##theo = FALSE; classic = FALSE
ConcaveFDR.fun = function(z){

  res = ConcaveFDR(z[1:d], "normal",plot=FALSE, verbose =FALSE, cutoff.method = 
  "smoothing")	
  fdr = res$lfdr.log


  tmp = approx(x = z[1:d], y= fdr, xout = s)$y  
  tmp[which(is.na(tmp))] = 0.0 #replace NAs by 0.0
  return ( c (tmp , res$param[3], res$param[5] ))
  
}

#test = ConcaveFDR.fun(z)
#plot(s,test[1:d])


HND.sigma.fdr.fun = function(z){
	
	param.out = hnd.fitnull(z[1:d],"normal")
	eta0 = param.out$eta0
     	sd = param.out$sd
	
	pval = 2-2*pnorm(abs(z/sd))
	

	### compute FDR values  
	az = qnorm(1-pval/2)
  	k = ConcaveFDR:::eta02k(eta0)

  	#Fdr = HND(az, k)
  	fdr = ConcaveFDR:::hnd(az, k)


tmp=approx(x = z[1:d], y= fdr , xout = s)$y
tmp[which(is.na(tmp))] = 0.0 #replace NAs by 0.0
return( c(tmp, eta0 , sd) )
}

#test = HND.sigma.fdr.fun(z)
#plot(s,test[1:d])


#external empirical null

HND.fdr.fun = function(z){
	
	 source("/home/bklaus/uni/Projects/fdrtool-dev/TESTConcave/log.fdr.R") 
    resFDR = log.fdr.fda(z[1:d], theo = FALSE, classic = FALSE)
     eta0 =  resFDR$eta0
     sd = resFDR$sd
     pval =resFDR$pval	

	# pv = pval.log
	### compute FDR values  
	az = qnorm(1-pval/2)
  k =  ConcaveFDR:::eta02k(eta0)

  	#Fdr = HND(az, k)
  	fdr = ConcaveFDR:::hnd(az, k)


tmp=approx(x = z[1:d], y= fdr , xout = s)$y
tmp[which(is.na(tmp))] = 0.0 #replace NAs by 0.0
return( c(tmp, eta0 , sd) )
}

#test = HND.fdr.fun(z)
#plot(s,test[1:d])

####################################





## ----SessionInfo, echo=FALSE---------------------------------------------

sessionInfo()




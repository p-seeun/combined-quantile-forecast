source("./functions.R")

library(fields)

## Load Data 
PM = read.csv("../Data/PM.csv", row.names = 1, check.names=F)
PM = as.matrix(PM)
prds = list(prd1=12697:12744, prd2=13417:13464, prd3=14161:14208, prd4=14905:14952, prd5=15625:15672, prd6=16369:16416,
           prd7=17089:17136, prd8=17833:17880, prd9=18577:18624, prd10=19273:19320, prd11=20017:20064, prd12=20737:20784)
tau.vec = c(0.1,0.3,0.5,0.7,0.9)

## Set h & period
h=1
prd.ind=1

prd = prds[[prd.ind]] 

## Find the number of quantile factors for each tau level
r.vec = c() # r.vec[i] : number of quantile factors of i-th tau level
for(l in 1:length(tau.vec)){
  rhat = c()
  for(j in seq(from=0,to=47,by=6)){
    icr = get_icr_q(PM[(prd[1]-500+j):(prd[1]-1+j),], range=5:15, tau=tau.vec[l])
    ind = which.min(icr$ic)
    rhat = c(rhat, icr$r[ind])
    print(c(l,j))
  }
  r.vec[i] = mean(rhat)
}
r.vec 

## The results of estimated 'r.vec' for each period are saved in data "r.vec_list" as a list of length 12. 
# r.vec_list = readRDS("../Data/r.vec_list.rds")
# r.vec = r.vec_list[[prd.ind]]


## Step 1. Estimate quantile levels of PM ## 

estPM = as.list(1:(length(prd)* 308))
dim(estPM) = c(length(prd), 308)

for(k in 1:308){
  for(i in prd){
    y = PM[(i-h-499):(i-h), k] # 500 * 1
    X = t( PM[(i-h-499):(i-h), -k] ) # 500 * 307
    
    estPM[[i-prd[1]+1,k]] = qfm_frcst(y, X, h=h, tau.vec=tau.vec, r.vec=r.vec) 
  }
  print(k)
}

## Step 2. Find weight for each k ##

# Find state S_t for each t and k
state.tot = matrix(nrow=length(prd), ncol=308)
for(k in 1:308){
 for(i in prd){
 est = estPM[[i-prd[1]+1 , k]]
 s = ifelse(PM[i,k] < est[1], 1, ifelse(PM[i,k] < est[2], 2, ifelse(PM[i,k] < est[4], 3, ifelse(PM[i,k] < est[5], 4, 5))))
 state.tot[i-prd[1]+1, k] = s
 }
}

weight = list()
loc = read.csv("../Data/station_loc.csv", row.names = 1, check.names=F)
# i-th row of dist.order lists stations in the order of distance from i-th station
dist.order<-matrix(nrow=308, ncol=307)
for(i in 1:308){
  dist.i <- drop(rdist.earth(x1=loc[i,], x2=loc, miles=FALSE)) 
  dist.order[i,]<-order(dist.i)[-1]
}

# Compute Markov chain probability
for(k in 1:308){
 close<-c(k, as.numeric(dist.order[k,1:9]))
 state = state.tot[,close]
 markov<-matrix(0, nrow=length(tau.vec), ncol=length(tau.vec))
 for(t in 1:(length(prd)-h)){
   for(r in 1:ncol(state)){
     pre = state[t,r]; post = state[t+h,r]
     markov[pre,post] = markov[pre,post] + 1
   }
 }
 mc<-matrix(0, nrow=length(tau.vec), ncol=length(tau.vec))
 for(i in 1:length(tau.vec)){
   mc[i,] = markov[i,]/sum(markov[i,])
 }
 weight[[k]] = mc
}

# Get combined estimator for each time and station
res = matrix(nrow=length(prd), ncol=308)
for(k in 1:308){
  for(i in prd){
    if(i-prd[1]<h){
      res[i-prd[1]+1, k] = estPM[[i-prd[1]+1, k]] %*% c(0.1,0.2,0.4,0.2,0.1) #combined est
    }
    else{
      s = state.tot[i-prd[1]+1-h, k]
      res[i-prd[1]+1, k] = estPM[[i-prd[1]+1, k]] %*% weight[[k]][s,] #combined est
    }
  }
}

# MAE for prd
mae(res, PM[prd,])



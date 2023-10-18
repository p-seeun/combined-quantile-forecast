source("./functions.R")
library(fields)

### Load Data 
PM = read.csv("../Data/PM.csv", row.names = 1, check.names=F)
PM = as.matrix(PM)
prds = list(prd1=12697:12744, prd2=13417:13464, prd3=14161:14208, prd4=14905:14952, prd5=15625:15672, prd6=16369:16416,
           prd7=17089:17136, prd8=17833:17880, prd9=18577:18624, prd10=19273:19320, prd11=20017:20064, prd12=20737:20784)

### Set h & period index
h = 2
prd.ind = 2

### Analysis starts
prd = prds[[prd.ind]] 

##### Naive, AR, ARIMA method #####

res.naive = matrix(nrow=length(prd), ncol=308)
res.ar = matrix(nrow=length(prd), ncol=308)
res.arima = matrix(nrow=length(prd), ncol=308)

for(k in 1:308){
  for(i in prd){
    y = PM[(i-h-499):(i-h), k]
    res.naive[i-prd[1]+1, k] = naive_frcst(y, h)
    res.ar[i-prd[1]+1, k] = ar_frcst(y, h)
    res.arima[i-prd[1]+1, k] = arima_frcst(y, h)
  }
  print(k)
}

# MAE results
mae(res.naive, PM[prd,])
mae(res.ar, PM[prd,])
mae(res.arima, PM[prd,])


##### SW(2002) method #####

# Estimate number of factors r
rhat = c()
for(i in seq(from=0,to=47,by=6)){
  icr = get_icr(PM[(prd[1]-499+i):(prd[1]+i),], 10:30)
  ind = which.min(icr$ic)
  rhat = c(rhat, icr$r[ind])
}
r = round(mean(rhat))

res.sw = matrix(nrow=length(prd), ncol=308)
for(k in 1:308){
  for(i in prd){
    y = PM[(i-h-499):(i-h), k]
    X = t( PM[(i-h-499):(i-h), -k] ) #n * T
    
    res.sw[i-prd[1]+1, k] = sw_frcst(y=y, X=X, h=h, r=r) 
  }
  print(k)
}

# MAE results
mae(res.sw, PM[prd,])


## 3. Near-stations 
# Load longitude, latitude of stations
loc = read.csv("../Data/station_loc.csv", row.names = 1, check.names=F)

# i-th row of dist.order lists stations in the order of distance from i-th station
dist.order<-matrix(nrow=308, ncol=307)
for(i in 1:308){
  dist.i <- drop(rdist.earth(x1=loc[i,], x2=loc, miles=FALSE)) 
  dist.order[i,]<-order(dist.i)[-1]
}

res.near = matrix(nrow=length(prd), ncol=308)

near_frcst <- function(y, fulldata, near, h){ # full data : Te*n matrix
  Te = nrow(fulldata)
  lm = lm(y[(1+h):Te] ~ fulldata[1:(Te-h), near])
  est = sum ( coef(lm) * c(1, fulldata[Te, near]), na.rm=T )
  return(est)
}


for(k in 1:308){
  near <- c(k, as.numeric(dist.order[k,1:r]) ) # Use the same number of r in SW(2002)
  for(i in prd){
    y = PM[(i-h-499):(i-h), k] 
    fulldata = PM[(i-h-499):(i-h),] # Te*n matrix
    res.near[i-prd[1]+1, k] = near_frcst(y=y, fulldata=as.matrix(fulldata), near=near, h=h)
  }
  print(k)
}

# MAE results
mae(res.near, PM[prd,])



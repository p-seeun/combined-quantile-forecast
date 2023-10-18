source("./functions.R")

set.seed(1)
iter = 500

# generating functions for three error distributions
rgen <- function(dist, n){ 
  if(dist=="norm"){ return(rnorm(n)) }
  if(dist=="t"){ return(rt(n, df=2)) }
  if(dist=="gamma"){ return(rgamma(n, shape=1, scale=5)) }
}

SIM.tot = list()  #save all results here
SIM = as.list(1:9) 
dim(SIM) = c(3,3)
dist.ind = 0

for(dist in c("norm", "t", "gamma")){
dist.ind = dist.ind + 1
for(n.ind in 1:3){
  for(t.ind in 1:3){
  
    n=50*n.ind; Te=50*t.ind 
    
    rsq = matrix(nrow=iter, ncol=3) # adjusted R^2 by quantile factors
    rsq.m = matrix(nrow=iter, ncol=3) # adjusted R^2 by mean factors
    
    result = matrix(nrow=iter, ncol=5) # MAE results by 5 methods
    colnames(result)<-c("naive","ar","arima","sw","proposed")
    
    for(u in 1:iter){
      ## Generate data
      # Factor
      factor<-matrix(nrow=3, ncol=Te)
      factor[,1]<-c(rnorm(1, mean=0, sd=(1-0.8^2)^(-0.5)), 
                    rnorm(1, mean=0, sd=(1-0.5^2)^(-0.5)), 
                    rnorm(1, mean=0, sd=(1-0.2^2)^(-0.5)) )
      for(i in 1:(Te-1)){
        factor[,i+1]<-c(0.8,0.5,0.2)*factor[,i] + rnorm(3)
      }
      #Lambda
      Lambda<-matrix(nrow=n, ncol=3)
      for(j in 1:3){
        Lambda[,j]<-rnorm(n)
      }
      
      #X
      X = matrix(nrow=n, ncol=Te)
      X = Lambda%*%factor + rgen(dist, n*Te)
      #y
      y<-c()
      y[1:Te] = colSums(factor)+rgen(dist,Te)
      real = sum(c(0.8,0.5,0.2)*factor[,Te]+rnorm(3)) #sum of f_{Te+1}
      y[Te+1] = real+rgen(dist,1)
      
      ## Experiment 1 (Table 1)
      # f^{PCA}
      Xf<-scale(t(X), scale=FALSE)
      svd<-svd(Xf)
      lamb<-t(as.matrix(sqrt(n)*t(svd$v)[1:3,]))
      fact<-(Xf %*% lamb)/n
      
      # f^{QFM}
      qfm <- do_qfm(Xf, r.vec=c(3,3,3,3,3), tau.vec = c(0.1,0.3,0.5,0.7,0.9))
      
      # Compute R^2
      lm<-lm(factor[1,]~qfm[[1]]$F[1:Te,]+qfm[[2]]$F[1:Te,]+qfm[[3]]$F[1:Te,]+qfm[[4]]$F[1:Te,]+qfm[[5]]$F[1:Te,])
      rsq[u,1]<-summary(lm)$adj.r.squared
      lm<-lm(factor[1,]~fact)
      rsq.m[u,1]<-summary(lm)$adj.r.squared
      
      lm<-lm(factor[2,]~qfm[[1]]$F[1:Te,]+qfm[[2]]$F[1:Te,]+qfm[[3]]$F[1:Te,]+qfm[[4]]$F[1:Te,]+qfm[[5]]$F[1:Te,])
      rsq[u,2]<-summary(lm)$adj.r.squared
      lm<-lm(factor[2,]~fact)
      rsq.m[u,2]<-summary(lm)$adj.r.squared
      
      lm<-lm(factor[3,]~qfm[[1]]$F[1:Te,]+qfm[[2]]$F[1:Te,]+qfm[[3]]$F[1:Te,]+qfm[[4]]$F[1:Te,]+qfm[[5]]$F[1:Te,])
      rsq[u,3]<-summary(lm)$adj.r.squared
      lm<-lm(factor[3,]~fact)
      rsq.m[u,3]<-summary(lm)$adj.r.squared
      
      ## Experiment 2 (Table 2)
      # Compute MAE
      comb_frcst = qfm_frcst(y[1:Te], X, h=1, tau.vec=c(0.1,0.3,0.5,0.7,0.9), r.vec=rep(3,5), qfm = qfm) %*% c(0.1,0.2,0.4,0.2,0.1)
      result[u,] = abs( c( naive_frcst(y[1:Te], h=1), ar_frcst(y[1:Te], h=1), arima_frcst(y[1:Te], h=1), sw_frcst(y[1:Te], X, h=1, r=3, fact = fact), 
                         comb_frcst) - y[Te+1] )
      print(u)
    }
    SIM[[n.ind,t.ind]] = list( rsq.pca = colMeans(rsq.m), rsq.qfm = colMeans(rsq), 
                               mae = colMeans(result), sd = apply(result,2,sd)/sqrt(iter) )
  }
}
  SIM.tot[[dist.ind]] = SIM
}



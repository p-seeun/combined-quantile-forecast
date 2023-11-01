source("./functions.R")

set.seed(111)
iter = 500

SIM.ext = as.list(1:9) #save all results here
dim(SIM.ext) = c(3,3)

for(n.ind in 1:3){
  for(t.ind in 1:3){
   
    n=50*n.ind; Te=50*t.ind 
   
    result = matrix(nrow=iter, ncol=6) # MAE results by 5 methods
    colnames(result)<-c("naive","ar","arima","sw","proposed","extended")
    
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
      X= Lambda%*%factor
      same = sample(n, n/2)
      X[same,] = X[same,] + rt(n*Te/2, df=2)
      X[-same,] = X[-same,] + rgamma(n*Te/2, shape=1, scale=5)
      #y
      y<-c()
      y[1:Te] = colSums(factor)+rgamma(Te, shape=1, scale=5)
      real = sum(c(0.8,0.5,0.2)*factor[,Te]+rnorm(3)) #sum of f_{Te+1}
      y[Te+1] = real+rgamma(1, shape=1, scale=5)
      
      # f^{PCA}
      Xf<-scale(t(X), scale=FALSE)
      svd<-svd(Xf)
      lamb<-t(as.matrix(sqrt(n)*t(svd$v)[1:3,]))
      fact<-(Xf %*% lamb)/n
      
      # f^{QFM}
      qfm <- do_qfm(Xf, r.vec=c(3,3,3,3,3), tau.vec = c(0.1,0.3,0.5,0.7,0.9))
      
      ## Compute MAE (Table 3)
      comb_frcst = qfm_frcst(y[1:Te], X, h=1, tau.vec=c(0.1,0.3,0.5,0.7,0.9), r.vec=rep(3,5), qfm = qfm) %*% c(0.1,0.2,0.4,0.2,0.1)
      ext_comb_frcst = ext_qfm_frcst(y[1:Te], X, h=1, tau.vec=c(0.1,0.3,0.5,0.7,0.9), r.vec=rep(3,5), qfm = qfm) %*% c(0.1,0.2,0.4,0.2,0.1)
      result[u,] = abs( c( naive_frcst(y[1:Te], h=1), ar_frcst(y[1:Te], h=1), arima_frcst(y[1:Te], h=1), sw_frcst(y[1:Te], X, h=1, r=3, fact = fact), 
                           comb_frcst, ext_comb_frcst ) - y[Te+1] )
    }
    SIM.ext[[n.ind,t.ind]] = list( mae = colMeans(result), sd = apply(result,2,sd)/sqrt(iter) )
  }
}





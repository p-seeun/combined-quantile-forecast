library(quantreg)
library(forecast)
library(hrqglas)

## Functions

# Check function
loss.fun <-function(tau, vec){
  return(mean(sapply(vec,function(x){return(ifelse(x>0, tau*x, (1-tau)*abs(x)))})))
}

# Estimate quantile factors
do_qfm <- function(X,r.vec,tau.vec){ #X should be centered #tau.vec needed 
  N = ncol(X)
  Te = nrow(X)
  temp <- svd(X)
  pred.quant.i = list()
  for(ind in 1:length(tau.vec)){
    print(paste("###########", ind, "th QFM estimation###########", sep=""))
    tau <- tau.vec[ind] 
    r = ifelse(length(r.vec)==1,r.vec,r.vec[ind])
    L.hat = sqrt(N)*(temp$v)[,1:r] #initial
    F.hat = 1/N * X %*% L.hat      #initial
    
    # mean(abs(X-F.hat%*%t(L.hat)))
    L.hat.new <- matrix(0, nrow=N, ncol=r)
    F.hat.new <- matrix(0, nrow=Te, ncol=r)
    
    iter <- 1
    while(1){
      for(i in 1:N){
        L.hat.new[i,] <- rq(X[,i] ~ F.hat-1, tau=tau)$coeff
      }
      for(t in 1:nrow(F.hat)){
        F.hat.new[t,] <- rq(X[t,] ~ L.hat.new-1, tau=tau)$coeff
      }
      q1 = qr(F.hat.new)
      q2 = qr(qr.R(q1)%*%t(L.hat.new))
      F.hat.new = sqrt(Te) * qr.Q(q1) %*% qr.Q(q2)
      L.hat.new = t(qr.R(q2))/sqrt(Te)
      
      if(iter>1 ){
        if(abs(loss.fun(tau,as.vector(X-F.hat.new %*% t(L.hat.new))) - loss.fun(tau,as.vector(X-F.hat %*% t(L.hat))))< 1e-3 ){
          print(loss.fun(tau,as.vector(X-F.hat.new %*% t(L.hat.new))) - loss.fun(tau,as.vector( X-F.hat %*% t(L.hat))))
          print(iter)
          break
        }
      }
      iter <- iter+1
      F.hat <- F.hat.new
      L.hat <- L.hat.new
      if(iter==200){
        print(paste("Max iteration at tau=", tau))
        break;
      }
    }
    pred.quant.i [[ind]] = list(F=F.hat.new, L=L.hat.new)
  }
  return(pred.quant.i)
}

# Naive method
naive_frcst <- function(y, h){
  return(tail(y,1))
}

# AR model
ar_frcst <- function(y, h){
  model_ar<-ar.ols(y, order.max = 6, demean=FALSE, intercept=TRUE)
  pred_ar<-predict(model_ar, n.ahead=h)
  return(pred_ar$pred[h])
}

# ARIMA model
arima_frcst <- function(y, h){
  order = arimaorder(auto.arima(y))
  model_arima <- try(Arima(y, order=order), silent=T)
  
  if ("try-error" %in% class(model_arima) ){
    model_arima <- try(Arima(y, order=order, method="ML"), silent=T)
    if("try-error" %in% class(model_arima) ){
      model_ar<-ar.ols(y, order.max = 6, demean=FALSE, intercept=TRUE)
      pred_ar<-predict(model_ar, n.ahead=h)
      return(pred_ar$pred[h])
    }
    pred_arima <- predict(model_arima, n.ahead=h)
    return(pred_arima$pred[h])
  }
  
  pred_arima <- predict(model_arima, n.ahead=h)
  return(pred_arima$pred[h])
}

# Estimate number of mean factors using IC
get_icr = function(X, range){
  X = scale(X, scale=F)
  svd = svd(X) 
  N = ncol(X); Te = nrow(X)
  icr = data.frame(matrix(nrow=0, ncol=2))
  for(r in range){
    lambda = t(as.matrix(sqrt(N)*t(svd$v)[1:r,])) 
    Factor = (X %*% lambda)/N
    ic = log(sum((X-Factor%*%t(lambda))^2)/(N*Te))+r*(N+Te)/(N*Te)*log(N*Te/(N+Te))
    icr = rbind(icr, c(r,ic))
  }
  colnames(icr) = c("r","ic")
  return(icr)
}

# SW(2002) method
sw_frcst <- function(y, X, h, r, fact = NA){ # X : n x T matrix
  
  N = nrow(X); Te = ncol(X)
  
  if(anyNA(fact)){ #If fact is not given
    Xf = scale(t(X), scale=FALSE)
    svd = svd(Xf)
    lamb = t(as.matrix(sqrt(N)*t(svd$v)[1:r,]))
    fact = (Xf %*% lamb)/N
  }
  
  lm = lm(y[(1+h):Te] ~ cbind(fact[1:(Te-h),], y[1:(Te-h)]))
  est = t(coef(lm)) %*% c(1, fact[Te,], y[Te])
  
  return(est)
}

# Estimate number of quantile factors using IC
get_icr_q<-function(X, range, tau){    
  X = scale(X, scale=F)
  N = ncol(X); Te = nrow(X)
  icr = data.frame(matrix(nrow=0, ncol=2))
  for(r in range){
    qfm = do_qfm(X, r.vec=r, tau.vec=tau)  
    lambda = qfm[[1]]$L #N*r
    Factor = qfm[[1]]$F #T*r
    ic = log(loss.fun(tau, as.vector(X-Factor%*%t(lambda))))+r*(N+Te)/(N*Te)*log(N*Te/(N+Te)) #icp1
    icr = rbind(icr, c(r,ic))
  }
  colnames(icr) = c("r","ic")
  return(icr)
}

# Proposed method (forecast each quantile level)
qfm_frcst <- function(y, X, h, tau.vec, r.vec, qfm = NA){ 
  Te <- ncol(X)
  m <- length(tau.vec)
  est <- c()
  if(anyNA(qfm)){
    qfm = do_qfm( scale(t(X), scale=FALSE), r.vec=r.vec, tau.vec=tau.vec )
  }
  for(l in 1:m){
    rq = rq(y[(1+h):Te] ~ cbind(qfm[[l]]$F[1:(Te-h),], y[1:(Te-h)]), tau = tau.vec[l])
    est[l] = t(coef(rq)) %*% c(1, qfm[[l]]$F[Te,], y[Te])
  }
  return(est)
}

# Extension of proposed method (forecast each quantile level)
ext_qfm_frcst <- function(y, X, h, tau.vec, r.vec, qfm = NA){
  Te = ncol(X)
  m = length(tau.vec)
  est <- c()
  if(anyNA(qfm)){
    qfm = do_qfm( scale(t(X), scale=FALSE), r.vec=r.vec, tau.vec=tau.vec )
  }
  x = y[1:Te]
  for(l in 1:m){
    x = cbind(x, qfm[[l]]$F[1:Te,])
  }
  for(l in 1:m){ # estimate each quantile
    g = c(1,rep(2:(m+1), r.vec))
    fit.cv = cv.hrq_glasso(x=x[1:(Te-h),], y=y[(h+1):Te], group.index=g, method="quantile", tau=tau.vec[l], loss="check")
    lambda = fit.cv$lambda.min
    fit = hrq_glasso(x=x[1:(Te-h),], y=y[(h+1):Te], group.index=g, lambda = lambda, method="quantile", tau=tau.vec[l])
    est[l] = t(fit$beta)%*%c(1, x[Te,])
  } 
  return(est)
}

# Compute MAE
mae <- function(a, b){
  mean(abs(a-b))
}

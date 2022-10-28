## For Covariance Observation (For Plane)
## Note that: 
## 1. For one dimensional signal, we can set it as a speial case of two dimensional signal
## 2. For k dimensional signal, we can modify the setting to distance based.

##=====================================
## Type 1: Multiple Observation
## This works for when there are multiple observations on each location
Cov.Est.Multi <- function(X,point,
                          corrmodel = "exponential"){
  n <- dim(X)[1]
  m <- dim(X)[2]
  X.demean <- apply(X,2,function(x){x-mean(x)})
  est.ind <- 1:m
  point.for.cov <- point[est.ind,]
  
  fit <- FitComposite(X.demean[1:(n-1),est.ind], coordx=point.for.cov,
                      corrmodel=corrmodel, likelihood="Full",
                      fixed = list(mean=0), 
                      type="Standard", replicate = n-1)
  
  param <- as.list(fit$param)
  
  cov.est <- Covmatrix(point[,1], point[,2], 
                       corrmodel=corrmodel,
                       param=param)
  Sigma.eps.est <- cov.est$covmatrix * n/(n-1)
  return(list(Sigma.eps.est=Sigma.eps.est,
              fit = fit))
}


## Type 2: Single Observation
## This works for when there is only one observations on each location.
## If use the spline basis as covariate, we require the smoothing signals.
## Initial value description:(Waiting to update, depends on its performance)
Cov.Est.Covariate <- function(X,point,
                              covariate=NULL,
                              corrmodel = "exponential"){
  if(is.vector(X)){
    X <- matrix(X,ncol=length(X))
  }
  if(!is.matrix(covariate)){
    stop("Please input covariate in matrix form.")
  }
  if(nrow(covariate)!=ncol(X)|ncol(X)!=nrow(point)){
    stop("The dimension doesn't match.")
  }
  m <- dim(X)[2]
  est.ind <- 1:m
  point.for.cov <- point[est.ind,]
  # covariate <- cbind(ns(point[,1],6),ns(point[,2],6))
  
  #========================================================================
  # The common start point of sill and nugget
  rem <- c(0.01,1,0)
  names(rem) <- c("scale","sill","nugget")
  
  #========================================================================
  # We consider two types of intial value for the coefficient of covariate
  # Zero Coefficient Start Point
  param.beta.1 <- c(mean(X),rep(0,ncol(covariate)))
  names(param.beta.1) <- c("mean",paste0("beta",1:ncol(covariate)))
  param.start.1 <- as.list(c(param.beta.1,rem))
  
  # LSE Start Point
  fit.lm <- lm(as.vector(X)~covariate)
  param.beta.2 <- fit.lm$coefficients
  names(param.beta.2) <- c("mean",paste0("beta",1:ncol(covariate)))
  param.start.2 <- as.list(c(param.beta.2,rem))
  #========================================================================
  
  #========================================================================
  # Optimization with two types of start point
  #X <- matrix(X[,est.ind],nrow=1)
  t1 <- Sys.time()
  t1
  fit1 <- FitComposite2(X[,est.ind], coordx=point.for.cov, covariate = covariate,
                        corrmodel=corrmodel, likelihood="Full",#fixed = list(mean=0), 
                        start = param.start.1,
                        type="Standard", replicate = n,optimizer = "L-BFGS-B")
  t2 <- Sys.time()
  t2
  fit2 <- FitComposite2(X[,est.ind], coordx=point.for.cov, covariate = covariate,
                        corrmodel=corrmodel, likelihood="Full",#fixed = list(mean=0), 
                        start = param.start.2,
                        type="Standard", replicate = n,optimizer = "L-BFGS-B")
  t3 <- Sys.time()
  param1 <- as.list(fit1$param)
  param1 <- param1[c("sill","nugget","scale")]
  cov.est1 <- Covmatrix(point[,1], point[,2], 
                       corrmodel=corrmodel,
                       param=param1)
  Sigma.eps.est <- cov.est1$covmatrix # * n/(n-1)
  return(Sigma.eps.est)
}


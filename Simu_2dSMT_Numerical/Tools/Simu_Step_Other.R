
#==========================
# Simu Step for covariance
#==========================
one_step_Cov <- function(h, detect.m = "rad.h", seed = 0,
                         data.type = "mvnorm",
                         mu,
                         r = 0.9,
                         const = 1,
                         estcov = F,
                         ...){
  set.seed(seed)
  print(paste("seed",seed))
  # Current DataType
  # "mvnorm", "Circle","SquareSpline","CircleMix"
  n <- 3
  pluszero <- mu>0
  X <- MASS::mvrnorm(n = n, mu = mu, 
                     Sigma = Sigma.eps.p)
  X <- matrix(X, nrow = n)
  #==== Estimate covariance
  if(estcov == T){
    X.demean <- apply(X,2,function(x){x-mean(x)})
    Sigma.eps.est <- NA
    control <- 1
    #varig <- EVariogram(X.demean, point,replicates = n,cloud=T)
    #est.ind <- sample(1:m,m/3)
    est.ind <- 1:m
    
    # rep(1,length(est.ind)) is to adjust for one_dimensional
    point.for.cov <- cbind(point[est.ind,],rep(1,length(est.ind)))
    corrmodel <- "exponential"
    
    print(paste0(c("seed",seed,":Start Full....")))
    time.flag.1.1 <- Sys.time()
    fit <- FitComposite(X.demean[1:(n-1),est.ind], coordx=point.for.cov,
                         corrmodel=corrmodel, likelihood="Full",fixed = list(mean=0), 
                         type="Standard", replicate = n-1)
    fit1 <- fit
    time.flag.1.2 <- Sys.time()
    print(paste0(c("seed",seed,":End Full.")))
    
    param <- as.list(fit$param)
    
    cov.est <- Covmatrix(cbind(point,rep(1,length(point))), 
                         corrmodel=corrmodel,
                         param=param)
    Sigma.eps.est.full <- cov.est$covmatrix * n/(n-1)
    
    print(paste0(c("seed",seed,":Start Taper....")))
    time.flag.2.1 <- Sys.time()
    fit <- FitComposite(X.demean[1:(n-1),est.ind], 
                        coordx=point.for.cov,
                        corrmodel=corrmodel,likelihood="Full",
                        type="Tapering",taper="Wendland1",maxdist=10,
                        start = list(scale=0.5,sill=0.5,nugget = 1),
                        fixed = list(mean=0), 
                        replicate = n-1)
    fit2 <- fit
    time.flag.2.2 <- Sys.time()
    print(paste0(c("seed",seed,":End Taper.")))
    param <- as.list(fit$param)
    
    cov.est <- Covmatrix(cbind(point,rep(1,length(point))), 
                         corrmodel=corrmodel,
                         param=param)
    Sigma.eps.est.taper <- cov.est$covmatrix * n/(n-1)
    print(Sigma.eps.est.full[1:5,1:5])
    print(Sigma.eps.est.taper[1:5,1:5])
    save(X, point.for.cov,fit1,fit2,
         file=paste0("Result/SimulationCov/Tmp/rep_mag_",magnitude,"_mu_",mu_type,"_cov_",Cov_type,"_seed_",seed,".RData"))
    
  }
  return(list(Sigma.full = Sigma.eps.est.full[1,],
              Sigma.taper = Sigma.eps.est.taper[1,],
              time.full = difftime(time.flag.1.1,time.flag.1.2,units = "secs"),
              time.taper = difftime(time.flag.2.1,time.flag.2.2,units = "secs")))
}

#==========================
# Simu Step for covariance (General Vecchia)
#==========================
one_step_Cov_Vecchia <- function(h, detect.m = "rad.h", seed = 0,
                         data.type = "mvnorm",
                         mu,
                         r = 0.9,
                         const = 1,
                         estcov = F,
                         ...){
  set.seed(seed)
  print(paste("seed",seed))
  # Current DataType
  # "mvnorm", "Circle","SquareSpline","CircleMix"
  n <- 3
  pluszero <- mu>0
  X <- MASS::mvrnorm(n = n, mu = mu, 
                     Sigma = Sigma.eps.p)
  X <- matrix(X, nrow = n)
  #==== Estimate covariance
  if(estcov == T){
    X.demean <- apply(X,2,function(x){x-mean(x)})
    covparms <- c(1.0, 0.1,0.5)# wish larger covariance
    
    vecchia.est=vecchia_estimate(X.demean[1:(n-1),],point,theta.ini=c(covparms,1))
    cov.est <- MaternFun(Dist.p,vecchia.est$theta.hat) * n/(n-1)
    cov.test<- MaternFun(Dist.p,c(covparms,0.1)) * n/(n-1)
    
    print(paste0(c("seed",seed,":Start Taper....")))
    time.flag.2.1 <- Sys.time()
    fit <- FitComposite(X.demean[1:(n-1),est.ind], 
                        coordx=point.for.cov,
                        corrmodel=corrmodel,likelihood="Full",
                        type="Tapering",taper="Wendland1",maxdist=10,
                        start = list(scale=0.5,sill=0.5,nugget = 1),
                        fixed = list(mean=0), 
                        replicate = n-1)
    fit2 <- fit
    time.flag.2.2 <- Sys.time()
    print(paste0(c("seed",seed,":End Taper.")))
    param <- as.list(fit$param)
    
    cov.est <- Covmatrix(cbind(point,rep(1,length(point))), 
                         corrmodel=corrmodel,
                         param=param)
    Sigma.eps.est.taper <- cov.est$covmatrix * n/(n-1)
    print(Sigma.eps.est.full[1:5,1:5])
    print(Sigma.eps.est.taper[1:5,1:5])
    save(X, point.for.cov,fit1,fit2,
         file=paste0("Result/SimulationCov/Tmp/rep_mag_",magnitude,"_mu_",mu_type,"_cov_",Cov_type,"_seed_",seed,".RData"))
    
  }
  return(list(Sigma.full = Sigma.eps.est.full[1,],
              Sigma.taper = Sigma.eps.est.taper[1,],
              time.full = difftime(time.flag.1.1,time.flag.1.2,units = "secs"),
              time.taper = difftime(time.flag.2.1,time.flag.2.2,units = "secs")))
}


#==========================
# Simu Step for covariance
# One Sample: Not directly worked!
#==========================
one_step_Cov_one_sample <- function(h, detect.m = "rad.h", seed = 0,
                                    data.type = "mvnorm",
                                    magnitude = 4,
                                    r = 0.9,
                                    const = 1,
                                    estcov = F,
                                    ...){
  set.seed(seed)
  print(paste("seed",seed))
  # Current DataType
  # "mvnorm", "Circle","SquareSpline","CircleMix"
  if(data.type == "mvnorm"){
    Sigma.mu.p <- Sigma.mu(m, rho.mu)
    mu.res <- mu.fix.gen.1D(point, "mvnorm", mu.mean = mu.mean, Sigma.mu = Sigma.mu.p, 
                            magnitude = magnitude)
  }else if(data.type == "Circle"){
    mu.res <- mu.fix.gen.1D(point, "uc.unif", magnitude)
  }else if(data.type == "SquareSpline"){
    mu.res <- mu.fix.gen.1D(point, "uc.spline", magnitude)
  }else if(data.type == "CircleMix"){
    mu.res <- mu.fix.gen.1D(point, "mixture", magnitude)
  }
  
  mu <- mu.res$mu
  pis.hat3 <- mu.res$pis
  
  pluszero <- mu>0
  X <- MASS::mvrnorm(n = n, mu = mu, 
                     Sigma = Sigma.eps.p)
  X <- matrix(X, nrow = n)
  #==== Estimate covariance
  if(estcov == T){
    #X.demean <- apply(X,2,function(x){x-mean(x)})
    Sigma.eps.est <- NA
    
    control <- 1
    #varig <- EVariogram(X.demean, point,replicates = n,cloud=T)
    #est.ind <- sample(1:m,m/3)
    est.ind <- 1:m
    
    # rep(1,length(est.ind)) is to adjust for one_dimensional
    point.for.cov <- cbind(point[est.ind,],rep(1,length(est.ind)))
    corrmodel <- "exponential"
    
    print(paste0(c("seed",seed,":Start Full....")))
    time.flag.1.1 <- Sys.time()
    fit <- FitComposite(X[1:(n),est.ind], coordx=point.for.cov,
                        corrmodel=corrmodel, likelihood="Full",#fixed = list(mean=0), 
                        type="Standard", replicate = n)
    time.flag.1.2 <- Sys.time()
    print(paste0(c("seed",seed,":End Full.")))
    
    param <- as.list(fit$param)
    
    cov.est <- Covmatrix(cbind(point,rep(1,length(point))), 
                         corrmodel=corrmodel,
                         param=param)
    Sigma.eps.est.full <- cov.est$covmatrix #* n/(n-1)
    
    print(paste0(c("seed",seed,":Start Taper....")))
    time.flag.2.1 <- Sys.time()
    fit <- FitComposite(X[1:n,est.ind], 
                        coordx=point.for.cov,
                        corrmodel=corrmodel,likelihood="Full",
                        type="Tapering",taper="Wendland1",maxdist=10,
                        start = list(scale=0.5,sill=0.5,nugget = 1),
                        #fixed = list(mean=0), 
                        replicate = n)
    
    time.flag.2.2 <- Sys.time()
    print(paste0(c("seed",seed,":End Taper.")))
    param <- as.list(fit$param)
    
    cov.est <- Covmatrix(cbind(point,rep(1,length(point))), 
                         corrmodel=corrmodel,
                         param=param)
    Sigma.eps.est.taper <- cov.est$covmatrix
    print(Sigma.eps.est.full[1:5,1:5])
    print(Sigma.eps.est.taper[1:5,1:5])
  }
  return(list(Sigma.full = Sigma.eps.est.full[1,],
              Sigma.taper = Sigma.eps.est.taper[1,],
              time.full = difftime(time.flag.1.1,time.flag.1.2,units = "secs"),
              time.taper = difftime(time.flag.2.1,time.flag.2.2,units = "secs")))
}

#==========================
# Simu Step for covariance
# One Sample: Adjusted by covariate
#==========================
one_step_Cov_covariate <- function(h, detect.m = "rad.h", seed = 0,
                                   data.type = "mvnorm",
                                   mu,
                                   r = 0.9,
                                   const = 1,
                                   estcov = F,
                                   ...){
  set.seed(seed)
  print(paste("seed",seed))
  # Current DataType
  # "mvnorm", "Circle","SquareSpline","CircleMix"
  
  pluszero <- mu>0
  X <- MASS::mvrnorm(n = n, mu = mu, 
                     Sigma = Sigma.eps.p)
  X <- matrix(X, nrow = n)
  #==== Estimate covariance
  if(estcov == T){
    #X.demean <- apply(X,2,function(x){x-mean(x)})
    Sigma.eps.est <- NA
    
    control <- 1
    #varig <- EVariogram(X.demean, point,replicates = n,cloud=T)
    #est.ind <- sample(1:m,m/3)
    est.ind <- 1:m
    
    # rep(1,length(est.ind)) is to adjust for one_dimensional
    point.for.cov <- cbind(point[est.ind,],rep(1,length(est.ind)))
    corrmodel <- "exponential"
    
    print(paste0(c("seed",seed,":Start Full....")))
    covariate = ns(point.for.cov[,1],6)
    time.flag.1.1 <- Sys.time()
    fit1 <- FitComposite2(X[1:(n),est.ind], coordx=point.for.cov, covariate = covariate,
                          corrmodel=corrmodel, likelihood="Full",#fixed = list(mean=0), 
                          start = list(beta1=0,beta2=0,beta3=0,beta4=0,beta5=0,beta6=0),
                          type="Standard", replicate = n,optimizer = "L-BFGS-B")
    time.flag.1.2 <- Sys.time()
    
    fit.lm <- lm(as.vector(X)~covariate-1)
    start2.pre <- fit.lm$coefficients
    names(start2.pre) <- paste0("beta",1:6)
    start2 <- as.list(start2.pre)
    fit2 <- FitComposite2(X[1:(n),est.ind], coordx=point.for.cov, covariate = covariate,
                          corrmodel=corrmodel, likelihood="Full",#fixed = list(mean=0), 
                          start = start2,
                          type="Standard", replicate = n,optimizer = "L-BFGS-B")
    
    
    time.flag.3.1 <- Sys.time()
    # Start with true parameter
    fit3 <- FitComposite2(X[1:(n),est.ind], coordx=point.for.cov, covariate = covariate,
                          corrmodel=corrmodel, likelihood="Full",#fixed = list(mean=0), 
                          start = list(beta1=0,beta2=0,beta3=0,beta4=0,beta5=0,beta6=0,
                                       scale=0.1,sill=0.5,nugget=0.5),
                          type="Standard", replicate = n,optimizer = "L-BFGS-B")
    time.flag.3.2 <- Sys.time()
    
    rem <- c(0.1,0.5,0.5);names(rem) <- c("scale","sill","nugget")
    rem <- c(2000,1,0);names(rem) <- c("scale","sill","nugget")
    rem <- c(1,0);names(rem) <- c("sill","nugget")
    start4 <- list(c(start2.pre,rem))
    fit4 <- FitComposite2(X[1:(n),est.ind], coordx=point.for.cov, covariate = covariate,
                          corrmodel=corrmodel, likelihood="Full",#fixed = list(mean=0), 
                          start = start4,
                          type="Standard", replicate = n,optimizer = "L-BFGS-B")
    
    if((fit1$logCompLik>fit2$logCompLik)&(fit1$logCompLik>fit3$logCompLik)&(fit1$logCompLik>fit4$logCompLik)){
      fit = fit1
    }else if((fit2$logCompLik>fit3$logCompLik)&(fit2$logCompLik>fit3$logCompLik)){
      fit = fit2
    }else if(fit3$logCompLik>fit4$logCompLik){
      fit = fit3
    }else{
      fit = fit4
    }
    fit = fit4
    print(paste0(c("seed",seed,":End Full.")))
    time.flag.2.1 <- Sys.time()
    param <- as.list(fit$param)
    paramcorr <- param[!strtrim(names(param),4)=="beta"]
    
    cov.est <- Covmatrix(cbind(point,rep(1,length(point))), 
                         corrmodel=corrmodel,
                         param=paramcorr)
    Sigma.eps.est.full <- cov.est$covmatrix #* n/(n-1)
    
    save(X, point.for.cov,fit,fit1,fit2,fit3,fit4,
         file=paste0("Result/SimulationCov/Tmp1/no_rep_mag_",magnitude,"_mu_",mu_type,"_cov_",Cov_type,"_seed_",seed,".RData"))
    print(paste0(c("seed",seed,":Start Taper....")))
    if(F){
      fit <- FitComposite(X[1:n,est.ind], 
                          coordx=point.for.cov,
                          corrmodel=corrmodel,likelihood="Full",
                          type="Tapering",taper="Wendland1",maxdist=10,
                          start = list(scale=0.5,sill=0.5,nugget = 1),
                          #fixed = list(mean=0), 
                          replicate = n)
    }
    time.flag.2.2 <- Sys.time()
    print(paste0(c("seed",seed,":End Taper.")))
    param <- as.list(fit$param)
    
    paramcorr <- param[!strtrim(names(param),4)=="beta"]
    
    cov.est <- Covmatrix(cbind(point,rep(1,length(point))), 
                         corrmodel=corrmodel,
                         #param=param,
                         param=paramcorr)
    
    Sigma.eps.est.taper <- cov.est$covmatrix
    print(Sigma.eps.est.full[1:5,1:5])
    print(Sigma.eps.est.taper[1:5,1:5])
    print(Sigma.eps.p[1:5,1:5])
  }
  return(list(Sigma.full = Sigma.eps.est.full[1,],
              Sigma.taper = Sigma.eps.est.taper[1,],
              time.full = difftime(time.flag.1.1,time.flag.1.2,units = "secs"),
              time.taper = difftime(time.flag.2.1,time.flag.2.2,units = "secs")))
}


#==========================
# Simu Step for pis
#==========================
one_step_pis <- function(h, detect.m = "rad.h", seed = 0,
                         data.type = "mvnorm",
                         magnitude = 4,
                         r = 0.9,
                         const = 1,
                         estcov = F,
                         ...){
  set.seed(seed)
  print(paste("seed",seed))
  # Current DataType
  # "mvnorm", "Circle","SquareSpline","CircleMix"
  if(data.type == "mvnorm"){
    Sigma.mu.p <- Sigma.mu(m, rho.mu)
    mu.res <- mu.fix.gen.1D(point, "mvnorm", mu.mean = mu.mean, Sigma.mu = Sigma.mu.p, 
                            magnitude = magnitude)
  }else if(data.type == "Circle"){
    mu.res <- mu.fix.gen.1D(point, "uc.unif", magnitude)
  }else if(data.type == "SquareSpline"){
    mu.res <- mu.fix.gen.1D(point, "uc.spline", magnitude)
  }else if(data.type == "CircleMix"){
    mu.res <- mu.fix.gen.1D(point, "mixture", magnitude)
  }
  
  mu <- mu.res$mu
  pis.hat3 <- mu.res$pis
  
  pluszero <- mu>0
  X <- MASS::mvrnorm(n = n, mu = mu, 
                     Sigma = Sigma.eps.p)
  X <- matrix(X, nrow = n)
  
  Sigma.eps.est <- Sigma.eps.p
  sgm <- sqrt(Sigma.eps.est[1,1])
  #==== Perform Algorithm
  dig <- 7
  T1 <- apply(X,2,function(x){sum(x)/sgm/sqrt(n)})
  T1 <- round(T1,dig)
  p.value <- 1 - pnorm(T1)
  
  #==== Calculate pis
  #==== In LAWS: kernel regression
  cat("seed",seed,"Start LAWS ==========\n")
  bh.th<-bh.func(p.value, 0.9)$th
  if(bh.th == 1){
    pis.hat <- rep(1,m)
  }else{
    pis.hat<-pis_1D.func.ker.reg(T1, point, tau=bh.th)
  }
  pis0.LAWS <- 1-pis.hat
  cat("seed",seed,"End LAWS ==========\n")
  
  #==== Calculate qhat in SABHA
  #==== We use the parameter in Ring and Li
  cat("seed",seed,"Start SABHA ==========\n")
  tau = 0.5;eps = 0.1; TV_bd = 2
  alpha_ADMM = 10^2; beta = 10^3; eta = 5; max_iters = 5000; converge_thr = 1e-4 # parameters for ADMM
  ADMM_params = c(alpha_ADMM,beta,eta,max_iters,converge_thr)
  
  qhat = Solve_q_TV_1dim(p.value, tau, eps, TV_bd, ADMM_params)
  pis0.SABHA <- qhat
  cat("seed",seed,"End SABHA ==========\n")
  
  #===========================================
  #==== AdaMT
  #==== Lei and Fithian 
  #==== Doesn't estimate pis. It adaptively updates the thresholds
  #===========================================
  #formula <- paste0("ns(x1, df=6)")
  #res_gam <- adapt_gam(x = data.frame(x1=point[,1]), pvals = p.value, 
  #                     pi_formulas = formula, mu_formulas = formula,alphas = q)
  #res_gam
  #===========================================
  #==== CAMT
  #==== Zhang and Chen
  #===========================================
  cat("seed",seed,"Start CAMT ==========\n")
  x1 = point[,1]
  CAMT.res <- camt.fdr(pvals = p.value, pi0.var = ns(x1,6), f1.var = ns(x1,6), 
                       alg.type = "EM", control.method = "knockoff+")
  pis0.CAMT <- CAMT.res$pi0
  #plot(c(1:900),CAMT.res$pi0)
  cat("seed",seed,"End CAMT ==========\n")
  
  #===========================================
  #==== FDRreg
  #==== Scott etc.
  #==== The prob is for non-null
  #===========================================
  cat("seed",seed,"Start FDRregT ==========\n")
  FDRregT.res <- FDRreg(T1, covars=ns(x1,6),nulltype= 'theoretical')
  cat("seed",seed,"End FDRregT ==========\n")
  pis0.FDRregT <- 1-FDRregT.res$priorprob
  if(F){
    cat("seed",seed,"Start FDRregE ==========\n")
    FDRregE.res <- FDRreg(T1, covars=ns(x1,6),nulltype= 'empirical')
    cat("seed",seed,"End FDRregE ==========\n")
    pis0.FDRregE <- 1-FDRregE.res$priorprob
  }
  return(list(pis0.LAWS = pis0.LAWS,
              pis0.SABHA = pis0.SABHA,
              pis0.CAMT = pis0.CAMT,
              pis0.FDRregT = pis0.FDRregT,
              pis0.True = mu.res$pis
              #pis0.FDRregE = pis0.FDRregE
  ))
}

#==========================
# Simu Step for covariance
# One Sample: Adjusted by covariate
#==========================
one_step_Cov_covariate_2D <- function(h, detect.m = "rad.h", seed = 0,
                                   data.type = "mvnorm",
                                   mu,
                                   r = 0.9,
                                   const = 1,
                                   estcov = F,
                                   ...){
  set.seed(seed)
  print(paste("seed",seed))
  # Current DataType
  # "mvnorm", "Circle","SquareSpline","CircleMix"
  
  pluszero <- mu>0
  X <- MASS::mvrnorm(n = n, mu = mu, 
                     Sigma = Sigma.eps.p)
  X <- matrix(X, nrow = n)
  
  covariate <- cbind(ns(point[,1],4),ns(point[,2],4))
  covar1 <- ns(point[,1],4)
  covar2 <- ns(point[,2],4)
  covariate <- matrix(NA,nrow=nrow(point),ncol = 16)
  #cbind(ns(point[,1],4),ns(point[,2],4))
  for(i in 1:4){
    for(j in 1:4){
      covariate[,(i-1)*4+j] <- covar1[,i]*covar2[,j]
    }
  }
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
  rem <- c(0.2,1,0)
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
  time.flag.1.1 <- Sys.time()
  
  if(F){
  fit1 <- FitComposite2(X[,est.ind], coordx=point.for.cov, covariate = covariate,
                        corrmodel=corrmodel, likelihood="Full",#fixed = list(mean=0), 
                        start = param.start.1,
                        type="Standard", replicate = n,optimizer = "L-BFGS-B")
  fit2 <- FitComposite2(X[,est.ind], coordx=point.for.cov, covariate = covariate,
                        corrmodel=corrmodel, likelihood="Full",#fixed = list(mean=0), 
                        start = param.start.2,
                        type="Standard", replicate = n,optimizer = "L-BFGS-B")
  
    
  }
  
  rem <- c(0.1,0.5,0.5)
  names(rem) <- c("scale","sill","nugget")
  corrmodel1 <- "gauss"
  fit1 <- FitComposite2(X[,est.ind], coordx=point.for.cov, covariate = covariate,
                        corrmodel=corrmodel1, likelihood="Full",#fixed = list(mean=0), 
                        start = as.list(rem),
                        type="Standard", replicate = n,optimizer = "L-BFGS-B")
  time.flag.1.2 <- Sys.time()
  
  rem <- c(0.1,0.8,0.2)
  names(rem) <- c("scale","sill","nugget")
  corrmodel2 <- "exponential"
  fit2 <- FitComposite2(X[,est.ind], coordx=point.for.cov, covariate = covariate,
                        corrmodel=corrmodel2, likelihood="Full",#fixed = list(mean=0), 
                        start = as.list(rem),
                        type="Standard", replicate = n,optimizer = "L-BFGS-B")
  
  time.flag.1.3 <- Sys.time()
  param1 <- as.list(fit1$param)
  param1 <- param1[c("sill","nugget","scale")]
  cov.est1 <- Covmatrix(point[,1], point[,2], 
                        corrmodel=corrmodel1,
                        param=param1)
  Sigma.eps.est.1 <- cov.est1$covmatrix # * n/(n-1)
  
  
  param1 <- as.list(fit2$param)
  param1 <- param1[c("sill","nugget","scale")]
  cov.est1 <- Covmatrix(point[,1], point[,2], 
                        corrmodel=corrmodel2,
                        param=param1)
  Sigma.eps.est.2 <- cov.est1$covmatrix # * n/(n-1)
  save.folder.name <- "SimulationCov"
  save.file.path <- paste0("Result/",save.folder.name,"/Tmp/mag_",
                           magnitude,"_mu_",mu_type,"_cov_",Cov_type,"_seed_",seed,".RData")
  
  save(list = ls(all.names = TRUE), file =save.file.path)
  return(list(Sigma.full = Sigma.eps.est.1[1,],
              Sigma.taper = Sigma.eps.est.2[1,],
              time.full = difftime(time.flag.1.1,time.flag.1.2,units = "secs"),
              time.taper = difftime(time.flag.1.2,time.flag.1.3,units = "secs")))
}
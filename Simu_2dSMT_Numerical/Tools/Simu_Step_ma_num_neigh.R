
#==================
# Simulation Step
# 1 dimension
# To evaluate the performance of adaptively choose neighorhood
#==================
one_step_1D_num_neigh <- function(h, grp,detect.m = "rad.h", seed = 0,
                        data.type = "mvnorm",
                        mu = 2,
                        const = 1,
                        ws.fun = function(x){1/x},
                        nbins.ihw = 5,
                        n.group.max = 5,
                        tau.tm = 1,
                        tau.ta = 1,
                        magnitude = 0,
                        mu_type = "Sparse",
                        Cov_type = "Weak",
                        save.folder.name = NULL,
                        EmpMethod = "NPEB",
                        ...){
  set.seed(seed)
  print(paste("seed",seed))
  # Current DataType
  # "mvnorm", "Circle","SquareSpline","CircleMix"
  #  if(data.type == "mvnorm"){
  #    Sigma.mu.p <- Sigma.mu(m, rho.mu)
  #    mu.res <- mu.fix.gen.1D(point, "mvnorm", mu.mean = mu.mean, Sigma.mu = Sigma.mu.p,
  #                            magnitude = magnitude)
  #  }else if(data.type == "Circle"){
  #    mu.res <- mu.fix.gen.1D(point, "uc.unif", magnitude)
  #  }else if(data.type == "SquareSpline"){
  #    mu.res <- mu.fix.gen.1D(point, "uc.spline", magnitude)
  #  }else if(data.type == "CircleMix"){
  #    mu.res <- mu.fix.gen.1D(point, "mixture", magnitude)
  #  }
  #save.folder.name <- "Simulation1D_unif"
  #save.folder.name <- "Simulation1D"
  #save.folder.name <- "Simulation1D_mv"
  save.file.path <- paste0("Result/",save.folder.name,"/Tmp/mag_",
                           magnitude,"_mu_",mu_type,"_cov_",Cov_type,"_seed_",seed,".RData")
  if(!file.exists(save.file.path)){
    pluszero <- mu>0
    #X <- user_mvrnorm(n = n, mu = mu,
    #                  Sigma = Sigma.eps.p)
    X <- MASS::mvrnorm(n = n +1, mu = mu,
                       Sigma = Sigma.eps.p)
    X.neigh <- matrix(X[n+1,], nrow = n)
    X <- matrix(X[1:n,], nrow = n)
    
    #======================================
    #==== Estimate covariance
    #======================================
    if(estcov==T){
      X.demean <- apply(X,2,function(x){x-mean(x)})
      Sigma.eps.est <- NA
      #est.ind <- sample(1:m,m/3)
      est.ind <- 1:m
      
      # Fit Composite only works for 2 dimensional spatial domain
      point.for.cov <- cbind(point[est.ind,],rep(1,length(est.ind)))
      corrmodel <- "exponential"
      
      #fit <- FitComposite(X.demean[1:(n-1),est.ind], coordx=point.for.cov,
      #                    corrmodel=corrmodel, likelihood="Full",fixed = list(mean=0),
      #                    type="Standard", replicate = n-1)
      
      # Tapering is more accurate, so we choose Tapering
      fit <- FitComposite(X.demean[1:(n-1),est.ind],
                          coordx=point.for.cov,
                          corrmodel=corrmodel,likelihood="Full",
                          type="Tapering",taper="Wendland1",maxdist=10,
                          start = list(scale=0.5,sill=0.5,nugget = 1),
                          fixed = list(mean=0),
                          replicate = n-1)
      param <- as.list(fit$param)
      
      # Generate the covariance matrix
      cov.est <- Covmatrix(cbind(point,rep(1,length(point))),
                           corrmodel=corrmodel,
                           param=param)
      Sigma.eps.est <- cov.est$covmatrix * n/(n-1)
      
      param <- as.list(fit$param)
      
      cov.est <- Covmatrix(cbind(point,rep(1,length(point))),
                           corrmodel=corrmodel,
                           param=param)
      Sigma.eps.est <- cov.est$covmatrix * n/(n-1)
    }else{
      Sigma.eps.est <- Sigma.eps.p
    }
    
    #==== End estimating covariance
    #=====================================
    
    sgm <- sqrt(diag(Sigma.eps.est))
    #==== Perform Algorithm
    dig <- 7
    Tm <- apply(X,2,function(x){sum(x)/sqrt(n)})/sgm
    Tm <- round(Tm,dig)
    
    fdp.res <-c()
    pow.res <-c()
    NP.res <- c()
    
    #==== Calculate pis
    print("Start calculating pis....")
    #==== Calculate pis2 (globalwise)
    lambda.test <- 0.5#median(Tm)#quantile(Tm,0.9)#
    pis.hata.tmp<- min(1,
                       mean(Tm<lambda.test)/pnorm(lambda.test))
    pis.hata <- rep(pis.hata.tmp,m)
    #pis.hata <- (Tm<lambda.test)/pnorm(lambda.test)
    
    
    #================================
    #==== Without spatial Info: One D
    #================================
    res.1D <- OneD_Detect(Tm, q,  pis = NULL,
                          const = const,
                          tau.tm = tau.tm)
    Tm.star <- res.1D$tm.min
    max.rej <- res.1D$max.rej
    fdp.res <- c(fdp.res,fdp(res.1D$selected,mu))
    pow.res <- c(pow.res,Pow(res.1D$selected,mu))
    NP.res <- c(NP.res,length(res.1D$selected))
    
    #===========================================
    #==== Without spatial Info: One D (with pis2)
    #===========================================
    res.1D.pis2 <- OneD_Detect(Tm, q, pis = pis.hata,
                               const = const,
                               tau.tm = tau.tm)
    Tm.star.pis2 <- res.1D.pis2$tm.min
    max.rej.pis2 <- res.1D.pis2$max.rej
    fdp.res <- c(fdp.res,fdp(res.1D.pis2$selected,mu))
    pow.res <- c(pow.res,Pow(res.1D.pis2$selected,mu))
    NP.res <- c(NP.res,length(res.1D.pis2$selected))
    
    
    fdp.res <- c(fdp(res.1D$selected,mu),
                 fdp(res.1D.pis2$selected,mu))
    pow.res <- c(Pow(res.1D$selected,mu),
                 Pow(res.1D.pis2$selected,mu))
    #=== Rename
    names(fdp.res) <- c("1D", "1D.pis2")
    #=== With Spatial Info
    #for(ii in 1:length(h)){
    #  hh <- h[ii]
    #---- Save Result
    #save.file.path <- paste0("Result/Simulation1D/Tmp/mag_",
    #                         magnitude,"_mu_",mu_type,"_cov_",Cov_type,"_seed_",seed,".RData")
    #save(list = ls(all.names = TRUE), file =save.file.path)
  }else{
    # Load the previous result and save time.
    load(save.file.path)
  }
  
  
    print(paste("seed:",seed)) 
    
    Tm.star <- round(Tm.star,dig)
    Tm.star.pis2 <- round(Tm.star.pis2,dig)
    
    Ta.star <- Inf
    Ta.star.pis2 <- Inf
    #===========================================
    #====== Fix Reigion
    #===========================================
    
    # Detect Neighbor
    # We consider three ways of choosing the number of neighbor
    hh <- 4
    hh.seq <- rep(4,m)
    hh.extra.seq <- hh.Sel(X.neigh)
    hh.current.seq <- hh.Sel(X)
    
    ##=========================
    ## Use hh.seq = rep(4,m)
    Neigh_Detect_res <- TwoDSMT::Neigh_Detect(hh = hh.seq,
                                              X = X,
                                              Dist = Dist.p,
                                              Sigma.eps = Sigma.eps.est,
                                              detect.m = detect.m)
    # mu = mu)
    Ta <- Neigh_Detect_res$Ta
    Va <- Neigh_Detect_res$Va
    VmVa.cov <- Neigh_Detect_res$VmVa.cov
    ind <- Neigh_Detect_res$ind
    mua <- Neigh_Detect_res$mua
    sig_ratio <- Neigh_Detect_res$sig_ratio
    
    #===========================================
    #== Run 2D Selection
    #===========================================
    print(paste("seed:",seed,"BH","hh: Type1 "))
    is.quick <- T
    res.2D <- TwoDSMT::Spatial_Detect_exact_grp_BH_down(Tm, Ta, Va,
                                                        VmVa.cov, ind,
                                                        q, max.rej,
                                                        Tm.star = Tm.star,
                                                        Ta.star = Ta.star,
                                                        const = const,
                                                        seed = seed,
                                                        mua = mua,
                                                        sig_ratio = sig_ratio,
                                                        EmpMethod = EmpMethod,
                                                        is.quick.stop = is.quick)
    selected.2D <- res.2D$selected
    tm <- res.2D$tm0
    ta <- res.2D$ta0
    print(paste("seed:",seed,"BH","hh: Type1 final.fdr",res.2D$final.fdr))
    
    print(paste("seed:",seed,"BH_pis2","hh: Type1 "))
    res.2D.pis2 <- TwoDSMT::Spatial_Detect_exact_grp_BH_down(Tm, Ta, Va,
                                                             VmVa.cov, ind,
                                                             q,
                                                             max.rej = max.rej.pis2,
                                                             pis = pis.hata,
                                                             Tm.star = Tm.star.pis2,
                                                             Ta.star = Ta.star.pis2,
                                                             const = const,
                                                             seed = seed,
                                                             mua = mua,
                                                             EmpMethod = EmpMethod,
                                                             is.quick.stop = is.quick)
    selected.2D.pis2 <- res.2D.pis2$selected
    tm <- res.2D.pis2$tm0
    ta <- res.2D.pis2$ta0
    print(paste("seed:",seed,"BH.pis2","hh: Type1 final.fdr",res.2D.pis2$final.fdr))
    
    
    
    ##=========================
    ## Use  hh.extra.seq <- hh.Sel(X.neigh)
    Neigh_Detect_res <- TwoDSMT::Neigh_Detect(hh = hh.extra.seq,
                                              X = X,
                                              Dist = Dist.p,
                                              Sigma.eps = Sigma.eps.est,
                                              detect.m = detect.m)
    # mu = mu)
    Ta <- Neigh_Detect_res$Ta
    Va <- Neigh_Detect_res$Va
    VmVa.cov <- Neigh_Detect_res$VmVa.cov
    ind <- Neigh_Detect_res$ind
    mua <- Neigh_Detect_res$mua
    sig_ratio <- Neigh_Detect_res$sig_ratio
    
    #===========================================
    #== Run 2D Selection (hh.extra.seq)
    #===========================================
    print(paste("seed:",seed,"BH","hh: Type2 "))
    is.quick <- T
    res.2D.hhType2 <- TwoDSMT::Spatial_Detect_exact_grp_BH_down(Tm, Ta, Va,
                                                        VmVa.cov, ind,
                                                        q, max.rej,
                                                        Tm.star = Tm.star,
                                                        Ta.star = Ta.star,
                                                        const = const,
                                                        seed = seed,
                                                        mua = mua,
                                                        sig_ratio = sig_ratio,
                                                        EmpMethod = EmpMethod,
                                                        is.quick.stop = is.quick)
    selected.2D.hhType2 <- res.2D.hhType2$selected
    tm <- res.2D.hhType2$tm0
    ta <- res.2D.hhType2$ta0
    print(paste("seed:",seed,"BH","hh: Type2 final.fdr",res.2D.hhType2$final.fdr))
    
    print(paste("seed:",seed,"BH_pis2","hh: Type2 "))
    res.2D.pis2.hhType2 <- TwoDSMT::Spatial_Detect_exact_grp_BH_down(Tm, Ta, Va,
                                                             VmVa.cov, ind,
                                                             q,
                                                             max.rej = max.rej.pis2,
                                                             pis = pis.hata,
                                                             Tm.star = Tm.star.pis2,
                                                             Ta.star = Ta.star.pis2,
                                                             const = const,
                                                             seed = seed,
                                                             mua = mua,
                                                             EmpMethod = EmpMethod,
                                                             is.quick.stop = is.quick)
    selected.2D.pis2.hhType2 <- res.2D.pis2.hhType2$selected
    tm <- res.2D.pis2.hhType2$tm0
    ta <- res.2D.pis2.hhType2$ta0
    print(paste("seed:",seed,"BH.pis2","hh: Type2 final.fdr",res.2D.pis2.hhType2$final.fdr))
    
    
    ##=========================
    ## Use  hh.current.seq <- hh.Sel(X)
    Neigh_Detect_res <- TwoDSMT::Neigh_Detect(hh = hh.current.seq,
                                              X = X,
                                              Dist = Dist.p,
                                              Sigma.eps = Sigma.eps.est,
                                              detect.m = detect.m)
    # mu = mu)
    Ta <- Neigh_Detect_res$Ta
    Va <- Neigh_Detect_res$Va
    VmVa.cov <- Neigh_Detect_res$VmVa.cov
    ind <- Neigh_Detect_res$ind
    mua <- Neigh_Detect_res$mua
    sig_ratio <- Neigh_Detect_res$sig_ratio
    
    #===========================================
    #== Run 2D Selection (hh.current.seq)
    #===========================================
    print(paste("seed:",seed,"BH","hh: Type3 "))
    is.quick <- T
    res.2D.hhType3 <- TwoDSMT::Spatial_Detect_exact_grp_BH_down(Tm, Ta, Va,
                                                        VmVa.cov, ind,
                                                        q, max.rej,
                                                        Tm.star = Tm.star,
                                                        Ta.star = Ta.star,
                                                        const = const,
                                                        seed = seed,
                                                        mua = mua,
                                                        sig_ratio = sig_ratio,
                                                        EmpMethod = EmpMethod,
                                                        is.quick.stop = is.quick)
    selected.2D.hhType3 <- res.2D.hhType3$selected
    tm <- res.2D.hhType3$tm0
    ta <- res.2D.hhType3$ta0
    print(paste("seed:",seed,"BH","hh: Type3 final.fdr",res.2D.hhType3$final.fdr))
    
    print(paste("seed:",seed,"BH_pis2","hh: Type3 "))
    res.2D.pis2.hhType3 <- TwoDSMT::Spatial_Detect_exact_grp_BH_down(Tm, Ta, Va,
                                                             VmVa.cov, ind,
                                                             q,
                                                             max.rej = max.rej.pis2,
                                                             pis = pis.hata,
                                                             Tm.star = Tm.star.pis2,
                                                             Ta.star = Ta.star.pis2,
                                                             const = const,
                                                             seed = seed,
                                                             mua = mua,
                                                             EmpMethod = EmpMethod,
                                                             is.quick.stop = is.quick)
    selected.2D.pis2.hhType3 <- res.2D.pis2.hhType3$selected
    tm <- res.2D.pis2.hhType3$tm0
    ta <- res.2D.pis2.hhType3$ta0
    print(paste("seed:",seed,"BH.pis2","hh: Type3 final.fdr",res.2D.pis2.hhType3$final.fdr))
    
    
    fdp.res <- c(fdp.res,
                 fdp(selected.2D,mu),
                 fdp(selected.2D.pis2,mu),
                 fdp(selected.2D.hhType2,mu),
                 fdp(selected.2D.pis2.hhType2,mu),
                 fdp(selected.2D.hhType3,mu),
                 fdp(selected.2D.pis2.hhType3,mu)
                 
    )
    
    pow.res <- c(pow.res,
                 Pow(selected.2D,mu),
                 Pow(selected.2D.pis2,mu),
                 Pow(selected.2D.hhType2,mu),
                 Pow(selected.2D.pis2.hhType2,mu),
                 Pow(selected.2D.hhType3,mu),
                 Pow(selected.2D.pis2.hhType3,mu)
    )
    
    names(fdp.res)[(length(fdp.res)-5):
                     length(fdp.res)] <- c("2D.Type1 ","2D.pis2.Type1 ",
                                                  "2D.Type2 ","2D.pis2.Type2 ",
                                                  "2D.Type3 ","2D.pis2.Type3 ")
    
  
  
  names(pow.res) <- names(fdp.res)
  print(fdp.res)
  print(pow.res)
  return(list(fdp.res = fdp.res,
              pow.res = pow.res,
              hh.extra.seq = hh.extra.seq,
              hh.current.seq = hh.current.seq,
              seed = seed))
}

##================================================
## Data adaptive neighbor selection method
##================================================
hh.Sel = function(X.sel){
  m = dim(X.sel)[2]
  Pmat=NULL
  ## We calculate the proportions 
  for(hhind in 1:6){  ## num of neighbor ranges from 2 to 7
    hh = hhind+1
    hh.seq <- rep(hh,m)
    
    pp0 <- Prop0(hh = hh.seq,
                 X = X.sel,
                 Dist = Dist.p,
                 Sigma.eps = Sigma.eps.est,
                 detect.m = detect.m)
    
    Pmat = cbind(Pmat,pp0$Prop.neigh)
  }
  Pcum = apply(Pmat==1,1,cumsum)
  ## make sure that the test statistics for all neighors are larger than 0
  compP = Pcum == matrix(1:6,nrow=dim(Pcum)[1],
                         ncol =dim(Pcum)[2],byrow=F)
  ind.init = colSums(compP)
  
  ind.init[ind.init==0] = 3
  ind.final = ind.init +1
  return(ind.final)
}

## We calculate the proportion of the neighors whose test statistics is greater than zero
Prop0 = function (hh, X, Dist, Sigma.eps, 
                  detect.m = "top.k", mu = NULL){
  m <- dim(X)[2]
  n <- dim(X)[1]
  if (detect.m == "top.k") {
    Prop.neigh <- lapply(1:m, function(i) {
      ind1 <- order(Dist[i, ])[2:(hh[i] + 1)]
      mean(colMeans(X)[ind1]>0)
    })
  }
  
  return(list(Prop.neigh=unlist(Prop.neigh)))
}



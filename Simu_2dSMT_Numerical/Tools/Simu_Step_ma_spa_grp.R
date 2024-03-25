#==================
# Simulation Step
# 1 dimension
# To evaluate the combination of spatial information and covariate information
#==================
one_step_1D_spa_grp <- function(h, grp,detect.m = "rad.h", seed = 0,
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
    X <- MASS::mvrnorm(n = n, mu = mu,
                       Sigma = Sigma.eps.p)
    X <- matrix(X, nrow = n)
    
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
    
    #==== Calculate pis2 (groupwise)
    pis.hata.grp  <- pis.hata
    for(grp.val in unique(grp)){
      ind = which(grp==grp.val)
      pis.hata.tmp<- min(1,
                         mean(mu[ind]<lambda.test)/pnorm(lambda.test))
      pis.hata.grp[ind] = pis.hata.tmp
    }
    
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
    
    #===========================================
    #==== Without spatial Info: One D (with groupwise pis2)
    #===========================================
    #==== Initilize weight function (For LAWS and SABHA respectively)
    ws.laws.fun <- function(x){(1-x)/x}
    ws.sabha.fun <- function(x){1/x}
    res.1D.pis2.Grp.SA <- OneD_Detect_w(Tm, q, pis = pis.hata.grp,
                                       const = const,
                                       ws.fun = ws.sabha.fun,
                                       tau.tm = tau.tm)
    tm.star.pis2.Grp.SA <- res.1D.pis2.Grp.SA$tm.min
    max.rej.pis2.Grp.SA <- res.1D.pis2.Grp.SA$max.rej
    fdp.res <- c(fdp.res,fdp(res.1D.pis2.Grp.SA$selected,mu))
    pow.res <- c(pow.res,Pow(res.1D.pis2.Grp.SA$selected,mu))
    NP.res <- c(NP.res,length(res.1D.pis2.Grp.SA$selected))
    
    
    res.1D.pis2.Grp.LAWS <- OneD_Detect_w(Tm, q, pis = pis.hata.grp,
                                       const = const,
                                       ws.fun = ws.laws.fun,
                                       tau.tm = tau.tm)
    tm.star.pis2.Grp.LAWS <- res.1D.pis2.Grp.LAWS$tm.min
    max.rej.pis2.Grp.LAWS <- res.1D.pis2.Grp.LAWS$max.rej
    fdp.res <- c(fdp.res,fdp(res.1D.pis2.Grp.LAWS$selected,mu))
    pow.res <- c(pow.res,Pow(res.1D.pis2.Grp.LAWS$selected,mu))
    NP.res <- c(NP.res,length(res.1D.pis2.Grp.LAWS$selected))
    
    
    fdp.res <- c(fdp(res.1D$selected,mu),
                 fdp(res.1D.pis2$selected,mu),
                 fdp(res.1D.pis2.Grp.SA$selected,mu),fdp(res.1D.pis2.Grp.LAWS$selected,mu))
    pow.res <- c(Pow(res.1D$selected,mu),
                 Pow(res.1D.pis2$selected,mu),
                 Pow(res.1D.pis2.Grp.SA$selected,mu),Pow(res.1D.pis2.Grp.LAWS$selected,mu))
    #=== Rename
    names(fdp.res) <- c("1D",
                        "1D.pis2","1D.pis2.Grp.SA","1D.pis2.Grp.LAWS")
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
  
  for(hh in h){
    print(paste("seed:",seed,"hh:",hh)) 
    
    #==== Set pre-chosen Tm, Ta
    # Tm.star
    Tm.star <- round(Tm.star,dig)
    #Tm.star.laws <- round(Tm.star.laws,dig)
    Tm.star.pis2 <- round(Tm.star.pis2,dig)
    # No need to be round because pws hasn't been round
    #tm.star.sabha <- round(tm.star.sabha,dig)
    #tm.star.ihw <- Tm.star.ihw
    #tm.star.ihw.null <- Tm.star.ihw.null
    # Ta.star
    Ta.star <- Inf
    Ta.star.pis2 <- Inf
    ta.star.pis2.Grp.SA <- 0
    ta.star.pis2.Grp.LAWS <- 0
    #===========================================
    #====== Fix Reigion
    #===========================================
    
    # Detect Neighbor
    hh.seq <- rep(hh,m)
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
    
    #== round numerical values to ensures the correctness of taking maximum
    Tm <- round(Tm,dig)
    Ta <- round(Ta,dig)
    Va <- round(Va,dig)
    VmVa.cov <- round(VmVa.cov,dig)
    
    #===========================================
    #== Run 2D Selection
    #===========================================
    print(paste("seed:",seed,"BH","hh:",hh))
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
    print(paste("seed:",seed,"BH","hh:",hh,"final.fdr",res.2D$final.fdr))
    
    print(paste("seed:",seed,"BH_pis2","hh:",hh))
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
    print(paste("seed:",seed,"BH.pis2","hh:",hh,"final.fdr",res.2D.pis2$final.fdr))
    
    print(paste("seed:",seed,"BH_pis2.Grp.SA","hh:",hh))
    res.2D.pis2.Grp.SA <- TwoDSMT::Spatial_Detect_exact_BH_down_reTm_reTa(Tm, Ta, Va,
                                                                   VmVa.cov, ind,
                                                                   q,
                                                                   max.rej = max.rej.pis2.Grp.SA,
                                                                   pis = pis.hata.grp,
                                                                   pws.tm.star = tm.star.pis2.Grp.SAs,
                                                                   pws.ta.star = ta.star.pis2.Grp.SA,
                                                                   const = const,
                                                                   ws.fun = ws.sabha.fun,
                                                                   seed = seed,
                                                                   mua = mua,
                                                                   EmpMethod = EmpMethod,
                                                                   is.quick.stop = is.quick)
    selected.2D.pis2.Grp.SA <- res.2D.pis2.Grp.SA$selected
    tm <- res.2D.pis2.Grp.SA$tm0
    ta <- res.2D.pis2.Grp.SA$ta0
    print(paste("seed:",seed,"BH_laws","hh:",hh,"final.fdr",res.2D.pis2.Grp.SA$final.fdr))
    
    print(paste("seed:",seed,"BH_pis2.Grp.LAWS","hh:",hh))
    res.2D.pis2.Grp.LAWS <- TwoDSMT::Spatial_Detect_exact_BH_down_reTm_reTa(Tm, Ta, Va,
                                                                          VmVa.cov, ind,
                                                                          q,
                                                                          max.rej = max.rej.pis2.Grp.LAWS,
                                                                          pis = pis.hata.grp,
                                                                          pws.tm.star = tm.star.pis2.Grp.LAWS,
                                                                          pws.ta.star = ta.star.pis2.Grp.LAWS,
                                                                          const = const,
                                                                          ws.fun = ws.laws.fun,
                                                                          seed = seed,
                                                                          mua = mua,
                                                                          EmpMethod = EmpMethod,
                                                                          is.quick.stop = is.quick)
    selected.2D.pis2.Grp.LAWS <- res.2D.pis2.Grp.LAWS$selected
    tm <- res.2D.pis2.Grp.LAWS$tm0
    ta <- res.2D.pis2.Grp.LAWS$ta0
    print(paste("seed:",seed,"BH_laws","hh:",hh,"final.fdr",res.2D.pis2.Grp.LAWS$final.fdr))
    
    
    fdp.res <- c(fdp.res,
                 fdp(selected.2D,mu),
                 fdp(selected.2D.pis2,mu),
                 fdp(selected.2D.pis2.Grp.SA,mu),
                 fdp(selected.2D.pis2.Grp.LAWS,mu)
    )
    
    pow.res <- c(pow.res,
                 Pow(selected.2D,mu),
                 Pow(selected.2D.pis2,mu),
                 Pow(selected.2D.pis2.Grp.SA,mu),
                 Pow(selected.2D.pis2.Grp.LAWS,mu)
    )
    
    names(fdp.res)[(length(fdp.res)-3):
                     length(fdp.res)] <- paste0(c("2D ","2D.pis2 ",
                                                  "2D.pis2.Grp.SA ","2D.pis2.Grp.LAWS "),
                                                hh)
    
  }
  
  names(pow.res) <- names(fdp.res)
  print(fdp.res)
  print(pow.res)
  return(list(fdp.res = fdp.res,
              pow.res = pow.res,
              seed = seed))
}





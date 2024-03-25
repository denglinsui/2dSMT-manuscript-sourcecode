#==================
# Simulation Step
# 1 dimension
#==================
one_step_1D <- function(h, detect.m = "rad.h", seed = 0,
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
    
    #==== Initilize weight function (For LAWS and SABHA respectively)
    ws.laws.fun <- function(x){(1-x)/x}
    ws.sabha.fun <- function(x){1/x}
    #save.image(paste0("mu_type:",mu_type," Cov_type",Cov_type,"seed",seed,".RData"))
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
    #=== Without spatial Info: BH
    p.value <- 1 - pnorm(Tm)
    result <- qvalue(p.value,pi0 = 1)
    selected.BH <- which(result$qvalues<=q)
    fdp(selected.BH,mu);Pow(selected.BH,mu)
    fdp.res <- fdp(selected.BH,mu)
    pow.res <- Pow(selected.BH,mu)
    NP.res <- length(selected.BH)
    
    #==== Calculate pis
    print("Start calculating pis(LAWS)....")
    bh.th<-bh.func(p.value, 0.9)$th
    if(bh.th == 1){
      pis.hat <- rep(1,m)
    }else{
      # Change to one-side case
      pis_1D.func <- function (x, tau = 0.1, h = 50) 
      {
        m <- length(x)
        s <- 1:m
        pval <- 1-pnorm(x)
        p.est <- rep(0, m)
        for (i in 1:m) {
          kht <- dnorm(s - i, 0, h)
          p.est[i] <- sum(kht[which(pval >= tau)])/((1 - tau) * 
                                                      sum(kht))
        }
        p.est[which(p.est > 1)] <- 1
        return(1 - p.est)
      }
      #pis.hat<-pis_1D.func.ker.reg(Tm, point, tau=bh.th)
      pis.hat <- pis_1D.func(Tm,tau=bh.th,h=50)
    }
    
    #=== With Spatial Info
    #==== Calculate qhat in SABHA
    #==== We use the parameter in Ring and Li
    print("Start calculating qhat(SABHA)....")
    tau = 0.5;eps = 0.1; TV_bd = 2
    alpha_ADMM = 10^2; beta = 10^3; eta = 5; max_iters = 5000; converge_thr = 1e-4 # parameters for ADMM
    ADMM_params = c(alpha_ADMM,beta,eta,max_iters,converge_thr)
    
    qhat = Solve_q_TV_1dim(p.value, tau, eps, TV_bd, ADMM_params)
    
    #==== Calculate pis2
    lambda.test <- 0.5#median(Tm)#quantile(Tm,0.9)#
    pis.hata.tmp<- min(1,
                       mean(Tm<lambda.test)/pnorm(lambda.test))
    pis.hata <- rep(pis.hata.tmp,m)
    #pis.hata <- (Tm<lambda.test)/pnorm(lambda.test)
    
    
    #== Round numerical values to ensures the correctness of taking maximum
    dig <- 7
    
    pis.hat <- round(pis.hat,dig)
    pis.hata <- round(pis.hata,dig)
    qhat <- round(qhat,dig)
    
    #================================
    #==== Without spatial Info: One D
    #==== Storey
    #================================
    print(paste0("seed",seed," Start Storey..........."))
    Storey.res <- qvalue(p.value,pi0 = pis.hata.tmp)
    Storey.selected <- which(Storey.res$qvalues<=q)
    fdp.res <- c(fdp.res,fdp(Storey.selected,mu))
    pow.res <- c(pow.res,Pow(Storey.selected,mu))
    NP.res <- c(NP.res,length(Storey.selected))
    
    #===========================================
    #==== LAWs
    #===========================================
    print(paste0("seed",seed," Start LAWS..........."))
    law.res <- law.func(pvs=p.value, pis.hat, q)
    law.selected <- which(law.res$de==1)
    fdp.res <- c(fdp.res,fdp(law.selected,mu))
    pow.res <- c(pow.res,Pow(law.selected,mu))
    NP.res <- c(NP.res,length(law.selected))
    
    #===========================================
    #==== SABHA
    #==== We use the parameter in Ring and Li
    #===========================================
    print(paste0("seed",seed," Start SABHA..........."))
    SABHA_method = function(pvals, qhat, alpha, tau){
      pvals[pvals>tau] = Inf
      khat=max(c(0,which(sort(qhat*pvals)<=alpha*(1:length(pvals))/length(pvals))))
      which(qhat*pvals<=alpha*khat/length(pvals))
    }
    
    sab.selected <- SABHA_method(p.value, qhat, q, tau)
    fdp.res <- c(fdp.res, fdp(sab.selected, mu))
    pow.res <- c(pow.res, Pow(sab.selected, mu))
    NP.res <- c(NP.res,length(sab.selected))
    
    #===========================================
    #==== AdaMT
    #==== Lei and Fithian
    #===========================================
    print(paste0("seed",seed," Start AdaMT..........."))
    formula <- paste0("ns(x1, df=6)")
    res_gam <- adapt_gam(x = data.frame(x1=point[,1]), pvals = p.value,
                         pi_formulas = formula, mu_formulas = formula,alphas = q)
    adaMT.selected <- res_gam$rejs
    fdp.res <- c(fdp.res, fdp(adaMT.selected, mu))
    pow.res <- c(pow.res, Pow(adaMT.selected, mu))
    NP.res <- c(NP.res,length(adaMT.selected))
    
    
    #===========================================
    #==== CAMT
    #==== Zhang and Chen
    #===========================================
    print(paste0("seed",seed," Start CAMT...........(Not Run)"))
    x1 = point[,1]
    #CAMT.res <- camt.fdr(pvals = p.value, pi0.var = ns(x1,6), f1.var = ns(x1,6),
    #                     alg.type = "EM", control.method = "knockoff+")
    #CAMT.selected <- which(CAMT.res$fdr<q)
    
    spl.num <- 6
    CAMT.res <- NULL
    while(is.null(CAMT.res)){
      CAMT.res <- tryCatch(camt.fdr(pvals = p.value, pi0.var = ns(x1,spl.num), 
                                    f1.var = ns(x1,spl.num),
                                    alg.type = "EM", control.method = "knockoff+"),
                           error = function(e){NULL})
      
      spl.num <- spl.num-1
    }
    
    CAMT.selected <- which(CAMT.res$fdr<q)
    #plot(c(1:900),CAMT.res$pi0)
    #CAMT.selected <- integer(0)
    fdp.res <- c(fdp.res, fdp(CAMT.selected, mu))
    pow.res <- c(pow.res, Pow(CAMT.selected, mu))
    NP.res <- c(NP.res,length(CAMT.selected))
    
    #===========================================
    #==== FDRreg
    #==== Scott etc.
    #===========================================
    print(paste0("seed",seed," Start FDRregT..........."))
    FDRregT.res <- FDRreg(Tm, covars=ns(x1,6),nulltype= 'theoretical')
    FDRregT.selected <- which(FDRregT.res$FDR<q)
    
    print(paste0("seed",seed," Start FDRregE..........."))
    #FDRregE.res <- FDRreg(Tm, covars=matrix(x1),nulltype= 'empirical')
    FDRregE.res <- FDRregT.res
    FDRregE.selected <- which(FDRregE.res$FDR<q)
    #pis0<-plot(c(1:900),FDRregT.res$postprob)
    fdp.res <- c(fdp.res, fdp(FDRregT.selected, mu))#, fdp(FDRregE.selected, mu))
    pow.res <- c(pow.res, Pow(FDRregT.selected, mu))#, Pow(FDRregE.selected, mu))
    NP.res <- c(NP.res,length(FDRregT.selected))#,length(FDRregE.selected))
    
    #===========================================
    #==== dBH
    #==== Fithian and Lei (One Order)
    #===========================================
    dBH.res <- dBH_mvgauss(zvals = Tm, Sigma = Sigma.eps.est,
                           side = "right", alpha = q,
                           gamma = 1, niter = 1, avals_type = "BH")
    dBH.selected <- dBH.res$rejs
    fdp.res <- c(fdp.res, fdp(dBH.selected, mu))
    pow.res <- c(pow.res, Pow(dBH.selected, mu))
    NP.res <- c(NP.res,length(dBH.selected))
    
    #===========================================
    #--- IHW: We use the location of IHW as covariate
    # (Without null proportion estimation)
    #===========================================
    ihw.res <- ihw(1-pnorm(Tm), point[,1], q,nbins = nbins.ihw)
    ihw.selected <- which(ihw.res@df$adj_pvalue<q)
    fdp.res <- c(fdp.res, fdp(ihw.selected, mu))
    pow.res <- c(pow.res, Pow(ihw.selected, mu))
    NP.res <- c(NP.res,length(ihw.selected))
    ihw.ws <- ihw.res@df$weight # Extract weights
    
    #===========================================
    #--- IHW: We use the location of IHW as covariate
    # (With null proportion estimation)
    #===========================================
    ihw.null.res <- ihw(1-pnorm(Tm), point[,1], q,nbins = nbins.ihw, null_proportion=T)
    ihw.null.selected <- which(ihw.null.res@df$adj_pvalue<q)
    fdp.res <- c(fdp.res, fdp(ihw.null.selected, mu))
    pow.res <- c(pow.res, Pow(ihw.null.selected, mu))
    NP.res <- c(NP.res,length(ihw.null.selected))
    ihw.null.ws <- ihw.null.res@df$weight # Extract weights
    
    
    pis.ihw <- ihw.ws/ihw.null.ws
    pis.ihw[is.na(pis.ihw)] <- 1 # NA, represents the null proportion is 1. Two weights are both 0.
    
    #==== Initilize for 2D detection ==========
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
    #==== Without spatial Info: One D (with pis)
    #==== LAWS based
    #===========================================
    res.1D.laws <- OneD_Detect_w(Tm, q, pis = 1-pis.hat,
                                 ws.fun = ws.laws.fun,
                                 const = const,
                                 tau.tm = tau.tm)
    tm.star.laws <- res.1D.laws$tm.min
    max.rej.laws <- res.1D.laws$max.rej
    fdp.res <- c(fdp.res,fdp(res.1D.laws$selected,mu))
    pow.res <- c(pow.res,Pow(res.1D.laws$selected,mu))
    NP.res <- c(NP.res,length(res.1D.laws$selected))
    
    
    
    #===========================================
    #==== Without spatial Info: One D (with qhat)
    #==== SABHA based
    #===========================================
    res.1D.sabha <- OneD_Detect_w(Tm, q, pis = qhat,
                                  ws.fun = ws.sabha.fun,
                                  const = const,
                                  tau.tm = tau.tm)
    tm.star.sabha <- res.1D.sabha$tm.min
    max.rej.sabha <- res.1D.sabha$max.rej
    fdp.res <- c(fdp.res,fdp(res.1D.sabha$selected,mu))
    pow.res <- c(pow.res,Pow(res.1D.sabha$selected,mu))
    NP.res <- c(NP.res,length(res.1D.sabha$selected))
    
    
    
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
    #==== Without spatial Info: One D (with IHW)
    #==== IHW based (Without null proportion)
    #===========================================
    res.1D.ihw <- OneD_Detect_w(Tm, q,
                                ws = ihw.ws,
                                const = const,
                                tau.tm = tau.tm)
    tm.star.ihw <- res.1D.ihw$tm.min
    max.rej.ihw <- res.1D.ihw$max.rej
    fdp.res <- c(fdp.res,fdp(res.1D.ihw$selected,mu))
    pow.res <- c(pow.res,Pow(res.1D.ihw$selected,mu))
    NP.sabha.res <- length(res.1D.ihw$selected)
    
    #===========================================
    #==== Without spatial Info: One D (with IHW)
    #==== IHW based (With null proportion)
    #===========================================
    res.1D.ihw.null <- OneD_Detect_w(Tm, q, pis = pis.ihw,
                                     ws = ihw.null.ws,
                                     const = const,
                                     tau.tm = tau.tm)
    tm.star.ihw.null <- res.1D.ihw.null$tm.min
    max.rej.ihw.null <- res.1D.ihw.null$max.rej
    fdp.res <- c(fdp.res,fdp(res.1D.ihw.null$selected,mu))
    pow.res <- c(pow.res,Pow(res.1D.ihw.null$selected,mu))
    NP.sabha.res <- length(res.1D.ihw.null$selected)
    
    fdp.res <- c(fdp(selected.BH,mu),fdp(Storey.selected,mu),fdp(law.selected,mu), fdp(sab.selected, mu), 
                 fdp(adaMT.selected, mu),fdp(CAMT.selected, mu), fdp(FDRregT.selected, mu),
                 fdp(dBH.selected, mu),fdp(ihw.selected, mu),fdp(ihw.null.selected, mu),
                 fdp(res.1D$selected,mu),fdp(res.1D.laws$selected,mu),fdp(res.1D.sabha$selected,mu),
                 fdp(res.1D.pis2$selected,mu),fdp(res.1D.ihw$selected,mu),fdp(res.1D.ihw.null$selected,mu))
    pow.res <- c(Pow(selected.BH,mu),Pow(Storey.selected,mu),Pow(law.selected,mu), Pow(sab.selected, mu), 
                 Pow(adaMT.selected, mu),Pow(CAMT.selected, mu), Pow(FDRregT.selected, mu),
                 Pow(dBH.selected, mu),Pow(ihw.selected, mu),Pow(ihw.null.selected, mu),
                 Pow(res.1D$selected,mu),Pow(res.1D.laws$selected,mu),Pow(res.1D.sabha$selected,mu),
                 Pow(res.1D.pis2$selected,mu),Pow(res.1D.ihw$selected,mu),Pow(res.1D.ihw.null$selected,mu))
    #=== Rename
    names(fdp.res) <- c("BH","Storey","LAWS","SABHA",
                        "AdaMT","CAMT","FDRreg(T)",
                        #"FDRreg(E)",
                        "dBH","IHW","IHW(NULL)",
                        "1D","1D.laws","1D.sabha", 
                        "1D.pis2","1D.ihw","1D.ihw.null")
    #=== With Spatial Info
    #for(ii in 1:length(h)){
    #  hh <- h[ii]
    #---- Save Result
    #save.file.path <- paste0("Result/Simulation1D/Tmp/mag_",
    #                         magnitude,"_mu_",mu_type,"_cov_",Cov_type,"_seed_",seed,".RData")
    save(list = ls(all.names = TRUE), file =save.file.path)
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
    ta.star.sabha <- 0
    ta.star.laws <- 0
    ta.star.ihw <- 0
    ta.star.ihw.null <- 0
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
    
    
    print(paste("seed:",seed,"BH.rect","hh:",hh))
    res.2D.rect <- TwoDSMT::Spatial_Detect_exact_grp_BH_down(Tm, Ta, Va,
                                                        VmVa.cov, ind,
                                                        q, max.rej,
                                                        Tm.star = Tm.star,
                                                        Ta.star = Ta.star,
                                                        const = const,
                                                        seed = seed,
                                                        mua = mua,
                                                        is.fixRectShape = T,
                                                        sig_ratio = sig_ratio,
                                                        EmpMethod = EmpMethod,
                                                        is.quick.stop = is.quick)
    selected.2D.rect <- res.2D.rect$selected
    tm <- res.2D.rect$tm0
    ta <- res.2D.rect$ta0
    print(paste("seed:",seed,"BH.rect","hh:",hh,"final.fdr",res.2D.rect$final.fdr))
    
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
    
    print(paste("seed:",seed,"BH_laws","hh:",hh))
    res.2D.laws <- TwoDSMT::Spatial_Detect_exact_BH_down_reTm_reTa(Tm, Ta, Va,
                                                                   VmVa.cov, ind,
                                                                   q,
                                                                   max.rej = max.rej.laws,
                                                                   pis = 1-pis.hat,
                                                                   pws.tm.star = tm.star.laws,
                                                                   pws.ta.star = ta.star.laws,
                                                                   const = const,
                                                                   ws.fun = ws.laws.fun,
                                                                   seed = seed,
                                                                   mua = mua,
                                                                   EmpMethod = EmpMethod,
                                                                   is.quick.stop = is.quick)
    selected.2D.laws <- res.2D.laws$selected
    tm <- res.2D.laws$tm0
    ta <- res.2D.laws$ta0
    print(paste("seed:",seed,"BH_laws","hh:",hh,"final.fdr",res.2D.laws$final.fdr))
    
    
    print(paste("seed:",seed,"BH_sabha","hh:",hh))
    res.2D.sabha <- TwoDSMT::Spatial_Detect_exact_BH_down_reTm_reTa(Tm, Ta, Va,
                                                                    VmVa.cov, ind,
                                                                    q,
                                                                    max.rej = max.rej.sabha,
                                                                    pis = qhat,
                                                                    pws.tm.star = tm.star.sabha,
                                                                    pws.ta.star = ta.star.sabha,
                                                                    const = const,
                                                                    ws.fun = ws.sabha.fun,
                                                                    n.group.max = n.group.max,
                                                                    seed = seed,
                                                                    mua = mua,
                                                                    EmpMethod = EmpMethod,
                                                                    is.quick.stop = is.quick)
    selected.2D.sabha <- res.2D.sabha$selected
    tm <- res.2D.sabha$tm0
    ta <- res.2D.sabha$ta0
    print(paste("seed:",seed,"BH_sabha","hh:",hh,"final.fdr",res.2D.sabha$final.fdr))
    
    
    print(paste("seed:",seed,"BH_ihw","hh:",hh))
    res.2D.ihw <- TwoDSMT::Spatial_Detect_exact_BH_down_reTm_reTa(Tm, Ta, Va,
                                                                  VmVa.cov, ind,
                                                                  q,
                                                                  max.rej = max.rej.ihw,
                                                                  pis = rep(1,m),
                                                                  pws.tm.star = tm.star.ihw,
                                                                  pws.ta.star = ta.star.ihw,
                                                                  const = const,
                                                                  ws = ihw.ws,
                                                                  n.group.max = n.group.max,
                                                                  seed = seed,
                                                                  mua = mua,
                                                                  EmpMethod = EmpMethod,
                                                                  is.quick.stop = is.quick)
    selected.2D.ihw <- res.2D.ihw$selected
    tm <- res.2D.ihw$tm0
    ta <- res.2D.ihw$ta0
    print(paste("seed:",seed,"BH_ihw","hh:",hh,"final.fdr",res.2D.ihw$final.fdr))
    
    
    print(paste("seed:",seed,"BH_ihw.null","hh:",hh))
    res.2D.ihw.null <- TwoDSMT::Spatial_Detect_exact_BH_down_reTm_reTa(Tm, Ta, Va,
                                                                       VmVa.cov, ind,
                                                                       q,
                                                                       max.rej = max.rej.ihw.null,
                                                                       pis = pis.ihw,
                                                                       pws.tm.star = tm.star.ihw.null,
                                                                       pws.ta.star = ta.star.ihw.null,
                                                                       const = const,
                                                                       ws = ihw.null.ws,
                                                                       n.group.max = n.group.max,
                                                                       seed = seed,
                                                                       mua = mua,
                                                                       EmpMethod = EmpMethod,
                                                                       is.quick.stop = is.quick)
    selected.2D.ihw.null <- res.2D.ihw.null$selected
    tm <- res.2D.ihw.null$tm0
    ta <- res.2D.ihw.null$ta0
    print(paste("seed:",seed,"BH_ihw.null","hh:",hh,"final.fdr",res.2D.ihw.null$final.fdr))
    
    
    fdp.res <- c(fdp.res,
                 fdp(selected.2D,mu),
                 fdp(selected.2D.laws,mu),
                 fdp(selected.2D.sabha,mu),
                 fdp(selected.2D.pis2,mu),
                 fdp(selected.2D.ihw,mu),
                 fdp(selected.2D.ihw.null,mu)
    )
    
    pow.res <- c(pow.res,
                 Pow(selected.2D,mu),
                 Pow(selected.2D.laws,mu),
                 Pow(selected.2D.sabha,mu),
                 Pow(selected.2D.pis2,mu),
                 Pow(selected.2D.ihw,mu),
                 Pow(selected.2D.ihw.null,mu)
    )
    
    names(fdp.res)[(length(fdp.res)-6):
                     length(fdp.res)] <- paste0(c("2D ","2D.rect ","2D.laws ","2D.sabha ",
                                                  "2D.pis2 ", "2D.ihw ","2D.ihw.null "),
                                                hh)
    
  }
  
  names(pow.res) <- names(fdp.res)
  print(fdp.res)
  print(pow.res)
  return(list(fdp.res = fdp.res,
              pow.res = pow.res,
              seed = seed))
}


#==================
# Simulation Step
# 2 dimension
#==================
one_step_2D <- function(h, detect.m = "rad.h", seed = 0,
                        data.type = "mvnorm",
                        const = 1,
                        ws.fun = function(x){1/x},
                        nbins.ihw = 5,
                        n.group.max = 5,
                        tau.tm = 1,
                        tau.ta = 1,
                        magnitude = 0,
                        mu_type = "Sparse",
                        Cov_type = "Weak",
                        mu = NULL,
                        EmpMethod = "NPEB",
                        ...){
  set.seed(seed)
  print(paste("seed",seed))
  
  save.folder.name <- "Simulation2D"
  
  tmpfoldername = paste0("Result/",save.folder.name,"/Tmp/")
  if(!dir.exists(tmpfoldername)){
    dir.create(tmpfoldername)
  }
  
  save.file.path <- paste0("Result/",save.folder.name,"/Tmp/mag_",
                           magnitude,"_mu_",mu_type,"_cov_",Cov_type,"_seed_",seed,".RData")
  if(!file.exists(save.file.path)){  
    pluszero <- mu>0
    
    X <- MASS::mvrnorm(n = n, mu = mu, 
                       Sigma = Sigma.eps.p)
    X <- matrix(X, nrow = n)
    
    sgm <- sqrt(diag(Sigma.eps.p))
    
    #==== Initilize weight function (For LAWS and SABHA respectively)
    ws.laws.fun <- function(x){(1-x)/x}
    ws.sabha.fun <- function(x){1/x}
    #==== Perform Algorithm
    dig <- 7
    Tm.acc <- apply(X,2,function(x){sum(x)/sqrt(n)})/sgm
    Tm.acc <- round(Tm.acc,dig)
    #==== Estimate covariance with two approach.
    #==== We only cares about the correlation and assume the variance is known to be one.
    if(estcov == T){
      corrmodel <- "exponential"
      
      #Type 1: Multiple Observation
      print(paste0("seed",seed,":Start Covariance Estimation (Full & Multi-Obs)...."))
      
      X.demean <- apply(X,2,function(x){x-mean(x)})
      est.ind <- 1:m
      point.for.cov <- point[est.ind,]
      rem <- c(0.1,0.5,0.5)
      names(rem) <- c("scale","sill","nugget")
      fit.m_obs <- FitComposite(X.demean[1:(n-1),est.ind], coordx=point.for.cov,
                                corrmodel=corrmodel, likelihood="Full",
                                fixed = list(mean=0), 
                                start = as.list(rem),
                                type="Standard", replicate = n-1)
      
      param.m_obs <- as.list(fit.m_obs$param)
      
      cov.est.m_obs <- Covmatrix(point[,1], point[,2], 
                                 corrmodel=corrmodel,
                                 param=param.m_obs)
      Sigma.eps.est.m_obs <- cov.est.m_obs$covmatrix * n/(n-1)
      
      print(paste0("seed",seed,":End Full & Multi-Obs."))
      
      #Type 2: Single Observation
      print(paste0("seed",seed,":Start Covariance Estimation (Full & Single-Obs)...."))
      # Generate covariate
      covar1 <- ns(point[,1],4)
      covar2 <- ns(point[,2],4)
      covariate <- matrix(NA,nrow=nrow(point),ncol = 16)
      
      #cbind(ns(point[,1],4),ns(point[,2],4))
      for(i in 1:4){
        for(j in 1:4){
          covariate[,(i-1)*4+j] <- covar1[,i]*covar2[,j]
        }
      }
      
      rem <- c(0.1,0.5,0.5)
      names(rem) <- c("scale","sill","nugget")
      fit.s_obs <- FitComposite2(matrix(Tm.acc[est.ind],nrow=1), coordx=point.for.cov, covariate = covariate,
                                 corrmodel=corrmodel, likelihood="Full",#fixed = list(mean=0), 
                                 start = as.list(rem),
                                 type="Standard", replicate = 1,optimizer = "L-BFGS-B")
      
      param.s_obs <- as.list(fit.s_obs$param)
      
      cov.est.s_obs <- Covmatrix(point[,1], point[,2], 
                                 corrmodel=corrmodel,
                                 param=param.s_obs)
      Sigma.eps.est.s_obs <- cov.est.s_obs$covmatrix
      
      print(paste0("seed",seed,":End Full & Single-Obs."))
    }else{
      Sigma.eps.est.m_obs <- Sigma.eps.p
      Sigma.eps.est.s_obs <- Sigma.eps.p
    }
    #=== Generate Tm for testing
    sgm.est <- sqrt(diag(Sigma.eps.est.m_obs))
    Tm <- apply(X,2,function(x){sum(x)/sqrt(n)})/sgm.est
    Tm <- round(Tm,dig)
    
    
    #=== Without spatial Info: BH
    p.value <- 1 - pnorm(Tm)
    result <- qvalue(p.value,pi0 = 1)
    selected.BH <- which(result$qvalues<=q)
    
    #==== Calculate pis
    bh.th<-bh.func(p.value, 0.9)$th
    if(bh.th == 1){
      pis.hat <- rep(1,m)
    }else{
      pis_2D.func<- function(x, tau=0.1, h=10)
      {
        ## pis_2D.func calculates the conditional proportions pis
        ## Arguments
        # x: a matrix of z-values
        # tau: the screening threshold, which can be prespecified or chosen adaptively
        # bdw: bandwidth
        ## Values
        # pis: conditional proportions
        
        dims<-dim(x)
        m<-dims[1]*dims[2]
        x.vec<-c(t(x))
        # pv.vec<-pnorm(-abs(x.vec)) #Change to one-side test
        pv.vec<-1-pnorm(x.vec)
        scr.idx<-which(pv.vec>=tau)
        p.est<-matrix(rep(0, m), dims[1], dims[2])  
        
        for (i in 1:dims[1]) 
        {
          for (j in 1:dims[2]) 
          {
            s<-c(i, j)
            dis.vec<-disvec.func(dims, s)
            kht<-dnorm(dis.vec, 0, h)
            p.est[i,j]<-min(1-1e-5, sum(kht[scr.idx])/((1-tau)*sum(kht)))
          }
        }
        pis.est<-1-c(p.est)
        return(pis.est)
      }
      # We require m.x = m.y and Tm is ordered by row or column
      pis.hat<-pis_2D.func(matrix(Tm,sqrt(m),sqrt(m)), tau=bh.th, h=4.5)
    }
    
    #=== With Spatial Info
    #==== Calculate qhat in SABHA
    #==== We use the parameter in Ring and Li
    print("Start calculating qhat(SABHA)....")
    tau = 0.5;eps = 0.1; TV_bd = 2
    alpha_ADMM = 10^2; beta = 10^3; eta = 5; max_iters = 5000; converge_thr = 1e-4 # parameters for ADMM
    ADMM_params = c(alpha_ADMM,beta,eta,max_iters,converge_thr)
    
    # We require m.x = m.y and Tm is ordered by row or column
    Pvals_reshape <- matrix(p.value,sqrt(m),sqrt(m))
    #qhat = Solve_q_TV_1dim(p.value, tau, eps, TV_bd, ADMM_params)
    qhat = Solve_q_TV_2dim(Pvals_reshape, tau, eps, TV_bd, ADMM_params)
    #==== Calculate pis2
    lambda.test <- 0.5#quantile(Tm,0.9)#
    pis.hata.tmp<- min(1,
                       mean(Tm<lambda.test)/pnorm(lambda.test))
    pis.hata <- rep(pis.hata.tmp,m)
    #pis.hata <- (Tm<lambda.test)/pnorm(lambda.test)
    #== Round numerical values to ensures the correctness of taking maximum
    dig <- 7
    
    pis.hat <- round(as.vector(pis.hat),dig)
    pis.hata <- round(pis.hata,dig)
    qhat <- round(as.vector(qhat),dig)
    
    #================================
    #==== Without spatial Info: One D
    #==== Storey
    #================================
    print(paste0("seed",seed," Start Storey..........."))
    Storey.res <- qvalue(p.value,pi0 = pis.hata.tmp)
    Storey.selected <- which(Storey.res$qvalues<=q)
    
    #===========================================
    #==== LAWs
    #===========================================
    print(paste0("seed",seed," Start LAWS..........."))
    law.res <- law.func(pvs=p.value, pis.hat, q)
    law.selected <- which(law.res$de==1)
    
    #===========================================
    #==== SABHA
    #==== We use the parameter in Ring and Li
    #===========================================
    print(paste0("seed",seed," Start SABHA..........."))
    SABHA_method = function(pvals, qhat, alpha, tau){
      pvals[pvals>tau] = Inf
      khat=max(c(0,which(sort(qhat*pvals)<=alpha*(1:length(pvals))/length(pvals))))
      which(qhat*pvals<=alpha*khat/length(pvals))
    }
    
    sab.selected <- SABHA_method(p.value, qhat, q, tau)
    
    #===========================================
    #==== AdaMT
    #==== Lei and Fithian
    #===========================================
    print(paste0("seed",seed," Start AdaMT..........."))
    formula <- c("ns(x1, df=6)","ns(x2, df=6)") # Time Consuming
    res_gam <- adapt_gam(x = data.frame(x1=point[,1],x2=point[,2]), pvals = p.value,
                         pi_formulas = formula, mu_formulas = formula,alphas = q)
    adaMT.selected <- res_gam$rejs
    
    #===========================================
    #==== CAMT
    #==== Zhang and Chen
    #===========================================
    print(paste0("seed",seed," Start CAMT...........(Not Run)"))
    covariate = cbind(ns(point[,1],6),ns(point[,2],6))
    CAMT.res <- camt.fdr(pvals = p.value, pi0.var = covariate,
                         f1.var = covariate,
                         alg.type = "EM", control.method = "knockoff+")
    CAMT.selected <- which(CAMT.res$fdr<q)
    
    #===========================================
    #==== FDRreg
    #==== Scott etc.
    #===========================================
    print(paste0("seed",seed," Start FDRregT..........."))
    FDRregT.res <- FDRreg(Tm, covars=covariate,nulltype= 'theoretical')
    FDRregT.selected <- which(FDRregT.res$FDR<q)
    
    print(paste0("seed",seed," Start FDRregE..........."))
    #FDRregE.res <- FDRreg(Tm, covars=matrix(x1),nulltype= 'empirical')
    FDRregE.res <- FDRregT.res
    FDRregE.selected <- which(FDRregE.res$FDR<q)
    
    #===========================================
    #==== dBH
    #==== Fithian and Lei (One Order)
    #==== Waiting to update (There are two types Sigma.est)
    #===========================================
    dBH.true.res <- dBH_mvgauss(zvals = Tm, Sigma = Sigma.eps.p,
                                side = "right", alpha = q,
                                gamma = 1, niter = 1, avals_type = "BH")
    dBH.true.selected <- dBH.true.res$rejs
    
    dBH.m_obs.res <- dBH_mvgauss(zvals = Tm, Sigma = Sigma.eps.est.m_obs,
                                 side = "right", alpha = q,
                                 gamma = 1, niter = 1, avals_type = "BH")
    dBH.m_obs.selected <- dBH.m_obs.res$rejs
    
    
    dBH.s_obs.res <- dBH_mvgauss(zvals = Tm, Sigma = Sigma.eps.est.s_obs,
                                 side = "right", alpha = q,
                                 gamma = 1, niter = 1, avals_type = "BH")
    dBH.s_obs.selected <- dBH.s_obs.res$rejs
    
    #===========================================
    #--- IHW: We use the location of IHW as covariate
    # (Without null proportion estimation)
    #===========================================
    ihw.res <- ihw(1-pnorm(Tm), point[,1], q,nbins = nbins.ihw)
    ihw.selected <- which(ihw.res@df$adj_pvalue<q)
    ihw.ws <- ihw.res@df$weight # Extract weights
    
    #===========================================
    #--- IHW: We use the location of IHW as covariate
    # (With null proportion estimation)
    #===========================================
    ihw.null.res <- ihw(1-pnorm(Tm), point[,1], q,nbins = nbins.ihw, null_proportion=T)
    ihw.null.selected <- which(ihw.null.res@df$adj_pvalue<q)
    ihw.null.ws <- ihw.null.res@df$weight # Extract weights
    
    pis.ihw <- ihw.ws/ihw.null.ws
    pis.ihw[is.na(pis.ihw)] <- 1 # NA, represents the null proportion is 1. Two weights are both 0.
    
    #==== Initilize for 2D detection ==========
    #================================
    #==== Without spatial Info: One D
    #================================
    res.1D <- OneD_Detect(Tm, q,  pis = NULL,
                          const = const,
                          tau.tm = tau.tm)
    Tm.star <- res.1D$tm.min
    max.rej <- res.1D$max.rej
    
    #===========================================
    #==== Without spatial Info: One D (with pis)
    #==== LAWS based
    #===========================================
    res.1D.laws <- OneD_Detect_w(Tm, q, pis = 1-pis.hat,
                                 ws.fun = ws.laws.fun,
                                 const = const,
                                 tau.tm = tau.tm)
    tm.star.laws <- res.1D.laws$tm.min
    max.rej.laws <- res.1D.laws$max.rej
    
    #===========================================
    #==== Without spatial Info: One D (with qhat)
    #==== SABHA based
    #===========================================
    res.1D.sabha <- OneD_Detect_w(Tm, q, pis = qhat,
                                  ws.fun = ws.sabha.fun,
                                  const = const,
                                  tau.tm = tau.tm)
    tm.star.sabha <- res.1D.sabha$tm.min
    max.rej.sabha <- res.1D.sabha$max.rej
    
    #===========================================
    #==== Without spatial Info: One D (with pis2)
    #===========================================
    res.1D.pis2 <- OneD_Detect(Tm, q, pis = pis.hata,
                               const = const,
                               tau.tm = tau.tm)
    Tm.star.pis2 <- res.1D.pis2$tm.min
    max.rej.pis2 <- res.1D.pis2$max.rej
    
    #===========================================
    #==== Without spatial Info: One D (with IHW)
    #==== IHW based (Without null proportion)
    #===========================================
    res.1D.ihw <- OneD_Detect_w(Tm, q,
                                ws = ihw.ws,
                                const = const,
                                tau.tm = tau.tm)
    tm.star.ihw <- res.1D.ihw$tm.min
    max.rej.ihw <- res.1D.ihw$max.rej
    
    #===========================================
    #==== Without spatial Info: One D (with IHW)
    #==== IHW based (With null proportion)
    #===========================================
    res.1D.ihw.null <- OneD_Detect_w(Tm, q, pis = pis.ihw,
                                     ws = ihw.null.ws,
                                     const = const,
                                     tau.tm = tau.tm)
    tm.star.ihw.null <- res.1D.ihw.null$tm.min
    max.rej.ihw.null <- res.1D.ihw.null$max.rej
    
    fdp.res <- c(fdp(selected.BH,mu),fdp(Storey.selected,mu),fdp(law.selected,mu), fdp(sab.selected, mu), 
                 fdp(adaMT.selected, mu),fdp(CAMT.selected, mu), fdp(FDRregT.selected, mu),
                 fdp(dBH.true.selected, mu),fdp(dBH.m_obs.selected, mu),fdp(dBH.s_obs.selected, mu),fdp(ihw.selected, mu),fdp(ihw.null.selected, mu),
                 fdp(res.1D$selected,mu),fdp(res.1D.laws$selected,mu),fdp(res.1D.sabha$selected,mu),
                 fdp(res.1D.pis2$selected,mu),fdp(res.1D.ihw$selected,mu),fdp(res.1D.ihw.null$selected,mu))
    pow.res <- c(Pow(selected.BH,mu),Pow(Storey.selected,mu),Pow(law.selected,mu), Pow(sab.selected, mu), 
                 Pow(adaMT.selected, mu),Pow(CAMT.selected, mu), Pow(FDRregT.selected, mu),
                 Pow(dBH.true.selected, mu),Pow(dBH.m_obs.selected, mu),Pow(dBH.s_obs.selected, mu),Pow(ihw.selected, mu),Pow(ihw.null.selected, mu),
                 Pow(res.1D$selected,mu),Pow(res.1D.laws$selected,mu),Pow(res.1D.sabha$selected,mu),
                 Pow(res.1D.pis2$selected,mu),Pow(res.1D.ihw$selected,mu),Pow(res.1D.ihw.null$selected,mu))
    #=== Rename
    names(fdp.res) <- c("BH","Storey","LAWS","SABHA",
                        "AdaMT","CAMT","FDRreg(T)",
                        #"FDRreg(E)",
                        "dBH(True)","dBH(m_obs)","dBH(s_obs)","IHW","IHW(NULL)",
                        "1D","1D.laws","1D.sabha", 
                        "1D.pis2","1D.ihw","1D.ihw.null")
    save(list = ls(all.names = TRUE), file =save.file.path)
  }else{
    # Load the previous result and save time.
    load(save.file.path)
  }
  
  #=== With Spatial Info
  for(Sigma.eps.type in c("Sigma.eps.p","Sigma.eps.est.m_obs")){
    #for(Sigma.eps.type in c("Sigma.eps.p","Sigma.eps.est.m_obs","Sigma.eps.est.s_obs")){
    Sigma.eps.est <- get(Sigma.eps.type)
    for(ii in 1:length(h)){
      hh <- h[ii]
      print(paste("seed:",seed,"hh:",hh))
      
      #==== Set pre-chosen Tm, Ta
      # Tm.star
      Tm.star <- round(Tm.star,dig)
      # Tm.star.pis <- round(Tm.star.pis,dig)
      Tm.star.pis2 <- round(Tm.star.pis2,dig)
      
      # Ta.star
      Ta.star <- Inf
      Ta.star.pis2 <- Inf
      ta.star.sabha <- 0
      ta.star.laws <- 0
      ta.star.ihw <- 0
      ta.star.ihw.null <- 0
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
      is.quick <- T
      print(paste("seed:",seed,"BH","hh:",hh))
      res.2D <- TwoDSMT::Spatial_Detect_exact_grp_BH_down(Tm, Ta, Va, 
                                                          VmVa.cov, ind,
                                                          q, max.rej,
                                                          Tm.star = Tm.star,
                                                          Ta.star = Ta.star,
                                                          const = const,
                                                          seed = seed,
                                                          mua = mua,
                                                          EmpMethod = EmpMethod,
                                                          is.quick.stop = is.quick)
      selected.2D <- res.2D$selected
      tm <- res.2D$tm0
      ta <- res.2D$ta0
      print(paste("seed:",seed,"BH","hh:",hh,"final.fdr",res.2D$final.fdr))
      
      print(paste("seed:",seed,"BHFixRect","hh:",hh))
      res.2D.rect <- TwoDSMT::Spatial_Detect_exact_grp_BH_down(Tm, Ta, Va, 
                                                               VmVa.cov, ind,
                                                               q, max.rej,
                                                               Tm.star = Tm.star,
                                                               Ta.star = Ta.star,
                                                               const = const,
                                                               seed = seed,
                                                               mua = mua,
                                                               is.fixRectShape = T,
                                                               EmpMethod = EmpMethod,
                                                               is.quick.stop = is.quick)
      selected.2D.rect <- res.2D.rect$selected
      tm <- res.2D.rect$tm0
      ta <- res.2D.rect$ta0
      print(paste("seed:",seed,"BHFixRect","hh:",hh,"final.fdr",res.2D.rect$final.fdr))
      
      
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
      
      print(paste("seed:",seed,"BH_laws","hh:",hh))
      res.2D.laws <- TwoDSMT::Spatial_Detect_exact_BH_down_reTm_reTa(Tm, Ta, Va,
                                                                     VmVa.cov, ind,
                                                                     q,
                                                                     max.rej = max.rej.laws,
                                                                     pis = 1-pis.hat,
                                                                     pws.tm.star = tm.star.laws,
                                                                     pws.ta.star = ta.star.laws,
                                                                     const = const,
                                                                     ws.fun = ws.laws.fun,
                                                                     seed = seed,
                                                                     mua = mua,
                                                                     EmpMethod = EmpMethod,
                                                                     is.quick.stop = is.quick)
      selected.2D.laws <- res.2D.laws$selected
      tm <- res.2D.laws$tm0
      ta <- res.2D.laws$ta0
      print(paste("seed:",seed,"BH_laws","hh:",hh,"final.fdr",res.2D.laws$final.fdr))
      
      
      print(paste("seed:",seed,"BH_sabha","hh:",hh))
      res.2D.sabha <- TwoDSMT::Spatial_Detect_exact_BH_down_reTm_reTa(Tm, Ta, Va,
                                                                      VmVa.cov, ind,
                                                                      q,
                                                                      max.rej = max.rej.sabha,
                                                                      pis = qhat,
                                                                      pws.tm.star = tm.star.sabha,
                                                                      pws.ta.star = ta.star.sabha,
                                                                      const = const,
                                                                      ws.fun = ws.sabha.fun,
                                                                      n.group.max = n.group.max,
                                                                      seed = seed,
                                                                      mua = mua,
                                                                      EmpMethod = EmpMethod,
                                                                      is.quick.stop = is.quick)
      selected.2D.sabha <- res.2D.sabha$selected
      tm <- res.2D.sabha$tm0
      ta <- res.2D.sabha$ta0
      print(paste("seed:",seed,"BH_sabha","hh:",hh,"final.fdr",res.2D.sabha$final.fdr))
      
      
      print(paste("seed:",seed,"BH_ihw","hh:",hh))
      res.2D.ihw <- TwoDSMT::Spatial_Detect_exact_BH_down_reTm_reTa(Tm, Ta, Va,
                                                                    VmVa.cov, ind,
                                                                    q,
                                                                    max.rej = max.rej.ihw,
                                                                    pis = rep(1,m),
                                                                    pws.tm.star = tm.star.ihw,
                                                                    pws.ta.star = ta.star.ihw,
                                                                    const = const,
                                                                    ws = ihw.ws,
                                                                    n.group.max = n.group.max,
                                                                    seed = seed,
                                                                    mua = mua,
                                                                    EmpMethod = EmpMethod,
                                                                    is.quick.stop = is.quick)
      selected.2D.ihw <- res.2D.ihw$selected
      tm <- res.2D.ihw$tm0
      ta <- res.2D.ihw$ta0
      print(paste("seed:",seed,"BH_ihw","hh:",hh,"final.fdr",res.2D.ihw$final.fdr))
      
      
      print(paste("seed:",seed,"BH_ihw.null","hh:",hh))
      res.2D.ihw.null <- TwoDSMT::Spatial_Detect_exact_BH_down_reTm_reTa(Tm, Ta, Va,
                                                                         VmVa.cov, ind,
                                                                         q,
                                                                         max.rej = max.rej.ihw.null,
                                                                         pis = pis.ihw,
                                                                         pws.tm.star = tm.star.ihw.null,
                                                                         pws.ta.star = ta.star.ihw.null,
                                                                         const = const,
                                                                         ws = ihw.null.ws,
                                                                         n.group.max = n.group.max,
                                                                         seed = seed,
                                                                         mua = mua,
                                                                         EmpMethod = EmpMethod,
                                                                         is.quick.stop = is.quick)
      selected.2D.ihw.null <- res.2D.ihw.null$selected
      tm <- res.2D.ihw.null$tm0
      ta <- res.2D.ihw.null$ta0
      print(paste("seed:",seed,"BH_ihw.null","hh:",hh,"final.fdr",res.2D.ihw.null$final.fdr))
      
      
      fdp.res <- c(fdp.res,
                   fdp(selected.2D,mu),
                   fdp(selected.2D.rect,mu),
                   fdp(selected.2D.laws,mu),
                   fdp(selected.2D.sabha,mu),
                   fdp(selected.2D.pis2,mu),
                   fdp(selected.2D.ihw,mu),
                   fdp(selected.2D.ihw.null,mu)
      )
      
      pow.res <- c(pow.res,
                   Pow(selected.2D,mu),
                   Pow(selected.2D.rect,mu),
                   Pow(selected.2D.laws,mu),
                   Pow(selected.2D.sabha,mu),
                   Pow(selected.2D.pis2,mu),
                   Pow(selected.2D.ihw,mu),
                   Pow(selected.2D.ihw.null,mu)
      )
      
      names(fdp.res)[(length(fdp.res)-6):
                       length(fdp.res)] <- paste0(c("2D ","2D.rect ","2D.laws ","2D.sabha ",
                                                    "2D.pis2 ", "2D.ihw ","2D.ihw.null "),
                                                  hh)
      
    }
  }
  names(pow.res) <- names(fdp.res)
  print(fdp.res)
  print(pow.res)
  return(list(fdp.res = fdp.res,
              pow.res = pow.res,
              seed = seed))
}
# 
# 
# #==================
# # Simulation Step
# # 1 dimension
# # Varying h
# #==================
# one_step_1D_h <- function(h, detect.m = "rad.h", seed = 0,
#                           data.type = "mvnorm",
#                           mu = mu,
#                           const = 1,
#                           nbins.ihw = 5,
#                           n.group.max = 5,
#                           tau.tm = 1,tau.ta =1,
#                           #ws.fun = function(x){1/x}, # It should depends on the method we use
#                           ...){
#   set.seed(seed)
#   print(paste("seed",seed))
#   # Current DataType
#   # "mvnorm", "Circle","SquareSpline","CircleMix"
#   
#   pluszero <- mu>0
#   #X <- user_mvrnorm(n = n, mu = mu, 
#   #                  Sigma = Sigma.eps.p)
#   X <- MASS::mvrnorm(n = n, mu = mu, 
#                      Sigma = Sigma.eps.p)
#   X <- matrix(X, nrow = n)
#   
#   #==== Initilize weight function (For LAWS and SABHA respectively)
#   ws.laws.fun <- function(x){(1-x)/x}
#   ws.sabha.fun <- function(x){1/x}
#   
#   #======================================
#   #==== Estimate covariance  
#   #======================================
#   if(estcov==T){
#     X.demean <- apply(X,2,function(x){x-mean(x)})
#     Sigma.eps.est <- NA
#     #est.ind <- sample(1:m,m/3)
#     est.ind <- 1:m
#     
#     # Fit Composite only works for 2 dimensional spatial domain
#     point.for.cov <- cbind(point[est.ind,],rep(1,length(est.ind)))
#     corrmodel <- "exponential"
#     
#     #fit <- FitComposite(X.demean[1:(n-1),est.ind], coordx=point.for.cov,
#     #                    corrmodel=corrmodel, likelihood="Full",fixed = list(mean=0), 
#     #                    type="Standard", replicate = n-1)
#     
#     # Tapering is more accurate, so we choose Tapering
#     fit <- FitComposite(X.demean[1:(n-1),est.ind], 
#                         coordx=point.for.cov,
#                         corrmodel=corrmodel,likelihood="Full",
#                         type="Tapering",taper="Wendland1",maxdist=10,
#                         start = list(scale=0.5,sill=0.5,nugget = 1),
#                         fixed = list(mean=0), 
#                         replicate = n-1)
#     param <- as.list(fit$param)
#     
#     # Generate the covariance matrix
#     cov.est <- Covmatrix(cbind(point,rep(1,length(point))), 
#                          corrmodel=corrmodel,
#                          param=param)
#     Sigma.eps.est <- cov.est$covmatrix * n/(n-1)
#     
#     param <- as.list(fit$param)
#     
#     cov.est <- Covmatrix(cbind(point,rep(1,length(point))), 
#                          corrmodel=corrmodel,
#                          param=param)
#     Sigma.eps.est <- cov.est$covmatrix * n/(n-1)
#   }else{
#     Sigma.eps.est <- Sigma.eps.p
#   }
#   
#   #==== End estimating covariance  
#   #=====================================
#   
#   sgm <- sqrt(diag(Sigma.eps.est))
#   #==== Perform Algorithm
#   dig <- 7
#   Tm <- apply(X,2,function(x){sum(x)/sqrt(n)})/sgm
#   Tm <- round(Tm,dig)
#   #=== Without spatial Info: BH
#   #=== (h=0)
#   p.value <- 1 - pnorm(Tm)
#   BH.res <- qvalue(p.value,pi0 = 1)
#   BH.selected <- which(BH.res$qvalues<=q)
#   #fdp(selected.BH,mu);Pow(selected.BH,mu)
#   
#   #================================
#   #==== Without spatial Info: One D
#   #================================
#   res.1D <- OneD_Detect(Tm, q, pis = NULL,
#                         const = const,
#                         tau.tm = tau.tm)
#   Tm.star <- res.1D$tm.min 
#   max.rej <- res.1D$max.rej
#   fdp.nopis.res <- fdp(res.1D$selected,mu)
#   pow.nopis.res <- Pow(res.1D$selected,mu)
#   NP.nopis.res <- length(res.1D$selected)
#   
#   #==== Calculate pis (LAWS)
#   print(paste0("seed:", seed, "========== Calculating pis =============="))
#   bh.th<-bh.func(p.value, 0.9)$th
#   if(bh.th == 1){
#     pis.hat <- rep(1,m)
#   }else{
#     pis.hat<-pis_1D.func.ker.reg(Tm, point, tau=bh.th)
#   }
#   
#   #==== Calculate qhat (SABHA)
#   print(paste0("seed:", seed, "========== Calculating qhat =============="))
#   tau = 0.5;eps = 0.1; TV_bd = 2
#   alpha_ADMM = 10^2; beta = 10^3; eta = 5; max_iters = 5000; converge_thr = 1e-4 # parameters for ADMM
#   ADMM_params = c(alpha_ADMM,beta,eta,max_iters,converge_thr)
#   
#   qhat = Solve_q_TV_1dim(p.value, tau, eps, TV_bd, ADMM_params)
#   
#   #==== Calculate pis2 (Storey)
#   lambda.test <- 0.5#median(Tm)#quantile(Tm,0.9)#
#   pis.hata.tmp<- min(1,
#                      mean(Tm<lambda.test)/pnorm(lambda.test))
#   pis.hata <- rep(pis.hata.tmp,m)
#   #pis.hata <- (Tm<lambda.test)/pnorm(lambda.test)
#   
#   #==============================================
#   #==== Run Other Algorithms: Storey,LAWS, SABHA, IHW
#   #==============================================
#   #--- Storey
#   print(paste0("seed",seed," Start Storey..........."))
#   Storey.res <- qvalue(p.value,pi0 = pis.hata.tmp)
#   Storey.selected <- which(Storey.res$qvalues<=q)
#   
#   #--- LAWS
#   print(paste0("seed",seed," Start LAWS..........."))
#   laws.res <- law.func(pvs=p.value, pis.hat, q)
#   laws.selected <- which(laws.res$de==1)
#   
#   #--- SABHA
#   print(paste0("seed",seed," Start SABHA..........."))
#   SABHA_method = function(pvals, qhat, alpha, tau){
#     pvals[pvals>tau] = Inf
#     khat=max(c(0,which(sort(qhat*pvals)<=alpha*(1:length(pvals))/length(pvals))))
#     which(qhat*pvals<=alpha*khat/length(pvals))
#   }
#   
#   sab.selected <- SABHA_method(p.value, qhat, q, tau)
#   
#   #--- IHW: We use the location of IHW as covariate
#   # (Without null proportion estimation)
#   print(paste0("seed",seed," Start IHW..........."))
#   ihw.res <- ihw(p.value, point[,1], q,nbins = nbins.ihw) 
#   ihw.selected <- which(ihw.res@df$adj_pvalue<q)
#   ihw.ws <- ihw.res@df$weight # Extract weights
#   
#   #--- IHW: We use the location of IHW as covariate
#   # (With null proportion estimation)
#   print(paste0("seed",seed," Start IHW..........."))
#   ihw.null.res <- ihw(p.value, point[,1], q,nbins = nbins.ihw, null_proportion=T) 
#   ihw.null.selected <- which(ihw.null.res@df$adj_pvalue<q)
#   ihw.null.ws <- ihw.null.res@df$weight # Extract weights
#   
#   pis.ihw <- ihw.ws/ihw.null.ws
#   pis.ihw[is.na(pis.ihw)] <- 1 # NA, represents the null proportion is 1. Two weights are both 0.
#   #--- Calculate fdp and power
#   fdp.res <- c(fdp(BH.selected,mu), fdp(Storey.selected, mu),
#                fdp(laws.selected,mu), fdp(sab.selected, mu), 
#                fdp(ihw.selected, mu),fdp(ihw.null.selected, mu))
#   pow.res <- c(Pow(BH.selected,mu), Pow(Storey.selected, mu),
#                Pow(laws.selected,mu), Pow(sab.selected, mu), 
#                Pow(ihw.selected, mu),Pow(ihw.null.selected, mu))
#   NP.res <- c(length(BH.selected), length(Storey.selected),
#               length(laws.selected),length(sab.selected),
#               length(ihw.selected),length(ihw.null.selected))
#   names(fdp.res) <- c("BH","Storey","LAWS","SABHA","IHW","IHW.null")
#   #== Round numerical values to ensures the correctness of taking maximum
#   dig <- 7
#   
#   pis.hat <- round(pis.hat,dig)
#   pis.hata <- round(pis.hata,dig)
#   qhat <- round(qhat,dig)
#   
#   #==== Initilize for 2D detection ==========
#   #===========================================
#   #==== Without spatial Info: One D (with pis)
#   #==== LAWS based
#   #===========================================
#   res.1D.laws <- OneD_Detect_w(Tm, q, pis = 1-pis.hat,
#                                ws.fun = ws.laws.fun,
#                                const = const,
#                                tau.tm = tau.tm)
#   tm.star.laws <- res.1D.laws$tm.min
#   max.rej.laws <- res.1D.laws$max.rej
#   fdp.laws.res <- fdp(res.1D.laws$selected,mu)
#   pow.laws.res <- Pow(res.1D.laws$selected,mu)
#   NP.laws.res <- length(res.1D.laws$selected)
#   
#   #===========================================
#   #==== Without spatial Info: One D (with qhat)
#   #==== SABHA based
#   #===========================================
#   res.1D.sabha <- OneD_Detect_w(Tm, q, pis = qhat,
#                                 const = const,
#                                 ws.fun = ws.sabha.fun,
#                                 tau.tm = tau.tm)
#   tm.star.sabha <- res.1D.sabha$tm.min
#   max.rej.sabha <- res.1D.sabha$max.rej
#   fdp.sabha.res <- fdp(res.1D.sabha$selected,mu)
#   pow.sabha.res <- Pow(res.1D.sabha$selected,mu)
#   NP.sabha.res <- length(res.1D.sabha$selected)
#   
#   
#   #===========================================
#   #==== Without spatial Info: One D (with IHW)
#   #==== IHW based (Without null proportion)
#   #===========================================
#   res.1D.ihw <- OneD_Detect_w(Tm, q,
#                               ws = ihw.ws,
#                               const = const,
#                               tau.tm = tau.tm)
#   tm.star.ihw <- res.1D.ihw$tm.min
#   max.rej.ihw <- res.1D.ihw$max.rej
#   fdp.ihw.res <- fdp(res.1D.ihw$selected,mu)
#   pow.ihw.res <- Pow(res.1D.ihw$selected,mu)
#   NP.ihw.res <- length(res.1D.ihw$selected)
#   
#   #===========================================
#   #==== Without spatial Info: One D (with IHW)
#   #==== IHW based (With null proportion)
#   #===========================================
#   res.1D.ihw.null <- OneD_Detect_w(Tm, q, pis = pis.ihw,
#                                    ws = ihw.null.ws,
#                                    const = const,
#                                    tau.tm = tau.tm)
#   tm.star.ihw.null <- res.1D.ihw.null$tm.min
#   max.rej.ihw.null <- res.1D.ihw.null$max.rej
#   fdp.ihw.null.res <- fdp(res.1D.ihw.null$selected,mu)
#   pow.ihw.null.res <- Pow(res.1D.ihw.null$selected,mu)
#   NP.ihw.null.res <- length(res.1D.ihw.null$selected)
#   
#   #===========================================
#   #==== Without spatial Info: One D (with pis2)
#   #---- The results from OneD_Detect and OneD_Detect_w differ
#   #===========================================
#   res.1D.pis2 <- OneD_Detect(Tm, q, pis = pis.hata,
#                              const = const,
#                              tau.tm = tau.tm)
#   Tm.star.pis2 <- res.1D.pis2$tm.min
#   max.rej.pis2 <- res.1D.pis2$max.rej
#   fdp.pis2.res <- fdp(res.1D.pis2$selected,mu)
#   pow.pis2.res <- Pow(res.1D.pis2$selected,mu)
#   NP.pis2.res <- length(res.1D.pis2$selected)
#   
#   #=== With Spatial Info
#   for(ii in 1:length(h)){
#     hh <- h[ii]
#     print(paste("seed:",seed,"hh:",hh))
#     
#     #==== Set pre-chosen Tm, Ta
#     # Tm.star
#     Tm.star <- round(Tm.star,dig)
#     #Tm.star.laws <- round(Tm.star.laws,dig)
#     Tm.star.pis2 <- round(Tm.star.pis2,dig)
#     # No need to be round because pws hasn't been round
#     #tm.star.sabha <- round(tm.star.sabha,dig)
#     
#     # Ta.star
#     Ta.star <- Inf
#     #Ta.star.laws <- Inf
#     Ta.star.pis2 <- Inf
#     ta.star.sabha <- 0
#     ta.star.laws <- 0
#     ta.star.ihw <- 0
#     ta.star.ihw.null <- 0
#     #===========================================
#     #====== Fix Reigion
#     #===========================================
#     # save(X, Tm,pis.hat,qhat,
#     #     file=paste0("Result/Simulation1D_h/Tmp/no_rep_mag_",magnitude,"_mu_",mu_type,"_cov_",Cov_type,"_seed_",seed,".RData"))
#     
#     # Detect Neighbor
#     hh.seq <- rep(hh,m)
#     Neigh_Detect_res <- Neigh_Detect(hh = hh.seq,
#                                      X = X, 
#                                      Dist = Dist.p,
#                                      Sigma.eps = Sigma.eps.est,
#                                      detect.m = detect.m)
#     Ta <- Neigh_Detect_res$Ta
#     Va <- Neigh_Detect_res$Va
#     VmVa.cov <- Neigh_Detect_res$VmVa.cov
#     ind <- Neigh_Detect_res$ind
#     
#     #== round numerical values to ensures the correctness of taking maximum
#     Tm <- round(Tm,dig)
#     Ta <- round(Ta,dig)
#     Va <- round(Va,dig)
#     VmVa.cov <- round(VmVa.cov,dig)
#     
#     #===========================================
#     #== Run 2D Selection
#     #===========================================
#     print(paste("seed:",seed,"BH","hh:",hh))
#     res.2D <- TwoDSMT::Spatial_Detect_exact_grp_BH_down(Tm, Ta, Va, 
#                                                         VmVa.cov, ind,
#                                                         q, max.rej,
#                                                         Tm.star = Tm.star,
#                                                         Ta.star = Ta.star,
#                                                         const = const,
#                                                         seed = seed)
#     selected.2D <- res.2D$selected
#     tm <- res.2D$tm0
#     ta <- res.2D$ta0
#     
#     print(paste("seed:",seed,"BH_pis2","hh:",hh))
#     res.2D.pis2 <- TwoDSMT::Spatial_Detect_exact_grp_BH_down(Tm, Ta, Va,
#                                                              VmVa.cov, ind,
#                                                              q, 
#                                                              max.rej = max.rej.pis2,
#                                                              pis = pis.hata,
#                                                              Tm.star = Tm.star.pis2,
#                                                              Ta.star = Ta.star.pis2,
#                                                              const = const,
#                                                              seed = seed)
#     selected.2D.pis2 <- res.2D.pis2$selected
#     tm <- res.2D.pis2$tm0
#     ta <- res.2D.pis2$ta0
#     
#     
#     print(paste("seed:",seed,"BH_laws","hh:",hh))
#     res.2D.laws <- TwoDSMT::Spatial_Detect_exact_BH_down_reTm_reTa(Tm, Ta, Va,
#                                                                    VmVa.cov, ind,
#                                                                    q, 
#                                                                    max.rej = max.rej.laws,
#                                                                    pis = 1-pis.hat,
#                                                                    pws.tm.star = tm.star.laws,
#                                                                    pws.ta.star = ta.star.laws,
#                                                                    const = const,
#                                                                    ws.fun = ws.laws.fun,
#                                                                    n.group.max = n.group.max,
#                                                                    seed = seed)
#     selected.2D.laws <- res.2D.laws$selected
#     tm <- res.2D.laws$tm0
#     ta <- res.2D.laws$ta0
#     
#     print(paste("seed:",seed,"BH_sabha","hh:",hh))
#     res.2D.sabha <- TwoDSMT::Spatial_Detect_exact_BH_down_reTm_reTa(Tm, Ta, Va,
#                                                                     VmVa.cov, ind,
#                                                                     q, 
#                                                                     max.rej = max.rej.sabha,
#                                                                     pis = qhat,
#                                                                     pws.tm.star = tm.star.sabha,
#                                                                     pws.ta.star = ta.star.sabha,
#                                                                     const = const,
#                                                                     ws.fun = ws.sabha.fun,
#                                                                     n.group.max = n.group.max,
#                                                                     seed = seed)
#     selected.2D.sabha <- res.2D.sabha$selected
#     tm <- res.2D.sabha$tm0
#     ta <- res.2D.sabha$ta0
#     
#     print(paste("seed:",seed,"BH_ihw","hh:",hh))
#     res.2D.ihw <- TwoDSMT::Spatial_Detect_exact_BH_down_reTm_reTa(Tm, Ta, Va,
#                                                                   VmVa.cov, ind,
#                                                                   q, 
#                                                                   max.rej = max.rej.ihw,
#                                                                   pis = rep(1,m),
#                                                                   pws.tm.star = tm.star.ihw,
#                                                                   pws.ta.star = ta.star.ihw,
#                                                                   const = const,
#                                                                   ws = ihw.ws,
#                                                                   n.group.max = n.group.max,
#                                                                   seed = seed)
#     selected.2D.ihw <- res.2D.ihw$selected
#     tm <- res.2D.ihw$tm0
#     ta <- res.2D.ihw$ta0
#     
#     
#     print(paste("seed:",seed,"BH_ihw.null","hh:",hh))
#     res.2D.ihw.null <- TwoDSMT::Spatial_Detect_exact_BH_down_reTm_reTa(Tm, Ta, Va,
#                                                                        VmVa.cov, ind,
#                                                                        q, 
#                                                                        max.rej = max.rej.ihw.null,
#                                                                        pis = pis.ihw,
#                                                                        pws.tm.star = tm.star.ihw.null,
#                                                                        pws.ta.star = ta.star.ihw.null,
#                                                                        const = const,
#                                                                        ws = ihw.null.ws,
#                                                                        n.group.max = n.group.max,
#                                                                        seed = seed)
#     selected.2D.ihw.null <- res.2D.ihw.null$selected
#     tm <- res.2D.ihw.null$tm0
#     ta <- res.2D.ihw.null$ta0
#     
#     fdp.nopis.res <- c(fdp.nopis.res, fdp(selected.2D,mu))
#     #fdp.pis.res <- c(fdp.pis.res,fdp(selected.2D.pis,mu))
#     fdp.pis2.res <- c(fdp.pis2.res,fdp(selected.2D.pis2,mu))
#     fdp.laws.res <- c(fdp.laws.res,fdp(selected.2D.laws,mu))
#     fdp.sabha.res <- c(fdp.sabha.res,fdp(selected.2D.sabha,mu))
#     fdp.ihw.res <- c(fdp.ihw.res,fdp(selected.2D.ihw,mu))
#     fdp.ihw.null.res <- c(fdp.ihw.null.res,fdp(selected.2D.ihw.null,mu))
#     
#     pow.nopis.res <- c(pow.nopis.res, Pow(selected.2D,mu))
#     #pow.pis.res <- c(pow.pis.res,Pow(selected.2D.pis,mu))
#     pow.pis2.res <- c(pow.pis2.res,Pow(selected.2D.pis2,mu))
#     pow.laws.res <- c(pow.laws.res,Pow(selected.2D.laws,mu))
#     pow.sabha.res <- c(pow.sabha.res,Pow(selected.2D.sabha,mu))
#     pow.ihw.res <- c(pow.ihw.res,Pow(selected.2D.ihw,mu))
#     pow.ihw.null.res <- c(pow.ihw.null.res,Pow(selected.2D.ihw.null,mu))
#     
#   }
#   print(rbind(fdp.nopis.res, fdp.pis2.res, fdp.laws.res, fdp.sabha.res, fdp.ihw.res, fdp.ihw.null.res))
#   print(rbind(pow.nopis.res, pow.pis2.res, pow.laws.res, pow.sabha.res, pow.ihw.res, pow.ihw.null.res))
#   return(list(fdp.res = fdp.res,
#               pow.res = pow.res,
#               fdp.nopis.res = fdp.nopis.res,
#               pow.nopis.res = pow.nopis.res,
#               fdp.pis2.res = fdp.pis2.res,
#               pow.pis2.res = pow.pis2.res,
#               fdp.laws.res = fdp.laws.res,
#               pow.laws.res = pow.laws.res,
#               fdp.sabha.res = fdp.sabha.res,
#               pow.sabha.res = pow.sabha.res,
#               fdp.ihw.res = fdp.ihw.res,
#               pow.ihw.res = pow.ihw.res,
#               fdp.ihw.null.res = fdp.ihw.null.res,
#               pow.ihw.null.res = pow.ihw.null.res,
#               seed = seed))
# }


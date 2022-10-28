#' The Robust Form
#' 
one_step_1D_Robust <- function(h, detect.m = "rad.h", seed = 0,
                          data.type = "mvnorm",
                          mu = mu,
                          const = 1,
                          nbins.ihw = 5,
                          n.group.max = 5,
                          tau.t1 = 1,
                          tau.t2 = 1,
                          #ws.fun = function(x){1/x}, # It should depends on the method we use
                          ...){
  set.seed(seed)
  print(paste("seed",seed))
  # Current DataType
  # "mvnorm", "Circle","SquareSpline","CircleMix"
  
  pluszero <- mu>0
  #X <- user_mvrnorm(n = n, mu = mu, 
  #                  Sigma = Sigma.eps.p)
  X <- MASS::mvrnorm(n = n, mu = mu, 
                     Sigma = Sigma.eps.p)
  X <- matrix(X, nrow = n)
  
  #==== Initilize weight function (For LAWS and SABHA respectively)
  ws.laws.fun <- function(x){(1-x)/x}
  ws.sabha.fun <- function(x){1/x}
  
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
  T1 <- apply(X,2,function(x){sum(x)/sqrt(n)})/sgm
  T1 <- round(T1,dig)
  #=== Without spatial Info: BH
  #=== (h=0)
  p.value <- 1 - pnorm(T1)
  BH.res <- qvalue(p.value,pi0 = 1)
  BH.selected <- which(BH.res$qvalues<=q)
  #fdp(selected.BH,mu);Pow(selected.BH,mu)
  
  #================================
  #==== Without spatial Info: One D
  #================================
  res.1D <- OneD_Detect(T1, q, pis = NULL,
                        const = const,
                        tau.t1 = tau.t1)
  T1.star <- res.1D$t1.min 
  max.rej <- res.1D$max.rej
  fdp.nopis.res <- fdp(res.1D$selected,mu)
  pow.nopis.res <- Pow(res.1D$selected,mu)
  NP.nopis.res <- length(res.1D$selected)
  
  #==== Calculate pis (LAWS)
  bh.th<-bh.func(p.value, 0.9)$th
  if(bh.th == 1){
    pis.hat <- rep(1,m)
  }else{
    pis.hat<-pis_1D.func.ker.reg(T1, point, tau=bh.th)
  }
  
  #==== Calculate qhat (SABHA)
  tau = 0.5;eps = 0.1; TV_bd = 2
  alpha_ADMM = 10^2; beta = 10^3; eta = 5; max_iters = 5000; converge_thr = 1e-4 # parameters for ADMM
  ADMM_params = c(alpha_ADMM,beta,eta,max_iters,converge_thr)
  
  qhat = Solve_q_TV_1dim(p.value, tau, eps, TV_bd, ADMM_params)
  
  #==== Calculate pis2 (Storey)
  lambda.test <- 0.5#median(T1)#quantile(T1,0.9)#
  pis.hat2.tmp<- min(1,
                     mean(T1<lambda.test)/pnorm(lambda.test))
  pis.hat2 <- rep(pis.hat2.tmp,m)
  #pis.hat2 <- (T1<lambda.test)/pnorm(lambda.test)
  
  #==============================================
  #==== Run Other Algorithms: Storey,LAWS, SABHA, IHW
  #==============================================
  #--- Storey
  print(paste0("seed",seed," Start Storey..........."))
  Storey.res <- qvalue(p.value,pi0 = pis.hat2.tmp)
  Storey.selected <- which(Storey.res$qvalues<=q)
  
  #--- LAWS
  print(paste0("seed",seed," Start LAWS..........."))
  laws.res <- law.func(pvs=p.value, pis.hat, q)
  laws.selected <- which(laws.res$de==1)
  
  #--- SABHA
  print(paste0("seed",seed," Start SABHA..........."))
  SABHA_method = function(pvals, qhat, alpha, tau){
    pvals[pvals>tau] = Inf
    khat=max(c(0,which(sort(qhat*pvals)<=alpha*(1:length(pvals))/length(pvals))))
    which(qhat*pvals<=alpha*khat/length(pvals))
  }
  
  sab.selected <- SABHA_method(p.value, qhat, q, tau)
  
  #--- IHW: We use the location of IHW as covariate
  # (Without null proportion estimation)
  print(paste0("seed",seed," Start IHW..........."))
  ihw.res <- ihw(p.value, point[,1], q,nbins = nbins.ihw) 
  ihw.selected <- which(ihw.res@df$adj_pvalue<q)
  ihw.ws <- ihw.res@df$weight # Extract weights
  
  #--- IHW: We use the location of IHW as covariate
  # (With null proportion estimation)
  print(paste0("seed",seed," Start IHW..........."))
  ihw.null.res <- ihw(p.value, point[,1], q,nbins = nbins.ihw, null_proportion=T) 
  ihw.null.selected <- which(ihw.null.res@df$adj_pvalue<q)
  ihw.null.ws <- ihw.null.res@df$weight # Extract weights
  
  pis.ihw <- ihw.ws/ihw.null.ws
  pis.ihw[is.na(pis.ihw)] <- 1 # NA, represents the null proportion is 1. Two weights are both 0.
  
  #--- Calculate fdp and power
  fdp.res <- c(fdp(BH.selected,mu), fdp(Storey.selected, mu),
               fdp(laws.selected,mu), fdp(sab.selected, mu), 
               fdp(ihw.selected, mu),fdp(ihw.null.selected, mu))
  pow.res <- c(Pow(BH.selected,mu), Pow(Storey.selected, mu),
               Pow(laws.selected,mu), Pow(sab.selected, mu), 
               Pow(ihw.selected, mu),Pow(ihw.null.selected, mu))
  NP.res <- c(length(BH.selected), length(Storey.selected),
              length(laws.selected),length(sab.selected),
              length(ihw.selected),length(ihw.null.selected))
  names(fdp.res) <- c("BH","Storey","LAWS","SABHA","IHW","IHW.null")
  #== Round numerical values to ensures the correctness of taking maximum
  dig <- 7
  
  pis.hat <- round(pis.hat,dig)
  pis.hat2 <- round(pis.hat2,dig)
  qhat <- round(qhat,dig)
  
  #==== Initilize for 2D detection ==========
  #===========================================
  #==== Without spatial Info: One D (with pis)
  #==== LAWS based
  #===========================================
  res.1D.laws <- OneD_Detect_w(T1, q, pis = 1-pis.hat,
                               ws.fun = ws.laws.fun,
                               const = const,
                               tau.t1 = tau.t1)
  t1.star.laws <- res.1D.laws$t1.min
  max.rej.laws <- res.1D.laws$max.rej
  fdp.laws.res <- fdp(res.1D.laws$selected,mu)
  pow.laws.res <- Pow(res.1D.laws$selected,mu)
  NP.laws.res <- length(res.1D.laws$selected)
  
  #===========================================
  #==== Without spatial Info: One D (with qhat)
  #==== SABHA based
  #===========================================
  res.1D.sabha <- OneD_Detect_w(T1, q, pis = qhat,
                                const = const,
                                ws.fun = ws.sabha.fun,
                                tau.t1 = tau.t1)
  t1.star.sabha <- res.1D.sabha$t1.min
  max.rej.sabha <- res.1D.sabha$max.rej
  fdp.sabha.res <- fdp(res.1D.sabha$selected,mu)
  pow.sabha.res <- Pow(res.1D.sabha$selected,mu)
  NP.sabha.res <- length(res.1D.sabha$selected)
  
  
  #===========================================
  #==== Without spatial Info: One D (with IHW)
  #==== IHW based (Without null proportion)
  #===========================================
  res.1D.ihw <- OneD_Detect_w(T1, q,
                              ws = ihw.ws,
                              const = const,
                              tau.t1 = tau.t1)
  t1.star.ihw <- res.1D.ihw$t1.min
  max.rej.ihw <- res.1D.ihw$max.rej
  fdp.ihw.res <- fdp(res.1D.ihw$selected,mu)
  pow.ihw.res <- Pow(res.1D.ihw$selected,mu)
  NP.ihw.res <- length(res.1D.ihw$selected)
  
  #===========================================
  #==== Without spatial Info: One D (with IHW)
  #==== IHW based (With null proportion)
  #===========================================
  res.1D.ihw.null <- OneD_Detect_w(T1, q, pis = pis.ihw,
                                   ws = ihw.null.ws,
                                   const = const,
                                   tau.t1 = tau.t1)
  t1.star.ihw.null <- res.1D.ihw.null$t1.min
  max.rej.ihw.null <- res.1D.ihw.null$max.rej
  fdp.ihw.null.res <- fdp(res.1D.ihw.null$selected,mu)
  pow.ihw.null.res <- Pow(res.1D.ihw.null$selected,mu)
  NP.ihw.null.res <- length(res.1D.ihw.null$selected)
  
  #===========================================
  #==== Without spatial Info: One D (with pis2)
  #---- The results from OneD_Detect and OneD_Detect_w differ
  #===========================================
  res.1D.pis2 <- OneD_Detect(T1, q, pis = pis.hat2,
                             const = const,
                             tau.t1 = tau.t1)
  T1.star.pis2 <- res.1D.pis2$t1.min
  max.rej.pis2 <- res.1D.pis2$max.rej
  fdp.pis2.res <- fdp(res.1D.pis2$selected,mu)
  pow.pis2.res <- Pow(res.1D.pis2$selected,mu)
  NP.pis2.res <- length(res.1D.pis2$selected)
  
  #=== With Spatial Info
  save(X,qhat, 
       file=paste0("Result/Simulation1D_Robust/Tmp/mag_",
                   magnitude,"_mu_",mu_type,"_cov_",Cov_type,"_seed_",seed,".RData")
  )
  for(est_cov_type in c("Weak","Medium","Strong")){
   Sigma.eps.est <-  Init_Setting_Cov(est_cov_type)
  for(ii in 1:length(h)){
    hh <- h[ii]
    print(paste("seed:",seed,"hh:",hh))
    
    #==== Set pre-chosen T1, T2
    # T1.star
    T1.star <- round(T1.star,dig)
    #T1.star.laws <- round(T1.star.laws,dig)
    T1.star.pis2 <- round(T1.star.pis2,dig)
    # No need to be round because pws hasn't been round
    #t1.star.sabha <- round(t1.star.sabha,dig)
    
    # T2.star
    T2.star <- Inf
    #T2.star.laws <- Inf
    T2.star.pis2 <- Inf
    t2.star.sabha <- 0
    t2.star.laws <- 0
    t2.star.ihw <- 0
    t2.star.ihw.null <- 0
    #===========================================
    #====== Fix Reigion
    #===========================================
    # save(X, T1,pis.hat,qhat,
    #     file=paste0("Result/Simulation1D_h/Tmp/no_rep_mag_",magnitude,"_mu_",mu_type,"_cov_",Cov_type,"_seed_",seed,".RData"))
    
    # Detect Neighbor
    hh.seq <- rep(hh,m)
    Neigh_Detect_res <- Neigh_Detect(hh = hh.seq,
                                     X = X, 
                                     Dist = Dist.p,
                                     Sigma.eps = Sigma.eps.est,
                                     detect.m = detect.m)
    T2 <- Neigh_Detect_res$T2
    V2 <- Neigh_Detect_res$V2
    V1V2.cov <- Neigh_Detect_res$V1V2.cov
    ind <- Neigh_Detect_res$ind
    
    #== round numerical values to ensures the correctness of taking maximum
    T1 <- round(T1,dig)
    T2 <- round(T2,dig)
    V2 <- round(V2,dig)
    V1V2.cov <- round(V1V2.cov,dig)
    
    #===========================================
    #== Run 2D Selection
    #===========================================
    print(paste("seed:",seed,"BH","hh:",hh))
    res.2D <- Spatial_Detect_exact_grp_BH_down(T1, T2, V2, 
                                               V1V2.cov, ind,
                                               q, max.rej,
                                               T1.star = T1.star,
                                               T2.star = T2.star,
                                               const = const,
                                               seed = seed,
                                               tau.t1 = tau.t1,
                                               tau.t2 = tau.t2)
    selected.2D <- res.2D$selected
    t1 <- res.2D$t10
    t2 <- res.2D$t20
    
    print(paste("seed:",seed,"BH_pis2","hh:",hh))
    res.2D.pis2 <- Spatial_Detect_exact_grp_BH_down(T1, T2, V2,
                                                    V1V2.cov, ind,
                                                    q, 
                                                    max.rej = max.rej.pis2,
                                                    pis = pis.hat2,
                                                    T1.star = T1.star.pis2,
                                                    T2.star = T2.star.pis2,
                                                    const = const,
                                                    seed = seed,
                                                    tau.t1 = tau.t1,
                                                    tau.t2 = tau.t2)
    selected.2D.pis2 <- res.2D.pis2$selected
    t1 <- res.2D.pis2$t10
    t2 <- res.2D.pis2$t20
    
    
    print(paste("seed:",seed,"BH_laws","hh:",hh))
    res.2D.laws <- Spatial_Detect_exact_BH_down_reT1_reT2(T1, T2, V2,
                                                          V1V2.cov, ind,
                                                          q, 
                                                          max.rej = max.rej.laws,
                                                          pis = 1-pis.hat,
                                                          pws.t1.star = t1.star.laws,
                                                          pws.t2.star = t2.star.laws,
                                                          const = const,
                                                          ws.fun = ws.laws.fun,
                                                          seed = seed,
                                                          tau.t1 = tau.t1,
                                                          tau.t2 = tau.t2)
    selected.2D.laws <- res.2D.laws$selected
    t1 <- res.2D.laws$t10
    t2 <- res.2D.laws$t20
    
    print(paste("seed:",seed,"BH_sabha","hh:",hh))
    res.2D.sabha <- Spatial_Detect_exact_BH_down_reT1_reT2(T1, T2, V2,
                                                           V1V2.cov, ind,
                                                           q, 
                                                           max.rej = max.rej.sabha,
                                                           pis = qhat,
                                                           pws.t1.star = t1.star.sabha,
                                                           pws.t2.star = t2.star.sabha,
                                                           const = const,
                                                           ws.fun = ws.sabha.fun,
                                                           n.group.max = n.group.max,
                                                           seed = seed,
                                                           tau.t1 = tau.t1,
                                                           tau.t2 = tau.t2)
    selected.2D.sabha <- res.2D.sabha$selected
    t1 <- res.2D.sabha$t10
    t2 <- res.2D.sabha$t20
    
    print(paste("seed:",seed,"BH_ihw","hh:",hh))
    res.2D.ihw <- Spatial_Detect_exact_BH_down_reT1_reT2(T1, T2, V2,
                                                         V1V2.cov, ind,
                                                         q, 
                                                         max.rej = max.rej.ihw,
                                                         pis = rep(1,m),
                                                         pws.t1.star = t1.star.ihw,
                                                         pws.t2.star = t2.star.ihw,
                                                         const = const,
                                                         ws = ihw.ws,
                                                         n.group.max = n.group.max,
                                                         seed = seed,
                                                         tau.t1 = tau.t1,
                                                         tau.t2 = tau.t2)
    selected.2D.ihw <- res.2D.ihw$selected
    t1 <- res.2D.ihw$t10
    t2 <- res.2D.ihw$t20
    
    
    print(paste("seed:",seed,"BH_ihw.null","hh:",hh))
    res.2D.ihw.null <- Spatial_Detect_exact_BH_down_reT1_reT2(T1, T2, V2,
                                                              V1V2.cov, ind,
                                                              q, 
                                                              max.rej = max.rej.ihw.null,
                                                              pis = pis.ihw,
                                                              pws.t1.star = t1.star.ihw.null,
                                                              pws.t2.star = t2.star.ihw.null,
                                                              const = const,
                                                              ws = ihw.null.ws,
                                                              n.group.max = n.group.max,
                                                              seed = seed,
                                                              tau.t1 = tau.t1,
                                                              tau.t2 = tau.t2)
    selected.2D.ihw.null <- res.2D.ihw.null$selected
    t1 <- res.2D.ihw.null$t10
    t2 <- res.2D.ihw.null$t20
    
    fdp.nopis.res <- c(fdp.nopis.res, fdp(selected.2D,mu))
    #fdp.pis.res <- c(fdp.pis.res,fdp(selected.2D.pis,mu))
    fdp.pis2.res <- c(fdp.pis2.res,fdp(selected.2D.pis2,mu))
    fdp.laws.res <- c(fdp.laws.res,fdp(selected.2D.laws,mu))
    fdp.sabha.res <- c(fdp.sabha.res,fdp(selected.2D.sabha,mu))
    fdp.ihw.res <- c(fdp.ihw.res,fdp(selected.2D.ihw,mu))
    fdp.ihw.null.res <- c(fdp.ihw.null.res,fdp(selected.2D.ihw.null,mu))
    
    pow.nopis.res <- c(pow.nopis.res, Pow(selected.2D,mu))
    #pow.pis.res <- c(pow.pis.res,Pow(selected.2D.pis,mu))
    pow.pis2.res <- c(pow.pis2.res,Pow(selected.2D.pis2,mu))
    pow.laws.res <- c(pow.laws.res,Pow(selected.2D.laws,mu))
    pow.sabha.res <- c(pow.sabha.res,Pow(selected.2D.sabha,mu))
    pow.ihw.res <- c(pow.ihw.res,Pow(selected.2D.ihw,mu))
    pow.ihw.null.res <- c(pow.ihw.null.res,Pow(selected.2D.ihw.null,mu))
    
  }
  }
  print(rbind(fdp.nopis.res, fdp.pis2.res, fdp.laws.res, fdp.sabha.res, fdp.ihw.res, fdp.ihw.null.res))
  print(rbind(pow.nopis.res, pow.pis2.res, pow.laws.res, pow.sabha.res, pow.ihw.res, pow.ihw.null.res))
  return(list(fdp.res = fdp.res,
              pow.res = pow.res,
              fdp.nopis.res = fdp.nopis.res,
              pow.nopis.res = pow.nopis.res,
              fdp.pis2.res = fdp.pis2.res,
              pow.pis2.res = pow.pis2.res,
              fdp.laws.res = fdp.laws.res,
              pow.laws.res = pow.laws.res,
              fdp.sabha.res = fdp.sabha.res,
              pow.sabha.res = pow.sabha.res,
              fdp.ihw.res = fdp.ihw.res,
              pow.ihw.res = pow.ihw.res,
              fdp.ihw.null.res = fdp.ihw.null.res,
              pow.ihw.null.res = pow.ihw.null.res,
              seed = seed))
}
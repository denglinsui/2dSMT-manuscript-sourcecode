mu.fix.gen <- function(point, method, magnitude,mu.mean,Sigma.mu){
  m <- nrow(point)
  if(method == "uc.unif"){
    location <- (point[,1]-1/2)^2+(point[,2]-1/2)^2<(1/4)^2
    pis <- numeric(m)
    pis[location] <- 0.9
    pis[!location] <- 0.01
    mu <- magnitude * rbinom(n=length(pis),size=1,prob=pis)
    pis <- as.numeric(mu<=0)
  }
  
  if(method == "uc.spline"){
    #=== Init Basis
    n.break <- 3
    n.order <- 4
    n.end <- n.order + n.break -2
    sp.basis <- create.bspline.basis(c(0.25,0.75), dropind=c(1,n.end),
                                     breaks = seq(0.25, 0.75, length.out = n.break), 
                                     norder = n.order)
    
    point.basis <- point %>%  as.data.frame() %>% mutate(mu.x = 0) %>% mutate(mu.y = 0)
    colnames(point.basis) <- c("x","y","mu.x","mu.y")
    
    #=== Basis product
    point.basis[point.basis$x>=0.25 & point.basis$x<=0.75,] <- 
      point.basis %>%
      filter(x>=0.25 & x<=0.75) %>% 
      mutate(mu.x = sqrt(n.end)/sqrt(n.end-2) * rowSums(eval.basis(x,basisobj = sp.basis)))
    
    point.basis[point.basis$y>=0.25 & point.basis$y<=0.75,] <- 
      point.basis %>%
      filter(y>=0.25 & y<=0.75) %>% 
      mutate(mu.y = sqrt(n.end)/sqrt(n.end-2) * rowSums(eval.basis(y,basisobj = sp.basis)))
    #sqrt(n.end)/sqrt(n.end-2) is for normalization
    point.basis <- point.basis %>% mutate(mu = mu.x*mu.y) %>% mutate(mu = mu)
    
    mu <- magnitude * point.basis$mu
    pis <- as.numeric(mu<=0)
    
  }
  
  if(method == "mvnorm"){
    mu <- magnitude*MASS::mvrnorm(n = 1, mu = rep(mu.mean,m), Sigma = Sigma.mu)
    pis <- as.numeric(mu<=0)
  }
  
  if(method == "mixture"){
    m <- nrow(point)
    pis <- numeric(m)
    location <- (point[,1]-1/2)^2+(point[,2]-1/2)^2<(1/4)^2
    pis[location] <- 0.9
    pis[-location] <- 0.1
    
    location <- rbinom(m,1,pis)
    mu <- location*magnitude
    pis <- 1-pis
  }
  return(list(mu=mu,
              pis=pis))
}


#==================
# Simulation Step
# 1 dimension
# Varying h
#==================
one_step_1D_h <- function(h, detect.m = "rad.h", seed = 0,
                          data.type = "mvnorm",
                          mu = mu,
                          const = 1,
                          nbins.ihw = 5,
                          n.group.max = 5,
                          tau.tm = 1,tau.ta =1,
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
  Tm <- apply(X,2,function(x){sum(x)/sqrt(n)})/sgm
  Tm <- round(Tm,dig)
  #=== Without spatial Info: BH
  #=== (h=0)
  p.value <- 1 - pnorm(Tm)
  BH.res <- qvalue(p.value,pi0 = 1)
  BH.selected <- which(BH.res$qvalues<=q)
  #fdp(selected.BH,mu);Pow(selected.BH,mu)
  
  #================================
  #==== Without spatial Info: One D
  #================================
  res.1D <- OneD_Detect(Tm, q, pis = NULL,
                        const = const,
                        tau.tm = tau.tm)
  Tm.star <- res.1D$tm.min 
  max.rej <- res.1D$max.rej
  fdp.nopis.res <- fdp(res.1D$selected,mu)
  pow.nopis.res <- Pow(res.1D$selected,mu)
  NP.nopis.res <- length(res.1D$selected)
  
  #==== Calculate pis (LAWS)
  print(paste0("seed:", seed, "========== Calculating pis =============="))
  bh.th<-bh.func(p.value, 0.9)$th
  if(bh.th == 1){
    pis.hat <- rep(1,m)
  }else{
    pis.hat<-pis_1D.func.ker.reg(Tm, point, tau=bh.th)
  }
  
  #==== Calculate qhat (SABHA)
  print(paste0("seed:", seed, "========== Calculating qhat =============="))
  tau = 0.5;eps = 0.1; TV_bd = 2
  alpha_ADMM = 10^2; beta = 10^3; eta = 5; max_iters = 5000; converge_thr = 1e-4 # parameters for ADMM
  ADMM_params = c(alpha_ADMM,beta,eta,max_iters,converge_thr)
  
  qhat = Solve_q_TV_1dim(p.value, tau, eps, TV_bd, ADMM_params)
  
  #==== Calculate pis2 (Storey)
  lambda.test <- 0.5#median(Tm)#quantile(Tm,0.9)#
  pis.hata.tmp<- min(1,
                     mean(Tm<lambda.test)/pnorm(lambda.test))
  pis.hata <- rep(pis.hata.tmp,m)
  #pis.hata <- (Tm<lambda.test)/pnorm(lambda.test)
  
  #==============================================
  #==== Run Other Algorithms: Storey,LAWS, SABHA, IHW
  #==============================================
  #--- Storey
  print(paste0("seed",seed," Start Storey..........."))
  Storey.res <- qvalue(p.value,pi0 = pis.hata.tmp)
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
  pis.hata <- round(pis.hata,dig)
  qhat <- round(qhat,dig)
  
  #==== Initilize for 2D detection ==========
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
  fdp.laws.res <- fdp(res.1D.laws$selected,mu)
  pow.laws.res <- Pow(res.1D.laws$selected,mu)
  NP.laws.res <- length(res.1D.laws$selected)
  
  #===========================================
  #==== Without spatial Info: One D (with qhat)
  #==== SABHA based
  #===========================================
  res.1D.sabha <- OneD_Detect_w(Tm, q, pis = qhat,
                                const = const,
                                ws.fun = ws.sabha.fun,
                                tau.tm = tau.tm)
  tm.star.sabha <- res.1D.sabha$tm.min
  max.rej.sabha <- res.1D.sabha$max.rej
  fdp.sabha.res <- fdp(res.1D.sabha$selected,mu)
  pow.sabha.res <- Pow(res.1D.sabha$selected,mu)
  NP.sabha.res <- length(res.1D.sabha$selected)
  
  
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
  fdp.ihw.res <- fdp(res.1D.ihw$selected,mu)
  pow.ihw.res <- Pow(res.1D.ihw$selected,mu)
  NP.ihw.res <- length(res.1D.ihw$selected)
  
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
  fdp.ihw.null.res <- fdp(res.1D.ihw.null$selected,mu)
  pow.ihw.null.res <- Pow(res.1D.ihw.null$selected,mu)
  NP.ihw.null.res <- length(res.1D.ihw.null$selected)
  
  #===========================================
  #==== Without spatial Info: One D (with pis2)
  #---- The results from OneD_Detect and OneD_Detect_w differ
  #===========================================
  res.1D.pis2 <- OneD_Detect(Tm, q, pis = pis.hata,
                             const = const,
                             tau.tm = tau.tm)
  Tm.star.pis2 <- res.1D.pis2$tm.min
  max.rej.pis2 <- res.1D.pis2$max.rej
  fdp.pis2.res <- fdp(res.1D.pis2$selected,mu)
  pow.pis2.res <- Pow(res.1D.pis2$selected,mu)
  NP.pis2.res <- length(res.1D.pis2$selected)
  
  #=== With Spatial Info
  for(ii in 1:length(h)){
    hh <- h[ii]
    print(paste("seed:",seed,"hh:",hh))
    
    #==== Set pre-chosen Tm, Ta
    # Tm.star
    Tm.star <- round(Tm.star,dig)
    #Tm.star.laws <- round(Tm.star.laws,dig)
    Tm.star.pis2 <- round(Tm.star.pis2,dig)
    # No need to be round because pws hasn't been round
    #tm.star.sabha <- round(tm.star.sabha,dig)
    
    # Ta.star
    Ta.star <- Inf
    #Ta.star.laws <- Inf
    Ta.star.pis2 <- Inf
    ta.star.sabha <- 0
    ta.star.laws <- 0
    ta.star.ihw <- 0
    ta.star.ihw.null <- 0
    #===========================================
    #====== Fix Reigion
    #===========================================
    # save(X, Tm,pis.hat,qhat,
    #     file=paste0("Result/Simulation1D_h/Tmp/no_rep_mag_",magnitude,"_mu_",mu_type,"_cov_",Cov_type,"_seed_",seed,".RData"))
    
    # Detect Neighbor
    hh.seq <- rep(hh,m)
    Neigh_Detect_res <- Neigh_Detect(hh = hh.seq,
                                     X = X, 
                                     Dist = Dist.p,
                                     Sigma.eps = Sigma.eps.est,
                                     detect.m = detect.m)
    Ta <- Neigh_Detect_res$Ta
    Va <- Neigh_Detect_res$Va
    VmVa.cov <- Neigh_Detect_res$VmVa.cov
    ind <- Neigh_Detect_res$ind
    
    #== round numerical values to ensures the correctness of taking maximum
    Tm <- round(Tm,dig)
    Ta <- round(Ta,dig)
    Va <- round(Va,dig)
    VmVa.cov <- round(VmVa.cov,dig)
    
    #===========================================
    #== Run 2D Selection
    #===========================================
    print(paste("seed:",seed,"BH","hh:",hh))
    res.2D <- Spatial_Detect_exact_grp_BH_down(Tm, Ta, Va, 
                                               VmVa.cov, ind,
                                               q, max.rej,
                                               Tm.star = Tm.star,
                                               Ta.star = Ta.star,
                                               const = const,
                                               seed = seed)
    selected.2D <- res.2D$selected
    tm <- res.2D$tm0
    ta <- res.2D$ta0
    
    print(paste("seed:",seed,"BH_pis2","hh:",hh))
    res.2D.pis2 <- Spatial_Detect_exact_grp_BH_down(Tm, Ta, Va,
                                                    VmVa.cov, ind,
                                                    q, 
                                                    max.rej = max.rej.pis2,
                                                    pis = pis.hata,
                                                    Tm.star = Tm.star.pis2,
                                                    Ta.star = Ta.star.pis2,
                                                    const = const,
                                                    seed = seed)
    selected.2D.pis2 <- res.2D.pis2$selected
    tm <- res.2D.pis2$tm0
    ta <- res.2D.pis2$ta0
    
    
    print(paste("seed:",seed,"BH_laws","hh:",hh))
    res.2D.laws <- Spatial_Detect_exact_BH_down_reTm_reTa(Tm, Ta, Va,
                                                          VmVa.cov, ind,
                                                          q, 
                                                          max.rej = max.rej.laws,
                                                          pis = 1-pis.hat,
                                                          pws.tm.star = tm.star.laws,
                                                          pws.ta.star = ta.star.laws,
                                                          const = const,
                                                          ws.fun = ws.laws.fun,
                                                          n.group.max = n.group.max,
                                                          seed = seed)
    selected.2D.laws <- res.2D.laws$selected
    tm <- res.2D.laws$tm0
    ta <- res.2D.laws$ta0
    
    print(paste("seed:",seed,"BH_sabha","hh:",hh))
    res.2D.sabha <- Spatial_Detect_exact_BH_down_reTm_reTa(Tm, Ta, Va,
                                                           VmVa.cov, ind,
                                                           q, 
                                                           max.rej = max.rej.sabha,
                                                           pis = qhat,
                                                           pws.tm.star = tm.star.sabha,
                                                           pws.ta.star = ta.star.sabha,
                                                           const = const,
                                                           ws.fun = ws.sabha.fun,
                                                           n.group.max = n.group.max,
                                                           seed = seed)
    selected.2D.sabha <- res.2D.sabha$selected
    tm <- res.2D.sabha$tm0
    ta <- res.2D.sabha$ta0
    
    print(paste("seed:",seed,"BH_ihw","hh:",hh))
    res.2D.ihw <- Spatial_Detect_exact_BH_down_reTm_reTa(Tm, Ta, Va,
                                                         VmVa.cov, ind,
                                                         q, 
                                                         max.rej = max.rej.ihw,
                                                         pis = rep(1,m),
                                                         pws.tm.star = tm.star.ihw,
                                                         pws.ta.star = ta.star.ihw,
                                                         const = const,
                                                         ws = ihw.ws,
                                                         n.group.max = n.group.max,
                                                         seed = seed)
    selected.2D.ihw <- res.2D.ihw$selected
    tm <- res.2D.ihw$tm0
    ta <- res.2D.ihw$ta0
    
    
    print(paste("seed:",seed,"BH_ihw.null","hh:",hh))
    res.2D.ihw.null <- Spatial_Detect_exact_BH_down_reTm_reTa(Tm, Ta, Va,
                                                              VmVa.cov, ind,
                                                              q, 
                                                              max.rej = max.rej.ihw.null,
                                                              pis = pis.ihw,
                                                              pws.tm.star = tm.star.ihw.null,
                                                              pws.ta.star = ta.star.ihw.null,
                                                              const = const,
                                                              ws = ihw.null.ws,
                                                              n.group.max = n.group.max,
                                                              seed = seed)
    selected.2D.ihw.null <- res.2D.ihw.null$selected
    tm <- res.2D.ihw.null$tm0
    ta <- res.2D.ihw.null$ta0
    
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

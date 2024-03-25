#==================
# Simulation Step
# 1 dimension
# To evaluate the combination of spatial information and covariate information
#==================
one_step_1D.CalTime <- function(h, grp,detect.m = "rad.h", seed = 0,
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
    
  }else{
    # Load the previous result and save time.
    load(save.file.path)
  }
  
  for(hh in h){
    print(paste("seed:",seed,"hh:",hh)) 
    
    #==== Set pre-chosen Tm, Ta
    # Tm.star
    Tm.star <- round(Tm.star,dig)
    Ta.star <- Inf
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
    print(paste("seed:",seed,"BH","hh:",hh))
    res.2D <- Spatial_Detect_exact_grp_BH_down.CalTime(Tm, Ta, Va,
                                                        VmVa.cov, ind,
                                                        q, max.rej,
                                                        Tm.star = Tm.star,
                                                        Ta.star = Ta.star,
                                                        const = const,
                                                        seed = seed,
                                                        mua = mua,
                                                        sig_ratio = sig_ratio,
                                                        EmpMethod = EmpMethod,
                                                        is.quick.stop = F)
    
    
  }
  
  return(res.2D)
}


### Functions for calculating time
Spatial_Detect_exact_grp_BH_down.CalTime <- function (Tm, Ta, Va, VmVa.cov, ind, q = 0.1, max.rej = NULL, 
                                                      pis = NULL, cutoff = NULL, Tm.star = Inf, Ta.star = Inf, 
                                                      const = 0, seed = 0, dig = 7, tau.tm = 1, tau.ta = 1, mua = NULL, 
                                                      is.quick.stop = F, stop.step = 10, is.fullset = F, is.fixRectShape = F, 
                                                      sig_ratio = NULL, EmpMethod = c("NPEB", "DeCor")) {
  if (is.null(pis)) {
    pis <- rep(1, length(Tm))
  }
  if (is.null(sig_ratio)) {
    sig_ratio <- rep(1, length(Tm))
  }
  EmpMethod <- match.arg(EmpMethod)
  if (is.fixRectShape) {
    if (EmpMethod == "NPEB") {
      avrg_c <- mean(sig_ratio)
    }
    else {
      avrg_c <- mean((sig_ratio - VmVa.cov)/(1 - VmVa.cov^2))
    }
  }
  ind.step <- 0
  m <- length(Tm)
  eta <- Ta
  if (EmpMethod == "NPEB") {
    if (is.fullset) {
      mm <- GLmix(x = eta[ind])
      normalized.prob <- mm$y/sum(mm$y)
      num.par <- 1
      for (ind.par in (ind[1] + 1):(ind[2] - 1)) {
        ind.cur <- ind + ind.par - 1
        ind.cur <- ind.cur[ind.cur <= m]
        mm.tmp <- GLmix(x = eta[ind.cur])
        mm$x <- c(mm$x, mm.tmp$x)
        normalized.prob <- c(normalized.prob, mm.tmp$y/sum(mm.tmp$y))
        num.par <- num.par + 1
      }
      normalized.prob <- normalized.prob/num.par
    }
    else {
      mm <- GLmix(x = eta[ind])
      normalized.prob <- mm$y/sum(mm$y)
    }
    if (!is.null(mua)) {
      mm.pre <- table(mua)
      mm$x <- as.numeric(names(mm.pre))
      mm$y <- mm.pre/sum(mm.pre)
      normalized.prob <- mm$y/sum(mm$y)
    }
    fdr.est <- function(tm, ta, NP.max, etype) {
      p <- length(Tm)
      L <- numeric(p)
      NP <- sum(Tm >= tm & Ta >= ta)
      NP <- ifelse(is.na(NP), 0, NP)
      if (NP == 0) {
        FDP <- 0
        FD <- 0
      }
      else if (NP < NP.max) {
        FDP <- NA
        FD <- NA
      }
      else {
        grp.val <- unique(cbind(Va, VmVa.cov))
        for (j in 1:nrow(grp.val)) {
          vVa <- grp.val[j, 1]
          vVmVa.cov <- grp.val[j, 2]
          row.ind <- which(Va == vVa & VmVa.cov == vVmVa.cov)
          L[row.ind] <- L.cal(tm, ta, mm, normalized.prob, 
                              vVa, vVmVa.cov)
        }
        if (etype == "FDR") {
          FD <- sum(pis * L) + const
          FDP <- FD/NP
        }
      }
      return(list(FDP = FDP, NP = NP, FD = FD))
    }
  }
  else if (EmpMethod == "DeCor") {
    Ta <- Ta - VmVa.cov * Tm
    fdr.est <- function(tm, ta, NP.max, etype) {
      p <- length(Tm)
      L <- numeric(p)
      NP <- sum(Tm >= tm & Ta >= ta)
      NP <- ifelse(is.na(NP), 0, NP)
      if (NP == 0) {
        FDP <- 0
        FD <- 0
      }
      else if (NP < NP.max) {
        FDP <- NA
        FD <- NA
      }
      else {
        L <- (1 - pnorm(tm)) * as.integer(Ta >= ta)
        if (etype == "FDR") {
          FD <- sum(pis * L) + const
          FDP <- FD/NP
        }
      }
      return(list(FDP = FDP, NP = NP, FD = FD))
    }
  }
  pv.tm <- 1 - pnorm(Tm)
  pv.ta <- 1 - pnorm(Ta)
  ## No prune
  Num.Step0 = m^2
  ## After step 1
  Tm.cutoff <- Tm[pv.tm <= 1 & pv.ta <= 1]
  Ta.cutoff <- Ta[pv.tm <= 1 & pv.ta <= 1]
  cutoff.rec <- cutoff.gen.rec(Tm.cutoff, Ta.cutoff, Inf, 
                               Inf)
  Num.Step1 = dim(cutoff.rec$cutoff)[1]
  
  if (is.null(cutoff)) {
    Tm.cutoff <- Tm[pv.tm <= tau.tm & pv.ta <= tau.ta]
    Ta.cutoff <- Ta[pv.tm <= tau.tm & pv.ta <= tau.ta]
    if (is.fixRectShape) {
      Tm.star <- Inf
      max.rej <- 0
    }
    cutoff.rec <- cutoff.gen.rec(Tm.cutoff, Ta.cutoff, Tm.star, 
                                 Ta.star)
    cutoff <- cutoff.rec$cutoff
    ind.Tm <- cutoff.rec$ind.Tm
    if (is.fixRectShape) {
      cutoff[, 2] <- cutoff[, 1] * avrg_c
    }
  }
  ## After Step 2
  Num.Step2 = dim(cutoff)[1]
  
  FDP <- NULL
  NP <- NULL
  tm.cand.set <- NULL
  ta.cand.set <- NULL
  NP.max <- max.rej
  index.Tm <- max.rej + 1
  m_prime <- length(ind.Tm) - 1
  
  ## Count Step 3
  Num.Step3 = 0
  Num.Step4 = 0 ## earling stop
  is.equal = T
  while ((index.Tm <= m_prime) & ((!is.quick.stop) | ind.step <= 
                                  stop.step)) {
    i.up <- ind.Tm[index.Tm + 1] - 1
    cur.ind.Tm <- ind.Tm[index.Tm]
    i.down <- cur.ind.Tm
    tm.down <- cutoff[i.down, 1]
    ta.down <- cutoff[i.down, 2]
    NP.down <- sum(Tm >= tm.down & Ta >= ta.down)
    i.down <- i.down + max(0, NP.max - NP.down)
    
   
    while (i.up >= i.down) {
      tm.down <- cutoff[i.down, 1]
      ta.down <- cutoff[i.down, 2]
      obj <- fdr.est(tm.down, ta.down, NP.max, etype = "FDR")
      if (is.na(obj$FDP)) {
        i.down <- i.down + max(0, NP.max - obj$NP)
      }
      else if (obj$FDP <= q) {
        ind.step <- 0
        NP <- c(NP, obj$NP)
        FDP <- c(FDP, obj$FDP)
        tm.cand.set <- c(tm.cand.set, tm.down)
        ta.cand.set <- c(ta.cand.set, ta.down)
        NP.max <- obj$NP
        tm0 <- tm.down
        ta0 <- ta.down
        i.down <- i.down + 1
        
        if(ind.step > stop.step){
          is.equal = F
        }
      }
      else {
        FD.down <- obj$FD
        NP.down <- obj$NP
        min.REJ <- ceiling(FD.down/q)
        i.down <- i.down + max(1, min.REJ - NP.down)
      }
      
      ## Search one point, plus one
      Num.Step3= Num.Step3+1
      if(ind.step <= stop.step){
        Num.Step4= Num.Step4+1
      }
      if (is.fixRectShape) {
        i.down <- i.up
      }
    }
    ind.step <- ind.step + 1
    index.Tm <- index.Tm + 1
  }
  if (is.null(FDP)) {
    selected <- NULL
  }
  else {
    NP.cand <- NP[NP == max(NP)]
    FDP.cand <- FDP[NP == max(NP)]
    tm.cand <- tm.cand.set[NP == max(NP)]
    ta.cand <- ta.cand.set[NP == max(NP)]
    tm0 <- tm.cand[which.min(FDP.cand)]
    ta0 <- ta.cand[which.min(FDP.cand)]
  }
  pos <- Tm >= tm0 & Ta >= ta0
  selected <- which(pos)
  final.obj <- fdr.est(tm0, ta0, NP.max, etype = "FDR")
  final.fdr <- final.obj$FDP
  return(c(Num.Step0=Num.Step0,Num.Step1=Num.Step1,Num.Step2=Num.Step2,
           Num.Step3=Num.Step3,Num.Step4=Num.Step4,is.equal=is.equal))
}





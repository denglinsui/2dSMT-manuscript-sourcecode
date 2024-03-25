#==================
# Simulation Step
# 2 dimension estimeting GMLE using the whole set
#==================

#==================
# Simulation Step
# 2 dimension
# Compare the convolutional density of GMLE estimated by the subset and the whole set
#==================
one_step_2D_Gm_tGm_comp <- function(h, detect.m = "rad.h", seed = 0,
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
  
  pluszero <- mu>0
  
  I_S <- Init_Setting_2D(mu_type = mu_type, 
                         Cov_type = Cov_type,
                         magnitude = magnitude,
                         #mu_gen_machine = "uc.unif",
                         point = point,
                         Dist.p = Dist.p)
  mu <- I_S$mu
  Sigma.eps.p <- I_S$Sigma.eps.p
  n <- 1
  
  X <- MASS::mvrnorm(n = n, mu = mu, 
                     Sigma = Sigma.eps.p)
  X <- matrix(X, nrow = n)
  
  sgm <- sqrt(diag(Sigma.eps.p))
  
  #=== With Spatial Info
  
  for(ii in 1:length(h)){
    hh <- h[ii]
    #===========================================
    #====== Fix Reigion
    #===========================================
    
    # Detect Neighbor
    hh.seq <- rep(hh,m)
    Neigh_Detect_res <- TwoDSMT::Neigh_Detect(hh = hh.seq,
                                              X = X, 
                                              Dist = Dist.p,
                                              Sigma.eps = Sigma.eps.p,
                                              detect.m = detect.m,
                                              mu = mu)
    Ta <- Neigh_Detect_res$Ta
    Va <- Neigh_Detect_res$Va
    VmVa.cov <- Neigh_Detect_res$VmVa.cov
    ind <- Neigh_Detect_res$ind
    mua <- Neigh_Detect_res$mua
    sig_ratio <- Neigh_Detect_res$sig_ratio
    
    eta = Ta
    ## True distribution
    mixCoeff.true.full <- rep(1/(m),(m))/(sum(rep(1/(m),(m)))+0.001)
    dist.true.full <- UnivarMixingDistribution( Dlist=lapply(mua,
                                                        function(x){Norm(mean=x)}), 
                                           mixCoeff=mixCoeff.true.full)
    ## True distribution for partial data
    m.sub <- length(ind)
    mua.sub <- mua[ind]
    mixCoeff.true.sub <- rep(1/(m.sub),(m.sub))/(sum( rep(1/(m.sub),(m.sub)))+0.001)
    dist.true.sub <- UnivarMixingDistribution( Dlist=lapply(mua.sub,
                                                        function(x){Norm(mean=x)}), 
                                           mixCoeff=mixCoeff.true.sub)
    
    ## Using fullset
    mm.full <- GLmix(x = eta)
    normalized.prob.full <- mm.full$y/sum(mm.full$y)
    normalized.prob.full <- normalized.prob.full/(sum(normalized.prob.full)+0.001)
    dist.full <- UnivarMixingDistribution( Dlist=lapply(mm.full$x,
                                                         function(x){Norm(mean=x)}), 
                                       mixCoeff=normalized.prob.full)
    
    #dist.full <- UnivarMixingDistribution( Norm(2,1),Norm(0,1), mixCoeff=c(0.6,0.4))
    ## Using partial set
    mm.sub <- GLmix(x = eta[ind])
    normalized.prob.sub <- mm.sub$y/sum(mm.sub$y)
    normalized.prob.sub <- normalized.prob.sub/(sum(normalized.prob.sub)+0.001)
    dist.sub <- UnivarMixingDistribution( Dlist=lapply(mm.sub$x,
                                                        function(x){Norm(mean=x)}), 
                                           mixCoeff=normalized.prob.sub)
    
    
    
  }
  HL.dist = c(
    HellingerDist(dist.true.full, dist.true.sub),
    HellingerDist(dist.true.full, dist.full),
    HellingerDist(dist.true.full, dist.sub),
    HellingerDist(dist.true.sub, dist.sub),
    m.sub)
  names(HL.dist) = c("EmG-EmtG","EmG-GMLEfull","EmpG-GMLEsub","EmtG-GMLEsub","mSub")
  HL.dist
  return(HL.dist)
}

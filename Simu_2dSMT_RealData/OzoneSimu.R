## Before run the simulation for Ozone data, please set the dictionary to 2dSMT 
#==== Run 2D Alg
library(IHW)
library(reshape2)
library(data.table)
library(dplyr)
library(ggplot2)
#library(ggmap)
#state <- map_data("state")
library(data.table)
library(backports)
library(TwoDSMT)
library(CompRandFld)
#-- For detect algorithms
library(qvalue)
#source('Simu_2dSMT_Numerical/Spatial_Detection_ma.R')
#source('Simu_2dSMT_Numerical/Tools/Hard_Group.R')
source("Simu_2dSMT_Numerical/Tools/All_q_est_functions.R")
load("Simu_2dSMT_RealData/Data/Ozone_Tmp.RData")

# Sequence for generating data.
generate_data <- function(beta, mu, year, sd.eps, Sigma.eps , seed=1){
  set.seed(seed)
  m <- nrow(point.for.cov)
  n <- length(year)
  Data <- matrix(mu,m,n,byrow=F)+matrix(beta,m,n,byrow=F)*matrix(year,m,n,byrow=T)
  Data <- Data + t(MASS::mvrnorm(n = n,mu=rep(0,m), Sigma = Sigma.eps))
}

# Get the estimators for the simulated data
linear_fit <- function(data,year){
  fit <- lm(data~year)
  mu.hat <- summary(fit)$coefficients[1,1]
  beta.hat <- summary(fit)$coefficients[2,1]
  sd.hat <- summary(fit)$sigma
  
  res.stat <- fit$residuals
  res.mat <- res.stat/sd.hat
  return(list(mu.hat=mu.hat,beta.hat=beta.hat,sd.hat=sd.hat,res.mat=res.mat))
}


#== Estimate fdr & power
fdp = function(selected,mu) sum(mu[selected] <= 0) / max(1, length(selected))
Pow = function(selected,mu){
  if(sum(mu>0)==0){
    return(0)
  }else{
    return(sum(mu[selected] > 0) / sum(mu > 0))
  }
}

#==== We run with 2 neighbor
Onestep <- function(seed,Corr.est.rd){
  print(paste0("seed:",seed,":start"))
  # Set dimension
  m <- nrow(point.for.cov)
  n <- length(year)
  
  year <- TT[,2]
  beta <- T2_Stat$beta_hat
  mu <- T2_Stat$mu_hat
  sd.eps <- Res_sd
  Sigma.eps <-Corr.est.rd* (Res_sd%*%t(Res_sd))
  
  Data <- generate_data(beta, mu, year, sd.eps, Sigma.eps , seed=seed)
  
  Coeff <- matrix(unlist(apply(Data,1,function(x){linear_fit(x,year=year)})),nrow=3+12)
  beta_hat <- Coeff[2,]
  Res_sd <- Coeff[3,]
  Res_mat <- Coeff[4:15,]
  
  est.ind <- 1:loc.num
  rep.ind <- 1:year.num
  rem <- c(10,0.5,0.5)
  names(rem) <- c("scale","sill","nugget")
  
  print(paste0("seed:",seed,":estcov"))
  if(T){
    fit <- FitComposite(Res_mat[rep.ind,est.ind], coordx=point.for.cov[est.ind,],
                        corrmodel=corrmodel, likelihood="Full",type="Standard", 
                        #corrmodel=corrmodel, likelihood='Marginal',type='Pairwise',#optimizer = "CG",  
                        fixed = list(mean=0), start = as.list(rem),
                        distance = "Geod",replicates = length(rep.ind))
  }
  param <- as.list(fit$param)
  
  #--- Estimating the covariance matrix
  cov.est <- Covmatrix(point.for.cov[,1], point.for.cov[,2], 
                       corrmodel=corrmodel,
                       distance = "Geod",
                       param=param)
  Corr.est <- cov.est$covmatrix 
  
  Corr.beta.est <- Corr.est*inv_tTT
  Sigma.beta.est <- Corr.beta.est* (Res_sd%*%t(Res_sd))
  
  fdp_res <- NULL
  pow_res <- NULL
  
  beta0.seq <- seq(-0.5,-0.1,by=0.1)
  
  for(beta0 in beta0.seq){
    # set mu
    print(paste0("seed:",seed,",beta0",beta0))
    minussbeta0 <-  beta<beta0
    
    # Now, T2 is for the simulated data.
    T2_Stat$beta_hat <- beta_hat
    T2_Stat$beta_sd <- sqrt(diag(Sigma.beta.est))
    T2_Stat$T2 <-  (beta_hat-beta0)/T2_Stat$beta_sd
    
    T2_Stat <- T2_Stat[,.(Latitude,Longitude,
                          beta_hat,beta_sd,
                          State.Code,
                          T2)]
    q <- 0.1
    const <- q
    tau.tm <- 1
    detect.m <- "top.k"
    Dist.prox <- -Corr.est #-Corr.est can play the role of dist
    sgm <- sqrt(diag(Sigma.beta.est))
    
    #==== Perform Algorithm
    dig <- 7
    Tm <- -T2_Stat$T2 # Take minus to make the direction of hypothesizes is consistent to our procedure
    p.value <- 1 - pnorm(Tm)
    
    #================================
    #==== Without spatial Info: One D
    
    #==== BH
    #================================
    result <- qvalue(p.value,pi0 = 1)
    BH.selected <- which(result$qvalues<=q)
    
    #==== Storey
    #================================
    #==== Calculate pis2
    lambda.test <- 0.5#quantile(Tm,0.9)#
    pis.hata.tmp<- min(1,
                       mean(Tm<lambda.test)/pnorm(lambda.test))
    pis.hata <- rep(pis.hata.tmp,m)
    
    Storey.res <- qvalue(p.value,pi0 = pis.hata.tmp)
    Storey.selected <- which(Storey.res$qvalues<=q)
    
    #===========================================
    #--- IHW: We use the location of IHW as covariate
    # (Without null proportion estimation)
    #===========================================
    print(paste0("seed:",seed," ihw"))
    nbins.ihw <- 5
    Lat.int <- seq(min(T2_Stat$Latitude),max(T2_Stat$Latitude)+0.1,length.out=4)
    Lon.int <- seq(min(T2_Stat$Longitude),max(T2_Stat$Longitude)+0.1,length.out=4)
    inter.ind <- function(x,y){max(which(x-y>=0))}
    Lon.ind <- sapply(T2_Stat$Longitude,function(x){inter.ind(x,y=Lon.int)})
    Lat.ind <- sapply(T2_Stat$Latitude,function(x){inter.ind(x,y=Lat.int)})
    loc.ind <- factor(paste(Lon.ind,Lat.ind))
    ihw.res <- ihw(1-pnorm(Tm),loc.ind, q)#,nbins = nbins.ihw)
    ihw.selected <- which(ihw.res@df$adj_pvalue<q)
    ihw.ws <- ihw.res@df$weight # Extract weights
    
    #===========================================
    #==== SABHA
    #==== We use the parameter in Ring and Li
    #===========================================
    print(paste0("seed:",seed," sabha"))
    tau = 0.5; eps = 0.1 # parameters for SABHA
    ADMM_params = c(10^2, 10^3, 2, 5000, 1e-3) # alpha_ADMM,beta,eta,max_iters,converge_thr
    qhat = Solve_q_block(p.value,tau,eps,loc.ind,ADMM_params)
    SABHA_method = function(pvals, qhat, alpha, tau){
      pvals[pvals>tau] = Inf
      khat=max(c(0,which(sort(qhat*pvals)<=alpha*(1:length(pvals))/length(pvals))))
      which(qhat*pvals<=alpha*khat/length(pvals))
    }
    
    sab.selected <- SABHA_method(p.value, qhat, q, tau)
    #=====================
    #=== Run 1D Alg
    #=====================
    print(paste0("seed:",seed," 1D"))
    
    res.1D <- OneD_Detect(Tm, q,  pis = NULL,
                          const = const,
                          tau.tm = tau.tm)
    Tm.star <- res.1D$tm.min
    max.rej <- res.1D$max.rej
    selected.1D <- res.1D$selected
    Ta.star <- Inf
    
    res.1D.pis2 <- OneD_Detect(Tm, q, pis = pis.hata,
                               const = const,
                               tau.tm = tau.tm)
    Tm.star.pis2 <- res.1D.pis2$tm.min
    max.rej.pis2 <- res.1D.pis2$max.rej
    selected.1D.pis2 <- res.1D.pis2$selected
    Ta.star.pis2 <- Inf
    
    res.1D.ihw <- OneD_Detect_w(Tm, q,
                                ws = ihw.ws,
                                const = const,
                                tau.tm = tau.tm)
    tm.star.ihw <- res.1D.ihw$tm.min
    max.rej.ihw <- res.1D.ihw$max.rej
    selected.1D.ihw <- res.1D.ihw$selected
    ta.star.ihw <- 0
    
    #==== SABHA based
    ws.sabha.fun <- function(x){1/x}
    res.1D.sabha <- OneD_Detect_w(Tm, q, pis = qhat,
                                  ws.fun = ws.sabha.fun,
                                  const = const,
                                  tau.tm = tau.tm)
    tm.star.sabha <- res.1D.sabha$tm.min
    max.rej.sabha <- res.1D.sabha$max.rej
    selected.1D.sabha <- res.1D.sabha$selected
    ta.star.sabha <- 0
    
    #=== Run 2D Alg
    # We only consider hh = 2 
    for(hh in 2){
      print(paste0("seed:",seed,":start 2D"))
      hh.seq <- rep(hh,m)
      # Detect Neighbor
      
      Neigh_Detect_res <- Neigh_Detect(hh = hh.seq,
                                       X = matrix(Tm,nrow = 1), 
                                       Dist = Dist.prox,
                                       Sigma.eps = Corr.est,
                                       detect.m = detect.m)
      Ta <- Neigh_Detect_res$Ta
      Va <- Neigh_Detect_res$Va
      VmVa.cov <- Neigh_Detect_res$VmVa.cov
      ind <- Neigh_Detect_res$ind
      mua <- Neigh_Detect_res$mua
      #--- 2D.BH
      res.2D <- TwoDSMT::Spatial_Detect_exact_grp_BH_down(Tm, Ta, Va,
                                                          VmVa.cov, ind,
                                                          q, max.rej,
                                                          Tm.star = Tm.star,
                                                          Ta.star = Ta.star,
                                                          const = const,
                                                          is.quick.stop=T)
      selected.2D <- res.2D$selected
      tm <- res.2D$tm0
      ta <- res.2D$ta0
      #--- 2D Storey
      res.2D.pis2 <- TwoDSMT::Spatial_Detect_exact_grp_BH_down(Tm, Ta, Va,
                                                               VmVa.cov, ind,
                                                               q,
                                                               max.rej = max.rej.pis2,
                                                               pis = pis.hata,
                                                               Tm.star = Tm.star.pis2,
                                                               Ta.star = Ta.star.pis2,
                                                               const = const,
                                                               is.quick.stop=T)
      
      selected.2D.pis2 <- res.2D.pis2$selected
      tm <- res.2D.pis2$tm0
      ta <- res.2D.pis2$ta0
      
      #--- 2D.IHW
      res.2D.ihw <-  TwoDSMT::Spatial_Detect_exact_BH_down_reTm_reTa(Tm, Ta, Va,
                                                                     VmVa.cov, ind,
                                                                     q,
                                                                     max.rej = max.rej.ihw,
                                                                     pis = rep(1,m),
                                                                     pws.tm.star = tm.star.ihw,
                                                                     pws.ta.star = ta.star.ihw,
                                                                     const = const,
                                                                     ws = ihw.ws,
                                                                     n.group.max = 5,
                                                                     is.quick.stop=T)
      selected.2D.ihw <- res.2D.ihw$selected
      tm <- res.2D.ihw$tm0
      ta <- res.2D.ihw$ta0
      
      res.2D.sabha <-  TwoDSMT::Spatial_Detect_exact_BH_down_reTm_reTa(Tm, Ta, Va,
                                                                       VmVa.cov, ind,
                                                                       q,
                                                                       max.rej = max.rej.sabha,
                                                                       pis = qhat,
                                                                       pws.tm.star = tm.star.sabha,
                                                                       pws.ta.star = ta.star.sabha,
                                                                       const = const,
                                                                       ws.fun = ws.sabha.fun,
                                                                       n.group.max = 5,
                                                                       is.quick.stop=T)
      selected.2D.sabha <- res.2D.sabha$selected
      tm <- res.2D.sabha$tm0
      ta <- res.2D.sabha$ta0
      
      fdp_res <- c(fdp_res, 
                   fdp(BH.selected,minussbeta0),fdp(Storey.selected,minussbeta0), fdp(ihw.selected,minussbeta0), fdp(sab.selected,minussbeta0),
                   fdp(selected.1D,minussbeta0),fdp(selected.1D.pis2,minussbeta0), fdp(selected.1D.ihw,minussbeta0), fdp(selected.1D.sabha,minussbeta0),
                   fdp(selected.2D,minussbeta0),fdp(selected.2D.pis2,minussbeta0), fdp(selected.2D.ihw,minussbeta0), fdp(selected.2D.sabha,minussbeta0)
      )
      pow_res <- c(pow_res, 
                   Pow(BH.selected,minussbeta0),Pow(Storey.selected,minussbeta0), Pow(ihw.selected,minussbeta0), Pow(sab.selected,minussbeta0),
                   Pow(selected.1D,minussbeta0),Pow(selected.1D.pis2,minussbeta0), Pow(selected.1D.ihw,minussbeta0), Pow(selected.1D.sabha,minussbeta0),
                   Pow(selected.2D,minussbeta0),Pow(selected.2D.pis2,minussbeta0), Pow(selected.2D.ihw,minussbeta0), Pow(selected.2D.sabha,minussbeta0))
      print(paste0("seed:",seed,":end"))
      print(fdp_res)
      print(pow_res)
    }
    #  save.image(paste0("Result/Ozone_h",hh,"_",beta0,".RData"))
  }
  
  names(fdp_res) <- paste(rep(beta0.seq,each = 12), 
                          rep(c("BH","ST","IHW","SABHA", 
                                "1D(BH)","1D(ST)","1D(IHW)","1D(SA)",
                                "2D(BH)","2D(ST)","2D(IHW)","2D(SA)"), times=5))
  
  return(list(fdp_res = fdp_res,
              pow_res = pow_res,
              seed = seed))
}

library(doMC)
cores <- 20
registerDoMC(cores = cores)

reptime <- 100
for(rd.cor.type in c("Exponential","Gaussian","Empirical")){
  if(rd.cor.type=="Exponential"){
    Corr.est.rd <- Corr.est
  }
  if(rd.cor.type=="Gaussian"){
    Corr.est.rd <- Corr.est.gau
  }
  if(rd.cor.type=="Empirical"){
    Corr.est.rd <- Corr.est.emp
  }
  m <- nrow(point.for.cov)
  n <- length(year)
  
  rr <- foreach(jj = 1:reptime,
                .combine = cbind) %dopar% Onestep(seed=jj,Corr.est.rd=Corr.est.rd)
  fdp_res <- NULL
  pow_res <- NULL
  seed <- numeric(reptime)
  for(i in 1:reptime){
    fdp_res <- rbind(fdp_res, rr[[3*i-2]])
    pow_res <- rbind(pow_res, rr[[3*i-1]])
    seed <- rr[[3*i]]
  }
  
  fdp_pow_print <- round(rbind(colMeans(fdp_res),colMeans(pow_res)),3)
  fdp_pow_sd_print <- round(rbind(apply(fdp_res,2,sd),apply(pow_res,2,sd)),3)
  
  rownames(fdp_pow_print) <- c("FDP","POWER")
  print(fdp_pow_print)
  print(fdp_pow_sd_print)
  
  sp_out <- function(xx){sapply(xx,function(x){sprintf("%.3f",x)})}
  
  fdp_print <- matrix(paste(sp_out(fdp_pow_print[1,]),"(",sp_out(fdp_pow_sd_print[1,]),")",sep=""),
                      nrow=5,byrow=T)
  
  pow_print <- matrix(paste(sp_out(fdp_pow_print[2,]),"(",sp_out(fdp_pow_sd_print[2,]),")",sep=""),
                      nrow=5,byrow=T)
  
  colnames(pow_print) <- c("BH","ST","IHW","SABHA", 
                           "1D(BH)","1D(ST)","1D(IHW)","1D(SA)",
                           "2D(BH)","2D(ST)","2D(IHW)","2D(SA)")
  rownames(pow_print) <- seq(-0.5,-0.1,by=0.1)
  
  colnames(fdp_print) <- colnames(pow_print)
  rownames(fdp_print) <- rownames(pow_print)
  
  apply(pow_print,1,function(x){paste(x[c(2:4,(10:12))],collapse ="&")})
  apply(fdp_print,1,function(x){paste(x[c(2:4,(10:12))],collapse ="&")})
  save.image(paste0("Simu_2dSMT_RealData/Result/Simu_Ozone_",rd.cor.type,".RData"))
  
}
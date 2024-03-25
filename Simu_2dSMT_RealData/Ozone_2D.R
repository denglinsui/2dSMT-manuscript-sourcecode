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
library(CompRandFld)
#-- For detect algorithms
library(qvalue)
# source('~/project/2DSpatial/Spatial_Detection_ma.R')
library(TwoDSMT)
source('~/project/2DSpatial/Tools/Hard_Group.R')
source("~/project/2DSpatial/Tools/All_q_est_functions.R")
load("~/project/2DSpatial_RealData/Data/Ozone_Tmp.RData")
#==== We run with 2 neighbor
beta0.seq <- seq(-0.5,-0.1,by=0.1)
for(beta0 in beta0.seq){
  T2_Stat <- merge(Beta_Stat,Sd_Stat) 
  #T2_Stat <- T2_Stat[,.(Latitude,Longitude,beta_hat,beta_sd,State.Code,T2 = (beta_hat-beta0)/beta_sd)]
  T2_Stat <- T2_Stat[,.(Latitude,Longitude,beta_hat,beta_sd,State.Code,T2 = (beta_hat-beta0)/beta_sd)]
  
  m <- nrow(T2_Stat)
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
  res.1D.pis2 <- OneD_Detect(Tm, q, pis = pis.hata,
                             const = const,
                             tau.tm = tau.tm)
  Tm.star.pis2 <- res.1D.pis2$tm.min
  max.rej.pis2 <- res.1D.pis2$max.rej
  Ta.star.pis2 <- Inf
  
  res.1D.ihw <- OneD_Detect_w(Tm, q,
                              ws = ihw.ws,
                              const = const,
                              tau.tm = tau.tm)
  tm.star.ihw <- res.1D.ihw$tm.min
  max.rej.ihw <- res.1D.ihw$max.rej
  ta.star.ihw <- 0
  
  #==== SABHA based
  ws.sabha.fun <- function(x){1/x}
  res.1D.sabha <- OneD_Detect_w(Tm, q, pis = qhat,
                                ws.fun = ws.sabha.fun,
                                const = const,
                                tau.tm = tau.tm)
  tm.star.sabha <- res.1D.sabha$tm.min
  max.rej.sabha <- res.1D.sabha$max.rej
  ta.star.sabha <- 0
  
  #=== Run 2D Alg
  # Detect Neighbor
  for(hh in 2){
    hh.seq <- rep(hh,m)
    
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
    
    #--- 2D Storey
    is.quick <- T
    res.2D.pis2 <- Spatial_Detect_exact_grp_BH_down(Tm, Ta, Va,
                                                    VmVa.cov, ind,
                                                    q,
                                                    max.rej = max.rej.pis2,
                                                    pis = pis.hata,
                                                    Tm.star = Tm.star.pis2,
                                                    Ta.star = Ta.star.pis2,
                                                    const = const,
                                                    is.quick.stop = is.quick)
    
    selected.2D.pis2 <- res.2D.pis2$selected
    tm <- res.2D.pis2$tm0
    ta <- res.2D.pis2$ta0
    
    #--- 2D.IHW
    res.2D.ihw <- Spatial_Detect_exact_BH_down_reTm_reTa(Tm, Ta, Va,
                                                         VmVa.cov, ind,
                                                         q,
                                                         max.rej = max.rej.ihw,
                                                         pis = rep(1,m),
                                                         pws.tm.star = tm.star.ihw,
                                                         pws.ta.star = ta.star.ihw,
                                                         const = const,
                                                         ws = ihw.ws,
                                                         n.group.max = 5,
                                                         is.quick.stop = is.quick)
    selected.2D.ihw <- res.2D.ihw$selected
    tm <- res.2D.ihw$tm0
    ta <- res.2D.ihw$ta0
    
    res.2D.sabha <- Spatial_Detect_exact_BH_down_reTm_reTa(Tm, Ta, Va,
                                                           VmVa.cov, ind,
                                                           q,
                                                           max.rej = max.rej.sabha,
                                                           pis = qhat,
                                                           pws.tm.star = tm.star.sabha,
                                                           pws.ta.star = ta.star.sabha,
                                                           const = const,
                                                           ws.fun = ws.sabha.fun,
                                                           n.group.max = 5,
                                                           is.quick.stop = is.quick)
                                                           
    selected.2D.sabha <- res.2D.sabha$selected
    tm <- res.2D.sabha$tm0
    ta <- res.2D.sabha$ta0
    
    # Discriminate Detect T1 and T2
    T2_Stat$Detect.ST <- "Null"
    T2_Stat$Detect.ST[Storey.selected] <- "ST"
    T2_Stat$Detect.ST[selected.2D.pis2] <- "2D(ST)"
    T2_Stat$Detect.ST[intersect(Storey.selected,selected.2D.pis2)] <- "ST & 2D(ST)"
    
    T2_Stat$Order.ST <- 0
    T2_Stat$Order.ST[Storey.selected] <- 2
    T2_Stat$Order.ST[selected.2D.pis2] <- 3
    T2_Stat$Order.ST[intersect(Storey.selected,selected.2D.pis2)] <- 1
    if(F){
      
    p.ST <- ggplot() +
      geom_map(
        data = state, map = state,
        aes(map_id = region),
        color = "black", fill = "lightgray", size = 0.1
      )+
      theme_void() +
      geom_point(aes(x=Longitude,y=Latitude,
                     color=Detect.ST,shape=Detect.ST,size=Detect.ST),
                 #size=1,
                 data=T2_Stat %>% arrange(Order.ST))+
      scale_colour_manual(values = c("Null" = "#999999",
                                     "ST" = "blue",
                                     "2D(ST)" = "red",
                                     "ST & 2D(ST)" = "dark green"))+
      scale_shape_manual(values = c("Null" = 3,
                                    "ST" = 16,
                                    "2D(ST)" = 17,
                                    "ST & 2D(ST)" = 2))+
      scale_size_manual(values = c("Null" = 1,
                                   "ST" = 2,
                                   "2D(ST)" =2,
                                   "ST & 2D(ST)" = 1))+
      theme(legend.position = "bottom",
            legend.title=element_text(size=12), 
            legend.text=element_text(size=12)) 
    #scale_colour_gradientn(colours=c("red","gray","blue"))+
    #
    print(p.ST)
    
    ggsave(paste0("Fig/Ozone(ST)_h",hh,"_",beta0,".eps"),width = 5,height=3.3)
    
    }
    T2_Stat$Detect.IHW <- "Null"
    T2_Stat$Detect.IHW[ihw.selected] <- "IHW"
    T2_Stat$Detect.IHW[selected.2D.ihw] <- "2D(IHW)"
    T2_Stat$Detect.IHW[intersect(ihw.selected,selected.2D.ihw)] <- "IHW & 2D(IHW)"
    
    T2_Stat$Order.IHW <- 0
    T2_Stat$Order.IHW[ihw.selected] <- 2
    T2_Stat$Order.IHW[selected.2D.ihw] <- 3
    T2_Stat$Order.IHW[intersect(ihw.selected,selected.2D.ihw)] <- 1
    if(F){
    p.IHW <- ggplot() +
      geom_map(
        data = state, map = state,
        aes(map_id = region),
        color = "black", fill = "lightgray", size = 0.1
      )+
      theme_void() +
      geom_point(aes(x=Longitude,y=Latitude,
                     color=Detect.IHW,shape=Detect.IHW,size=Detect.IHW),
                 #size=1,
                 data=T2_Stat %>% arrange(Order.IHW))+
      scale_colour_manual(values = c("Null" = "#999999",
                                     "IHW" = "blue",
                                     "2D(IHW)" = "red",
                                     "IHW & 2D(IHW)" = "dark green"))+
      scale_shape_manual(values = c("Null" = 3,
                                    "IHW" = 16,
                                    "2D(IHW)" = 17,
                                    "IHW & 2D(IHW)" = 2))+
      scale_size_manual(values = c("Null" = 1,
                                   "IHW" = 2,
                                   "2D(IHW)" = 2,
                                   "IHW & 2D(IHW)" = 1))+
      theme(legend.position = "bottom",
            legend.title=element_text(size=12), 
            legend.text=element_text(size=12))
    #scale_colour_gradientn(colours=c("red","gray","blue"))+
    #
    print(p.IHW)
    
    ggsave(paste0("Fig/Ozone(IHW)_h",hh,"_",beta0,".eps"),width = 5,height=3.3)
    }
    
    T2_Stat$Detect.SA <- "Null"
    T2_Stat$Detect.SA[sab.selected] <- "SA"
    T2_Stat$Detect.SA[selected.2D.sabha] <- "2D(SA)"
    T2_Stat$Detect.SA[intersect(sab.selected,selected.2D.sabha)] <- "SA & 2D(SA)"
    
    T2_Stat$Order.SA <- 0
    T2_Stat$Order.SA[sab.selected] <- 2
    T2_Stat$Order.SA[selected.2D.sabha] <- 3
    T2_Stat$Order.SA[intersect(sab.selected,selected.2D.sabha)] <- 1
    if(F){
      p.SA <- ggplot() +
        geom_map(
          data = state, map = state,
          aes(map_id = region),
          color = "black", fill = "lightgray", size = 0.1
        )+
        theme_void() +
        geom_point(aes(x=Longitude,y=Latitude,
                       color=Detect.SA,shape=Detect.SA,size=Detect.SA),
                   #size=1,
                   data=T2_Stat %>% arrange(Order.SA))+
        scale_colour_manual(values = c("Null" = "#999999",
                                       "SA" = "blue",
                                       "2D(SA)" = "red",
                                       "SA & 2D(SA)" = "dark green"))+
        scale_shape_manual(values = c("Null" = 3,
                                      "SA" = 16,
                                      "2D(SA)" = 17,
                                      "SA & 2D(SA)" = 2))+
        scale_size_manual(values = c("Null" = 1,
                                     "SA" = 2,
                                     "2D(SA)" = 2,
                                     "SA & 2D(SA)" = 1))+
        theme(legend.position = "bottom",
              legend.title=element_text(size=12), 
              legend.text=element_text(size=12))
      #scale_colour_gradientn(colours=c("red","gray","blue"))+
      #
      print(p.SA)
      
      ggsave(paste0("Fig/Ozone(SA)_h",hh,"_",beta0,".eps"),width = 5,height=3.3)
    }
    
    save.image(paste0("Fig/Ozone_h",hh,"_",beta0,".RData"))
  }
  
}

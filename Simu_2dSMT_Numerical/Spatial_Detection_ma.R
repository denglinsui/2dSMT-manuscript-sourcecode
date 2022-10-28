##' Adaptive Neighbor
##' 
Adapt_Neigh <- function(Dist, X, Hh){
  m <- dim(Dist)[1]
  # Hh: the upper bound
  
  nei.order <- matrix(apply(Dist,1,order)[2:(Hh+1),],nrow=Hh)
  
  order.X <- matrix(ncol = m, nrow = Hh)
  for(i in 1:m){
    order.X[,i] <-  X[nei.order[,i]]
  }
  
  signif <- matrix(apply(order.X, 2, 
                         function(x){
                           cumsum(x)/sqrt(1:Hh)
                         }),nrow=Hh) 
  opt.h <- apply(signif, 2, which.max)
  
  return(opt.h)
}


##' Neibor Detection
Neigh_Detect <- function(hh, X, Dist, Sigma.eps, detect.m = "top.k",mu=NULL){
  if(detect.m == "top.k"){
    #The neighbor is determined by the top hh nearest point
    S.neigh <- lapply(1:m,
                      function(i){
                        ind1 <- order(Dist[i,])[2:(hh[i]+1)]
                        ind1}) 
  }
  
  if(detect.m == "rad.h"){
    #The neighbor is chosen based on the cutoff distance
    S.neigh <- lapply(1:m,
                      function(i){
                        ind1 <- which(Dist[i,]<=hh[i])
                        ind <- ind1[ind1!=i]
                        ind})
  }
  
  ind = Neigh_Partition(S.neigh)
  #ind = 1:m
  
  Va <- sapply(1:m,
               function(i){
                 ind.Va = S.neigh[[i]]
                 sqrt(sum(Sigma.eps[ind.Va,ind.Va]))})
  
  Ta <- sapply(1:m,
               function(i){
                 ind.Ta = S.neigh[[i]]
                 sum(X[,ind.Ta]) /(sqrt(n))})
  Ta <- Ta/Va
  if(!is.null(mu)){
    mua <- sapply(1:m,
                 function(i){
                   ind.Ta = S.neigh[[i]]
                   sum(mu[ind.Ta]) /(sqrt(n))})
    mua <- mua/Va
  }else{
    mua <-NULL
  }
  
  VmVa.cov <- sapply(1:m,
                     function(i){
                       ind.Cov = S.neigh[[i]]
                       sum(Sigma.eps[i,ind.Cov])})
  VmVa.cov <- VmVa.cov/Va/sqrt(diag(Sigma.eps))
  Va <- Va/Va
  return(list(Ta = Ta,
              Va = Va,
              VmVa.cov = VmVa.cov,
              ind = ind,
              mua = mua))
}

#===============================================================
#' Neighbor Partition
#' To make the neighbors of different location do not intersect.

Neigh_Partition <- function(S.neigh){
  S.len <- lapply(S.neigh,length)
  #S.len.order <- order(unlist(S.len))
  S.len.order <- 1:length(unlist(S.len))
  #S.len.order <- sample(seq(1,length(unlist(S.len)),3))
  
  Union <- NULL
  ind <- NULL
  for(ind.order in S.len.order){
    S.neigh.cur <- S.neigh[[ind.order]]
    if(length(intersect(S.neigh.cur,Union))==0){
      ind <- c(ind,ind.order)
      Union <- union(Union,S.neigh.cur)
    }
  }
  
  return(ind)
}

#Cufoff for exact: not identical to BH like
cutoff.gen <- function(Tm, Ta, Tm.star = Inf, Ta.star = Inf){
  # No ties are allowed!
  index <- order(Tm, decreasing = T)
  cutoff <- matrix(c(Inf, Inf), ncol = 2)
  Ta.tmp <- c()
  for(i in 1:length(Tm)){
    index.Tm <- index[i]
    Ta.tmp <- c(Ta.tmp,Ta[index.Tm])
    if(Tm[index.Tm]<=Tm.star){
      Ta.choose <- sort(Ta.tmp[Ta.tmp<=min(Ta[index.Tm],Ta.star)])
      #index.Ta.pre <- index[1:i]
      #index.Ta <- index[which(Ta[index.Ta.pre]<=Ta[index.Tm])]
      
      #tmp <- matrix(c(rep(Tm[index.Tm], length(index.Ta)),
      #                sort(Ta[index.Ta])), 
      #              ncol = 2) 
      tmp <- matrix(c(rep(Tm[index.Tm], length(Ta.choose)),
                      Ta.choose), 
                    ncol = 2) 
      # OneStep Optimize
      #tmp <- tmp[tmp[,1]<=Tm.star & tmp[,2]<=Ta.star,]
      cutoff <- rbind(cutoff, tmp)
    }
  }
  
  return(cutoff)
}

cutoff.gen.rec <- function(Tm, Ta, Tm.star = Inf, Ta.star = Inf){
  #=== Ties are allowed
  #=== We generate the corresponding replicates
  index <- order(Tm, decreasing = T)
  cutoff <- matrix(c(Inf, Inf), ncol = 2)
  Ta.tmp <- c()
  index.Tm.change <- c(1)
  index.Tm.change.end <- 2
  for(i in 1:length(Tm)){
    index.Tm.change[i+1] <- index.Tm.change.end
    index.Tm <- index[i]
    Ta.tmp[i] <- Ta[index.Tm]
    if(Tm[index.Tm]<=Tm.star){
      Ta.choose <- sort(Ta.tmp[Ta.tmp<=min(Ta[index.Tm],Ta.star)],
                        decreasing = T)
      #index.Ta.pre <- index[1:i]
      #index.Ta <- index[which(Ta[index.Ta.pre]<=Ta[index.Tm])]
      
      #tmp <- matrix(c(rep(Tm[index.Tm], length(index.Ta)),
      #                sort(Ta[index.Ta])), 
      #              ncol = 2) 
      Ta.choose.len <- length(Ta.choose)
      tmp <- matrix(c(rep(Tm[index.Tm], Ta.choose.len),
                      Ta.choose), 
                    ncol = 2) 
      #record the position of Tm
      index.Tm.change.end <- index.Tm.change.end + Ta.choose.len
      # OneStep Optimize
      #tmp <- tmp[tmp[,1]<=Tm.star & tmp[,2]<=Ta.star,]
      cutoff <- rbind(cutoff, tmp)
    }
  }
  index.Tm.change[length(Tm)+2] <- index.Tm.change.end
  
  return(list(cutoff = cutoff,
              ind.Tm = index.Tm.change))
}


cutoff.pw.gen.rec <- function(pws.tm, pws.ta, pws.tm.star = 0, pws.ta.star = 0){
  #=== Ties are allowed
  #=== The cutoffs are generated adjusting for the features of weighted p-values
  #=== A large amount of weighted p-value is 1, the replicates are useless in this case.
  
  pws.tm.rev <- -pws.tm
  pws.ta.rev <- -pws.ta
  
  index <- order(pws.tm.rev, decreasing = T)
  cutoff <- matrix(c(0, 0), ncol = 2) # rejects none.
  ta.tmp <- c()
  index.tm.change <- c(1)
  index.tm.change.end <- 2
  for(i in 1:length(pws.tm.rev)){
    index.tm.change[i+1] <- index.tm.change.end
    index.tm <- index[i]
    ta.tmp[i] <- pws.ta.rev[index.tm]
    if(pws.tm.rev[index.tm]<=-pws.ta.star){
      ta.choose <- sort(ta.tmp[ta.tmp<=min(pws.ta.rev[index.tm],pws.ta.star)],
                        decreasing = T)
      
      #--- Remove replicates with 1 with respect to pws.ta
      if(sum(ta.choose==-1)!=0){
        ta.choose <- c(ta.choose[ta.choose>-1],-1)
      }
      
      #index.Ta.pre <- index[1:i]
      #index.Ta <- index[which(Ta[index.Ta.pre]<=Ta[index.Tm])]
      
      #tmp <- matrix(c(rep(Tm[index.Tm], length(index.Ta)),
      #                sort(Ta[index.Ta])), 
      #              ncol = 2) 
      ta.choose.len <- length(ta.choose)
      tmp <- matrix(c(rep(pws.tm.rev[index.tm], ta.choose.len),
                      ta.choose), 
                    ncol = 2) 
      #record the position of Tm
      index.tm.change.end <- index.tm.change.end + ta.choose.len
      # OneStep Optimize
      #tmp <- tmp[tmp[,1]<=Tm.star & tmp[,2]<=Ta.star,]
      cutoff <- rbind(cutoff, tmp)
    }
    
    #=--- Remove replicates with 1 with respect to pws.tm
    if(pws.tm.rev[index.tm]==-1) break
  }
  index.tm.change[length(index.tm.change)+1] <- index.tm.change.end
  
  return(list(cutoff = -cutoff,
              ind.tm = index.tm.change))
}


##'Spatial Detect Algorithm
##'Preliminary
L.cal <- function(tm,ta,
                  mm,normalized.prob,
                  vVa,vVmVa.cov){
  #=== Estimating L
  x1 <- -tm 
  x2 <- tm 
  y1 <- -ta- mm$x #* vVa
  y2 <- ta- mm$x #* vVa
  A1 <- pbivnorm::pbivnorm(x = x1, y = y1, rho = vVmVa.cov)
  A2 <- pbivnorm::pbivnorm(x = x2, y = y1, rho = vVmVa.cov)
  A3 <- pbivnorm::pbivnorm(x = x2, y = y2, rho = vVmVa.cov)
  A4 <- pbivnorm::pbivnorm(x = x1, y = y2, rho = vVmVa.cov)
  
  B1 <- pnorm(x1)
  B2 <- pnorm(x2)
  C1 <- pnorm(y1)
  C2 <- pnorm(y2)
  
  if(F){
    L <- sum(normalized.prob * (1 + A1 + B1 + C1 + A3 - A2 - B2 - C2 - A4))
  }
  #One side
  L <- sum(normalized.prob *(1-B2 -C2 + A3)) 
  
  return(L)
}


##'One Dimension Detection
OneD_Detect <- function(Tm,
                        q=0.1,
                        pis = NULL,
                        const = 0,
                        tau.tm = 1){
  if(is.null(pis)){
    pis <- rep(1,length(Tm))
  }
  tau.Tm <- qnorm(1-tau.tm)
  p <- length(Tm)
  Tm.ind <- rank(-Tm)
  
  FDP <- (sum(pis)*(1-pnorm(Tm))+const)/Tm.ind
  
  if(sum(FDP<=q)!=0){
    threshold <- max(tau.Tm,min(Tm[FDP<=q]))
    selected <- which(Tm>=threshold)
    tm.min <- threshold
    max.rej <- length(selected)
  }else{
    selected <- NULL
    tm.min <- Inf
    max.rej <- 0
  }
  return(list(selected = selected,
              tm.min = tm.min,
              max.rej = max.rej))
}

##'One Dimension Detection for Ta
OneD_Detect_Ta <- function(Ta,
                           q=0.1,
                           ind,
                           pis = NULL){
  if(is.null(pis)){
    pis <- rep(1,length(Ta))
  }
  
  eta <- Ta
  #=== Perform Non-parametric Emprical Bayesian
  mm <- REBayes::GLmix(x = eta[ind])
  normalized.prob <- mm$y / sum(mm$y)
  
  p <- length(Ta)
  Ta.ind <- rank(-Ta)
  
  FD <- pnorm(Ta)
  
  FD.individual <- sapply(Ta, function(ta){1 - sum(pnorm(ta- mm$x)*normalized.prob)})
  FDP <- sum(pis)*FD.individual/Ta.ind
  
  if(sum(FDP<=q)!=0){
    threshold <- min(Ta[FDP<=q])
    selected <- which(Ta>threshold)
    ta.min <- threshold
    max.rej <- length(selected)
  }else{
    selected <- NULL
    ta.min <- Inf
    max.rej <- 0
  }
  return(list(selected = selected,
              ta.min = ta.min,
              max.rej = max.rej))
}


#'Equivalent Expression for weighted BH procedure
#'With weight and pis assigned separately
OneD_Detect_w <- function(Tm,
                          q=0.1,
                          pis = NULL,
                          const = 0,
                          #ws.fun = function(pis){(1-pis)/(pis)},
                          ws.fun = NULL,
                          ws = NULL,
                          tol = 1e-9,
                          tau.tm = 1){
  #Initialization
  if(is.null(pis)){
    pis <- rep(1,length(Tm))
  }
  
  # Automatically calculate the weight
  # pis is given previously or is assigned by its default value:1.
  if(is.null(ws)&is.null(ws.fun)){
    ws <- rep(1,m)
  }else if(is.null(ws)&(!is.null(ws.fun))){
    nu <- 10e-5
    pis.new <- pis
    pis.new[which(pis.new<nu)] <- nu # stabilization
    pis.new[which(pis.new>1-nu)] <- 1-nu # stabilization
    ws <- ws.fun(pis.new)
  }
  
  nu <- 1e-5
  ws[ws<nu] <- nu # Make sure no NaN generate in the 2D process
  
  pv.tm <- 1-pnorm(Tm)
  pws<-pmin(rep(1,m),pv.tm/ws)
  for(tm.search in sort(pws,decreasing = T)){
    FD <- sum(pis * (1-pnorm(qnorm(1-pmin(rep(tau.tm,m),ws*tm.search)))))+const
    Discover <- sum(Tm>=qnorm(1-pmin(rep(tau.tm,m),ws*tm.search))-tol)
    #print(Discover)
    fdp.hat <- FD/max(Discover,1)
    if(fdp.hat<q) break
  }
  
  
  if(sum(fdp.hat<=q)!=0){
    selected <- which(Tm>=qnorm(1-pmin(rep(tau.tm,m),ws*tm.search))-tol)
    max.rej <- length(selected)
    tm.min <- tm.search
  }else{
    selected <- NULL
    tm.min <- 0
    max.rej <- 0
    #'Equivalent Expression for LAWS or SABHA
    
  }
  return(list(selected = selected,
              tm.min = tm.min,
              max.rej = max.rej))
}

##'Exact Version
Spatial_Detect_exact_grp <- function(Tm, Ta, Va, VmVa.cov, ind,
                                     q=0.1, max.rej = NULL,
                                     pis = NULL, cutoff = NULL,
                                     Tm.star = Inf, Ta.star = Inf){
  if(is.null(pis)){
    pis <- rep(1,length(Tm))
  }
  
  #Initialization
  m <- length(Tm)
  eta <- Ta
  #=== Perform Non-parametric Emprical Bayesian
  mm <- REBayes::GLmix(x = eta[ind])
  normalized.prob <- mm$y / sum(mm$y)
  
  #=== Estimate false discovery
  fdr.est <- function (tm, ta, NP.max, etype){
    p <- length(Tm)
    L <- numeric(p)
    NP <- sum(Tm >= tm & Ta >= ta)
    NP <- ifelse(is.na(NP), 0, NP)
    
    if (NP == 0){
      FDP <- 0
    }else if(NP<NP.max){
      FDP <- NA
    }else{
      #group_by (Ta,VmVa.cov)
      grp.val <- unique(cbind(Va,VmVa.cov))
      for(j in 1:nrow(grp.val)){
        vVa <- grp.val[j,1]
        vVmVa.cov <- grp.val[j,2]
        row.ind <- which(Va==vVa & VmVa.cov==vVmVa.cov)
        L[row.ind] <- L.cal(tm,ta,
                            mm,normalized.prob,
                            vVa,vVmVa.cov)
      }
      
      if (etype == 'FDR') {
        # When pi0.est, this is simply summation
        FD <- sum(pis * L)
        FDP <- FD / NP
        #print(paste(FDP,p*(1-B2)/NP))
      }
    } 
    return(list(FDP = FDP, 
                NP = NP))	
  }
  
  #=== Initialize the searching grid
  if(is.null(cutoff)){
    cutoff <- cutoff.gen(Tm,Ta,Tm.star,Ta.star)
  }
  
  cutoff.rej <- apply(cutoff,1,function(t){sum(Tm >= t[1] & Ta >= t[2])})
  
  #== Reordering
  #Add sort also according to cutoff
  cutoff.rej.order <- order(-cutoff.rej,cutoff[,1],cutoff[,2],decreasing = T)
  cutoff.rej <- cutoff.rej[cutoff.rej.order]
  cutoff <- cutoff[cutoff.rej.order,]
  
  #Index Def
  cutoff.rej.ind <- sapply(0:m,
                           function(i){min(which(cutoff.rej>=i))})
  #cutoff.rej.ind <- c(cutoff.rej.ind,length(cutoff.rej))
  
  FDP <- NULL
  NP <- NULL
  tmta <- NULL
  NP.max <- 0
  conse <- 0
  cut.nrow <- nrow(cutoff)
  if(is.null(max.rej)){
    i <- 1
  }else{
    # We detect the rejection in one-dimension and begin with a pre-chosen starting point.
    i <- cutoff.rej.ind[max.rej+1]
  }
  while(i <= cut.nrow){
    tm <- cutoff[i,1]
    ta <- cutoff[i,2]
    
    #In this case, NP.max doesn't work.
    obj <- fdr.est(tm, ta, NP.max, etype="FDR")
    
    # We use the order that rejection increases
    # The first item ensures the rejection is meaningful
    if(is.na(obj$FDP)){
      #print(obj$NP,NP.max)
    }
    if(obj$FDP<=q){
      #Reset conse
      conse <- 0
      NP <- c(NP, obj$NP)
      FDP <- c(FDP, obj$FDP)
      tmta <- c(tmta, paste(tm, ta))
      NP.max <- obj$NP
      
      tm0 <- tm
      ta0 <- ta
      
      #print(paste(i,obj$NP,obj$FDP))
      # We don't consider the cases that 
      #   the rejection is lower than current maximum rejection num.
      if(obj$NP<m){
        i <- cutoff.rej.ind[NP.max+2]
      }else{
        break
      }
    }else{
      # We believe the FDR increases with NP increasing
      # Thus, once at a contour FDP doesn't below target level, we stop.
      if(i==(cutoff.rej.ind[obj$NP+2]-1)){
        conse <- conse+1
        print(paste("conse", conse, i,obj$NP,obj$FDP))
      } 
    }
    if(conse>1){
      i.seq <- (cutoff.rej.ind[NP.max+1]):(cutoff.rej.ind[NP.max+2]-1)
      
      FDP.est <- sapply(i.seq,
                        function(i){
                          res <- fdr.est(cutoff[i,1], cutoff[i,2], NP.max-1, etype="FDR")
                          return(res$FDP)
                        }
      )
      ind <- i.seq[which.min(FDP.est)]
      tm0 <- cutoff[ind,1]
      ta0 <- cutoff[ind,2]
      break
    }
    i <- i+1
  }
  
  pos <- Tm >= tm0 & Ta >= ta0
  
  selected <- which(pos)
  
  return(list(selected = selected,
              tm0 = tm0,
              ta0 = ta0))
}




##'Exact Version
Spatial_Detect_exact_grp_norder <- function(Tm, Ta, Va, VmVa.cov, ind,
                                            q=0.1, max.rej = NULL,
                                            pis = NULL, cutoff = NULL,
                                            Tm.star = Inf, Ta.star = Inf){
  if(is.null(pis)){
    pis <- rep(1,length(Tm))
  }
  
  #Initialization
  m <- length(Tm)
  eta <- Ta
  #=== Perform Non-parametric Emprical Bayesian
  mm <- REBayes::GLmix(x = eta[ind])
  normalized.prob <- mm$y / sum(mm$y)
  
  #=== Estimate false discovery
  fdr.est <- function (tm, ta, NP.max, etype){
    p <- length(Tm)
    L <- numeric(p)
    NP <- sum(Tm >= tm & Ta >= ta)
    NP <- ifelse(is.na(NP), 0, NP)
    
    if (NP == 0){
      FDP <- 0
    }else if(NP<NP.max){
      FDP <- NA
    }else{
      #group_by (Ta,VmVa.cov)
      grp.val <- unique(cbind(Va,VmVa.cov))
      for(j in 1:nrow(grp.val)){
        vVa <- grp.val[j,1]
        vVmVa.cov <- grp.val[j,2]
        row.ind <- which(Va==vVa & VmVa.cov==vVmVa.cov)
        L[row.ind] <- L.cal(tm,ta,
                            mm,normalized.prob,
                            vVa,vVmVa.cov)
      }
      
      if (etype == 'FDR') {
        # When pi0.est, this is simply summation
        FD <- sum(pis * L)
        FDP <- FD / NP
        #print(paste(FDP,p*(1-B2)/NP))
      }
    } 
    return(list(FDP = FDP, 
                NP = NP))	
  }
  
  #=== Initialize the searching grid
  if(is.null(cutoff)){
    cutoff <- cutoff.gen(Tm,Ta,Tm.star,Ta.star)
  }
  
  FDP <- NULL
  NP <- NULL
  tmta <- NULL
  NP.max <- 0
  conse <- 0
  i <- 0
  if(!is.null(max.rej)){
    NP.max <- max.rej
  } else{
    NP.max <- 0
  }
  cut.nrow <- nrow(cutoff)
  while(i <= cut.nrow){
    tm <- cutoff[i,1]
    ta <- cutoff[i,2]
    
    #In this case, NP.max doesn't work.
    obj <- fdr.est(tm, ta, NP.max, etype="FDR")
    
    # We use the order that rejection increases
    # The first item ensures the rejection is meaningful
    if(is.na(obj$FDP)){
      print(obj$NP)
      print(NP.max)
      print(paste0(obj$NP,NP.max))
    } else if(obj$FDP<=q){
      #Reset conse
      conse <- 0
      NP <- c(NP, obj$NP)
      FDP <- c(FDP, obj$FDP)
      tmta <- c(tmta, paste(tm, ta))
      NP.max <- obj$NP
      
      tm0 <- tm
      ta0 <- ta
      
    }
    i <- i + 1
  }
  pos <- Tm >= tm0 & Ta >= ta0
  
  selected <- which(pos)
  
  return(list(selected = selected))
}


##'Exact Version
##'Use One Dimension BH to accelarate
Spatial_Detect_exact_grp_BH <- function(Tm, Ta, Va, VmVa.cov, ind,
                                        q=0.1, max.rej = NULL,
                                        pis = NULL, cutoff = NULL,
                                        Tm.star = Inf, Ta.star = Inf,
                                        seed = 0,
                                        dig=7){
  if(is.null(pis)){
    pis <- rep(1,length(Tm))
  }
  
  #Initialization
  m <- length(Tm)
  eta <- Ta
  #=== Perform Non-parametric Emprical Bayesian
  mm <- REBayes::GLmix(x = eta[ind])
  normalized.prob <- mm$y / sum(mm$y)
  
  #=== Estimate false discovery
  fdr.est <- function (tm, ta, NP.max, etype){
    p <- length(Tm)
    L <- numeric(p)
    NP <- sum(Tm >= tm & Ta >= ta)
    NP <- ifelse(is.na(NP), 0, NP)
    
    if (NP == 0){
      FDP <- 0
      FD <- 0
    }else if(NP<NP.max){
      FDP <- NA
      FD <- NA
    }else{
      #group_by (Ta,VmVa.cov)
      grp.val <- unique(cbind(Va,VmVa.cov))
      for(j in 1:nrow(grp.val)){
        vVa <- grp.val[j,1]
        vVmVa.cov <- grp.val[j,2]
        row.ind <- which(Va==vVa & VmVa.cov==vVmVa.cov)
        L[row.ind] <- L.cal(tm,ta,
                            mm,normalized.prob,
                            vVa,vVmVa.cov)
      }
      
      if (etype == 'FDR') {
        # When pi0.est, this is simply summation
        FD <- sum(pis * L)
        FDP <- FD / NP
        #print(paste(FD,FDP))
      }
    } 
    return(list(FDP = FDP, 
                NP = NP,
                FD = FD))	
  }
  
  #=== Initialize the searching grid
  if(is.null(cutoff)){
    cutoff.rec <- cutoff.gen.rec(Tm,Ta,Tm.star,Ta.star)
    cutoff <- cutoff.rec$cutoff
    ind.Tm <- cutoff.rec$ind.Tm
  }
  
  FDP <- NULL
  NP <- NULL
  #tmta <- NULL
  tm.cand.set <- NULL
  ta.cand.set <- NULL
  NP.max <- max.rej
  
  index.Tm <- max.rej+1
  
  while(index.Tm <= m+1){
    #print(index.Tm)
    i.up <- ind.Tm[index.Tm+1]-1
    cur.ind.Tm <- ind.Tm[index.Tm]
    i.down <- cur.ind.Tm
    
    # We only consider the rej > current max rej
    tm.down <- cutoff[i.down,1]
    ta.down <- cutoff[i.down,2]
    
    NP.down <- sum(Tm >= tm.down & Ta >= ta.down)
    i.down <- i.down + max(0, NP.max - NP.down)
    
    while(i.up >= i.down){
      print(paste(i.up,i.down))
      #=== Search Up
      tm <- cutoff[i.up,1]
      ta <- cutoff[i.up,2]
      
      #In this case, NP.max doesn't work.
      obj <- fdr.est(tm, ta, NP.max, etype="FDR")
      
      # We use the order that rejection increases
      # The first item ensures the rejection is meaningful
      if(is.na(obj$FDP)){
        # Impossible happen
        # obj$FDP == NA <=> obj$NP is smaller than maximum rejection
        # There is no need to search with current Tm
        break
      } else if(obj$FDP<=q){
        # Save to find minimum FDP
        NP <- c(NP, obj$NP)
        FDP <- c(FDP, obj$FDP)
        tm.cand.set <- c(tm.cand.set, tm)
        ta.cand.set <- c(ta.cand.set, ta)
        #tmta <- c(tmta, paste(tm, ta))
        NP.max <- obj$NP
        tm0 <- tm
        ta0 <- ta
        
        break
      } else{
        # obj$FDP>q
        # We find a way to accelerate the calculation
        i.up <- i.up - 1
      }
      #=== Search Down
      # If not stand, we change i.down to accelerate
      tm.down <- cutoff[i.down,1]
      ta.down <- cutoff[i.down,2]
      obj <- fdr.est(tm.down, ta.down, NP.max, etype="FDR")
      
      if(is.na(obj$FDP)){
        # In case of ties
        i.down <- i.down + max(0, NP.max - obj$NP)
      } else if(obj$FDP<=q){
        # Save to find minimum FDP
        NP <- c(NP, obj$NP)
        FDP <- c(FDP, obj$FDP)
        tm.cand.set <- c(tm.cand.set, tm.down)
        ta.cand.set <- c(ta.cand.set, ta.down)
        #tmta <- c(tmta, paste(tm.down, ta.down))
        NP.max <- obj$NP
        tm0 <- tm.down
        ta0 <- ta.down
        
        i.down <- i.down +1
      } else{
        # obj$FDP>q
        # Actually, we could alse use this strategy on the x-axis to accelarate
        # We find a way to accelerate the calculation
        FD.down <- obj$FD
        NP.down <- obj$NP
        min.REJ <- ceiling(FD.down/q)
        # Actually, min.REJ>NP.down
        # Current rejecting number is not enough
        #i.down <- i.down + max(0, min.REJ-NP.down)
        i.down <- i.down + min.REJ - NP.down
      }
    }
    index.Tm <- index.Tm + 1
  }
  
  if(is.null(FDP)){
    selected <- NULL
  }else{
    NP.cand <- NP[NP==max(NP)]
    FDP.cand <- FDP[NP==max(NP)]
    #tmta.cand <- tmta[NP==max(NP)]
    tm.cand <- tm.cand.set[NP==max(NP)]
    ta.cand <- ta.cand.set[NP==max(NP)]
    tm0 <- tm.cand[which.min(FDP.cand)]
    ta0 <- ta.cand[which.min(FDP.cand)]
    #tmta <- unique(tmta.cand[which.min(FDP.cand)])
    #tm0 <- as.numeric(unlist(strsplit(tmta,split=' '))[1])
    #ta0 <- as.numeric(unlist(strsplit(tmta,split=' '))[2])
  }
  print(paste0("seed",seed,"start"))
  pos <- Tm >= tm0 & Ta >= ta0
  print(paste0("seed",seed,"end"))
  
  selected <- which(pos)
  
  return(list(selected = selected,
              tm0 = tm0,
              ta0 = ta0))
}


##'Exact Version
##'Use One Dimension BH to accelarate
Spatial_Detect_exact_grp_BH_down <- function(Tm, Ta, Va, VmVa.cov, ind,
                                             q=0.1, max.rej = NULL,
                                             pis = NULL, cutoff = NULL,
                                             Tm.star = Inf, Ta.star = Inf,
                                             const = 0,
                                             seed = 0,
                                             dig=7,
                                             tau.tm = 1,
                                             tau.ta = 1,
                                             mua = NULL){
  if(is.null(pis)){
    pis <- rep(1,length(Tm))
  }
  
  #Initialization
  m <- length(Tm)
  eta <- Ta
  #=== Perform Non-parametric Emprical Bayesian
  mm <- REBayes::GLmix(x = eta[ind])
  normalized.prob <- mm$y / sum(mm$y)
  if(!is.null(mua)){
    mm.pre <- table(mua)
    mm$x <- as.numeric(names(mm.pre))
    mm$y <- mm.pre/sum(mm.pre)
    normalized.prob <- mm$y / sum(mm$y)
  }
  #mm$x <- rep(0,length(mm$x))
  #=== Estimate false discovery
  fdr.est <- function (tm, ta, NP.max, etype){
    p <- length(Tm)
    L <- numeric(p)
    NP <- sum(Tm >= tm & Ta >= ta)
    NP <- ifelse(is.na(NP), 0, NP)
    
    if (NP == 0){
      FDP <- 0
      FD <- 0
    }else if(NP<NP.max){
      FDP <- NA
      FD <- NA
    }else{
      #group_by (Ta,VmVa.cov)
      grp.val <- unique(cbind(Va,VmVa.cov))
      for(j in 1:nrow(grp.val)){
        vVa <- grp.val[j,1]
        vVmVa.cov <- grp.val[j,2]
        row.ind <- which(Va==vVa & VmVa.cov==vVmVa.cov)
        L[row.ind] <- L.cal(tm,ta,
                            mm,normalized.prob,
                            vVa,vVmVa.cov)
      }
      
      if (etype == 'FDR') {
        # When pi0.est, this is simply summation
        FD <- sum(pis * L) + const
        FDP <- FD / NP
        #print(paste(FD,NP,FDP))
      }
    } 
    return(list(FDP = FDP, 
                NP = NP,
                FD = FD))	
  }
  
  pv.tm <- 1-pnorm(Tm)
  pv.ta <- 1-pnorm(Ta)
  
  #=== Initialize the searching grid
  if(is.null(cutoff)){
    Tm.cutoff <- Tm[pv.tm<=tau.tm & pv.ta<=tau.ta] # Add tau-censoring
    Ta.cutoff <- Ta[pv.tm<=tau.ta & pv.ta<=tau.ta] 
    cutoff.rec <- cutoff.gen.rec(Tm.cutoff,Ta.cutoff,Tm.star,Ta.star)
    cutoff <- cutoff.rec$cutoff
    ind.Tm <- cutoff.rec$ind.Tm
  }
  
  FDP <- NULL
  NP <- NULL
  #tmta <- NULL
  tm.cand.set <- NULL
  ta.cand.set <- NULL
  NP.max <- max.rej
  
  index.Tm <- max.rej+1
  
  while(index.Tm <= m+1){
    #print(index.Tm)
    i.up <- ind.Tm[index.Tm+1]-1
    cur.ind.Tm <- ind.Tm[index.Tm]
    i.down <- cur.ind.Tm
    
    # We only consider the rej > current max rej
    tm.down <- cutoff[i.down,1]
    ta.down <- cutoff[i.down,2]
    
    NP.down <- sum(Tm >= tm.down & Ta >= ta.down)
    i.down <- i.down + max(0, NP.max - NP.down)
    
    while(i.up >= i.down){
      #print(paste(i.up,i.down))
      
      #=== Search Down
      # If not stand, we change i.down to accelerate
      tm.down <- cutoff[i.down,1]
      ta.down <- cutoff[i.down,2]
      obj <- fdr.est(tm.down, ta.down, NP.max, etype="FDR")
      
      if(is.na(obj$FDP)){
        # In case of ties
        i.down <- i.down + max(0, NP.max - obj$NP)
      } else if(obj$FDP<=q){
        # Save to find minimum FDP
        NP <- c(NP, obj$NP)
        FDP <- c(FDP, obj$FDP)
        tm.cand.set <- c(tm.cand.set, tm.down)
        ta.cand.set <- c(ta.cand.set, ta.down)
        #print(paste(tm.down,ta.down,obj$FD,obj$NP,obj$FDP))
        #tmta <- c(tmta, paste(tm.down, ta.down))
        NP.max <- obj$NP
        tm0 <- tm.down
        ta0 <- ta.down
        
        i.down <- i.down +1
      } else{
        # obj$FDP>q
        # Actually, we could alse use this strategy on the x-axis to accelarate
        # We find a way to accelerate the calculation
        #print(paste(tm.down,obj$FD,obj$NP,obj$FDP))
        FD.down <- obj$FD
        NP.down <- obj$NP
        min.REJ <- ceiling(FD.down/q)
        # Actually, min.REJ>NP.down
        # Current rejecting number is not enough
        i.down <- i.down + max(1, min.REJ-NP.down)
        #i.down <- i.down + min.REJ - NP.down
      }
    }
    index.Tm <- index.Tm + 1
  }
  
  if(is.null(FDP)){
    selected <- NULL
  }else{
    NP.cand <- NP[NP==max(NP)]
    FDP.cand <- FDP[NP==max(NP)]
    #tmta.cand <- tmta[NP==max(NP)]
    tm.cand <- tm.cand.set[NP==max(NP)]
    ta.cand <- ta.cand.set[NP==max(NP)]
    tm0 <- tm.cand[which.min(FDP.cand)]
    ta0 <- ta.cand[which.min(FDP.cand)]
    #tmta <- unique(tmta.cand[which.min(FDP.cand)])
    #tm0 <- as.numeric(unlist(strsplit(tmta,split=' '))[1])
    #ta0 <- as.numeric(unlist(strsplit(tmta,split=' '))[2])
  }
  #print(paste0("seed",seed,"start"))
  pos <- Tm >= tm0 & Ta >= ta0
  #print(paste0("seed",seed,"end"))
  
  selected <- which(pos)
  
  final.obj <- fdr.est(tm0, ta0, NP.max, etype="FDR")
  final.fdr <- final.obj$FDP
  return(list(selected = selected,
              tm0 = tm0,
              ta0 = ta0,
              final.fdr = final.fdr))
}

##'Exact Version
##'Use One Dimension BH to accelarate
##'Reorder p-values with pi(s) (Actually, we do this in the form of Tm)
Spatial_Detect_exact_BH_down_reTm <- function(Tm, Ta, Va, VmVa.cov, ind,
                                              q=0.1, max.rej = NULL,
                                              pis = NULL, cutoff = NULL,
                                              Tm.star = Inf, Ta.star = Inf,
                                              const = 0,
                                              seed = 0,
                                              dig=7){
  if(is.null(pis)){
    pis <- rep(1,length(Tm))
  }
  
  #Initialization
  m <- length(Tm)
  eta <- Ta
  #=== Perform Non-parametric Emprical Bayesian
  mm <- REBayes::GLmix(x = eta[ind])
  normalized.prob <- mm$y / sum(mm$y)
  #mm$x <- rep(0,length(mm$x))
  #=== Estimate false discovery
  fdr.est <- function (tm, ta, pws, NP.max, etype){
    p <- length(Tm)
    L <- numeric(p)
    ind.tmp <- ws*tm>=1
    NP <- sum((Tm>= qnorm(pmin(rep(1,m),pmax(rep(0,m),1-ws*tm)))) & Ta>= ta)
    #NP <- sum(Tm >= tm & Ta >= ta)
    NP <- ifelse(is.na(NP), 0, NP)
    
    if (NP == 0){
      FDP <- 0
      FD <- 0
    }else if(NP<NP.max){
      FDP <- NA
      FD <- NA
    }else{
      #group_by (Ta,VmVa.cov)
      #No longer to be grouped
      for(j in (1:p)[!ind.tmp]){
        L[j] <- L.cal(qnorm(1-ws[j]*tm),ta,
                      mm,normalized.prob,
                      Va[j],VmVa.cov[j])
      }
      for(j in (1:p)[ind.tmp]){
        L[j] <- 1-sum(normalized.prob * pnorm(ta,mean = mm$x*Va[j]))
      }
      
      if (etype == 'FDR') {
        # When pi0.est, this is simply summation
        FD <- sum(pis * L) + const
        FDP <- FD / NP
        #print(paste(FDP,p*(1-B2)/NP))
      }
    }
    return(list(FDP = FDP,
                NP = NP,
                FD = FD))
  }
  
  #=== Initialize the searching grid
  if(is.null(cutoff)){
    cutoff.rec <- cutoff.gen.rec(Tm,Ta,Tm.star,Ta.star)
    cutoff <- cutoff.rec$cutoff
    ind.Tm <- cutoff.rec$ind.Tm
  }
  
  #=== Initialize reordered p-values
  nu <- 10e-5
  pis.new <- pis
  pis.new[which(pis.new<nu)] <- nu # stabilization
  pis.new[which(pis.new>1-nu)] <- 1-nu # stabilization
  ws <- (1-pis.new)/pis.new
  pv.tm <- 1-pnorm(Tm)
  pws<-pv.tm/ws
  
  #=== Searching for the optimal threshold
  FDP <- NULL
  NP <- NULL
  #tmta <- NULL
  tm.cand.set <- NULL
  ta.cand.set <- NULL
  NP.max <- max.rej
  
  index.Tm <- max.rej+1
  
  while(index.Tm <= m+1){
    #print(index.Tm)
    i.up <- ind.Tm[index.Tm+1]-1
    cur.ind.Tm <- ind.Tm[index.Tm]
    i.down <- cur.ind.Tm
    
    # We only consider the rej > current max rej
    tm.down <- cutoff[i.down,1]
    ta.down <- cutoff[i.down,2]
    
    NP.down <- sum(Tm >= tm.down & Ta >= ta.down)
    i.down <- i.down + max(0, NP.max - NP.down)
    
    while(i.up >= i.down){
      #print(paste(i.up,i.down))
      #print(i.down)
      #=== Search Down
      # If not stand, we change i.down to accelerate
      tm.down <- cutoff[i.down,1]
      ta.down <- cutoff[i.down,2]
      obj <- fdr.est(tm.down, ta.down, pws, NP.max, etype="FDR")
      
      if(is.na(obj$FDP)){
        # In case of ties
        i.down <- i.down + max(0, NP.max - obj$NP)
      } else if(obj$FDP<=q){
        # Save to find minimum FDP
        NP <- c(NP, obj$NP)
        FDP <- c(FDP, obj$FDP)
        tm.cand.set <- c(tm.cand.set, tm.down)
        ta.cand.set <- c(ta.cand.set, ta.down)
        #tmta <- c(tmta, paste(tm.down, ta.down))
        NP.max <- obj$NP
        tm0 <- tm.down
        ta0 <- ta.down
        
        i.down <- i.down +1
      } else{
        # obj$FDP>q
        # Actually, we could alse use this strategy on the x-axis to accelarate
        # We find a way to accelerate the calculation
        FD.down <- obj$FD
        NP.down <- obj$NP
        min.REJ <- ceiling(FD.down/q)
        # Actually, min.REJ>NP.down
        # Current rejecting number is not enough
        #i.down <- i.down + max(0, min.REJ-NP.down)
        i.down <- i.down + min.REJ - NP.down
      }
    }
    index.Tm <- index.Tm + 1
  }
  
  if(is.null(FDP)){
    selected <- NULL
  }else{
    NP.cand <- NP[NP==max(NP)]
    FDP.cand <- FDP[NP==max(NP)]
    #tmta.cand <- tmta[NP==max(NP)]
    tm.cand <- tm.cand.set[NP==max(NP)]
    ta.cand <- ta.cand.set[NP==max(NP)]
    tm0 <- tm.cand[which.min(FDP.cand)]
    ta0 <- ta.cand[which.min(FDP.cand)]
    #tmta <- unique(tmta.cand[which.min(FDP.cand)])
    #tm0 <- as.numeric(unlist(strsplit(tmta,split=' '))[1])
    #ta0 <- as.numeric(unlist(strsplit(tmta,split=' '))[2])
  }
  #print(paste0("seed",seed,"start"))
  ind.tmp <- which(ws*tm0>1)
  pos <- (Tm>= qnorm(pmin(rep(0,m),pmax(rep(0,m),1-ws*tm0)))) & Ta>= tp <- list()
  p1 <- list()
  for(mag in c(0,1,2,3)){
    I_S <- Init_Setting_1D(mu_type = mu_type,
                           Cov_type = Cov_type,
                           magnitude = mag,
                           point = point)
    mu <- I_S$mu
    Sigma.eps.p <- I_S$Sigma.eps.p
    X <- MASS::mvrnorm(n = n, mu = mu,
                       Sigma = Sigma.eps.p)
    X <- matrix(X, nrow = n)
    est.ind <- 1:m
    corrmodel <- "exponential"
    corrmodel <- "matern"
    corrmodel <- "stable"
    null <- which(mu==0)
    #est.ind <- which(X<=max(X[null]))
    for(k in 1:10){
      Tm.pre <- apply(X,2,function(x){sum(x)/sqrt(n)})
      Tm <- Tm.pre#/sqrt(diag(Sigma.eps.est.full))
      qval <- p.adjust(2*(1-pnorm(abs(Tm))), method = "BY")
      sel.ind <- which((qval<=0.1))
      est.ind <- which(!(1:m %in% unique(as.vector(sapply(sel.ind,
                                                          function(x){
                                                            pmin(rep(m,21),
                                                                 pmax(rep(0,21),
                                                                      x+(-10:10)))})))))
      #est.ind <- (1:m)[-c(sel.ind)]
      length(est.ind)
      #est.ind <- NULL
      #diff(sel.ind)
      print(length(est.ind))
      used_in_est <- rep("Neigh",m)
      used_in_est[est.ind] <- "False"
      used_in_est[sel.ind] <- "Rej"
      p[[mag+1]] <-
        ggplot(data.frame(x=point[,1],y=X,
                          used_in_est=used_in_est,
                          is_null = (mu==0),
                          mu=mu))+
        geom_point(aes(x=x,y=X,color=used_in_est))+
        geom_line(aes(x=x,y=mu),alpha=0.8)+
        ggtitle(paste0("BY(0.1),magnitude=",mag))
      print(p[[mag+1]])
      point.for.cov <- cbind(point[est.ind,],rep(1,length(est.ind)))
      print(paste0(c("seed",seed,":Start Full....")))
      fit <- FitComposite(#X[1:(n),est.ind]-matrix(fit.glm$fitted[est.ind],nrow=1),
        X[1:(n),est.ind],
        coordx=point.for.cov,
        corrmodel=corrmodel, likelihood="Full",fixed = list(mean=0),
        type="Standard",
        #type="Tapering",taper="Wendland1",maxdist=10,
        start = list(scale=0.5,sill=0.5,nugget = 1),
        replicate = n)
      param <- as.list(c(fit$param,mean=0))
      cov.est <- Covmatrix(cbind(point,rep(1,length(point))),
                           corrmodel=corrmodel,
                           param=param)
      Sigma.eps.est.full <- cov.est$covmatrix
      print(Sigma.eps.est.full[1,1:5])
    }
    library(reshape2)
    data1 <- data.frame(x=point[,1],
                        Cor_True=Sigma.eps.p[1,],
                        Cor_Est=Sigma.eps.est.full[1,]) %>%
      filter(x<1) %>%
      melt(id=c("x"))
    p1[[mag+1]] <- ggplot(data=data1,aes(x=x,y=value,color=variable))+
      geom_point()+ggtitle(paste0("Estimated Covariance,magnitude=",mag))
    p1[[mag+1]]
  }
  p <- list()
  p1 <- list()
  for(mag in c(0,1,2,3)){
    I_S <- Init_Setting_1D(mu_type = mu_type,
                           Cov_type = Cov_type,
                           magnitude = mag,
                           point = point)
    mu <- I_S$mu
    Sigma.eps.p <- I_S$Sigma.eps.p
    X <- MASS::mvrnorm(n = n, mu = mu,
                       Sigma = Sigma.eps.p)
    X <- matrix(X, nrow = n)
    est.ind <- 1:m
    corrmodel <- "exponential"
    corrmodel <- "matern"
    corrmodel <- "stable"
    null <- which(mu==0)
    #est.ind <- which(X<=max(X[null]))
    for(k in 1:10){
      Tm.pre <- apply(X,2,function(x){sum(x)/sqrt(n)})
      Tm <- Tm.pre#/sqrt(diag(Sigma.eps.est.full))
      qval <- p.adjust(2*(1-pnorm(abs(Tm))), method = "BY")
      sel.ind <- which((qval<=0.1))
      est.ind <- which(!(1:m %in% unique(as.vector(sapply(sel.ind,
                                                          function(x){
                                                            pmin(rep(m,21),
                                                                 pmax(rep(0,21),
                                                                      x+(-10:10)))})))))
      #est.ind <- (1:m)[-c(sel.ind)]
      length(est.ind)
      #est.ind <- NULL
      #diff(sel.ind)
      print(length(est.ind))
      used_in_est <- rep("Neigh",m)
      used_in_est[est.ind] <- "False"
      used_in_est[sel.ind] <- "Rej"
      p[[mag+1]] <-
        ggplot(data.frame(x=point[,1],y=X,
                          used_in_est=used_in_est,
                          is_null = (mu==0),
                          mu=mu))+
        geom_point(aes(x=x,y=X,color=used_in_est))+
        geom_line(aes(x=x,y=mu),alpha=0.8)+
        ggtitle(paste0("BY(0.1),magnitude=",mag))
      print(p[[mag+1]])
      point.for.cov <- cbind(point[est.ind,],rep(1,length(est.ind)))
      print(paste0(c("seed",seed,":Start Full....")))
      fit <- FitComposite(#X[1:(n),est.ind]-matrix(fit.glm$fitted[est.ind],nrow=1),
        X[1:(n),est.ind],
        coordx=point.for.cov,
        corrmodel=corrmodel, likelihood="Full",fixed = list(mean=0),
        type="Standard",
        #type="Tapering",taper="Wendland1",maxdist=10,
        start = list(scale=0.5,sill=0.5,nugget = 1),
        replicate = n)
      param <- as.list(c(fit$param,mean=0))
      cov.est <- Covmatrix(cbind(point,rep(1,length(point))),
                           corrmodel=corrmodel,
                           param=param)
      Sigma.eps.est.full <- cov.est$covmatrix
      print(Sigma.eps.est.full[1,1:5])
    }
    library(reshape2)
    data1 <- data.frame(x=point[,1],
                        Cor_True=Sigma.eps.p[1,],
                        Cor_Est=Sigma.eps.est.full[1,]) %>%
      filter(x<1) %>%
      melt(id=c("x"))
    p1[[mag+1]] <- ggplot(data=data1,aes(x=x,y=value,color=variable))+
      geom_point()+ggtitle(paste0("Estimated Covariance,magnitude=",mag))
    p1[[mag+1]]
  }
  20
  #print(paste0("seed",seed,"end"))
  
  selected <- which(pos)
  
  return(list(selected = selected,
              tm0 = tm0,
              ta0 = ta0))
}



##'Exact Version
##'Use One Dimension BH to accelerate
##'Reorder p-values with pi(s) (Actually, we do this in the form of Tm and Ta)
##'Instead of calculating weight function via Assume weight and pis are given respectively.
Spatial_Detect_exact_BH_down_reTm_reTa <- function(Tm, Ta, Va, VmVa.cov, ind,
                                                   q=0.1, max.rej = NULL,
                                                   pis = NULL, # The prior probability of being null.
                                                   ws = NULL, # weights for reordering.
                                                   cutoff = NULL,
                                                   pws.tm.star = 0, pws.ta.star = 0,
                                                   const = 0,
                                                   seed = 0,
                                                   dig=7,
                                                   ws.fun = function(pis){(1-pis)/pis},
                                                   n.group.max = 10,
                                                   is.exact.group = F,
                                                   tol=1e-9,
                                                   tau.tm=1,
                                                   tau.ta=1,
                                                   adj.up = 0.9,adj.ud = 1.1,
                                                   mua = NULL){
  if(tau.tm<0 | tau.tm>1|tau.ta<0 | tau.ta>1){
    stop("Please check the input: tau must between 0 and 1.")
  }
  
  eps.in <- 10^{-dig}
  
  #==== Initialization
  m <- length(Tm)
  eta <- Ta
  
  #---- Initialize pis and ws
  if(is.null(pis)){
    pis <- rep(1,m)
  }
  
  if(is.null(ws)&is.null(ws.fun)){
    ws <- rep(1,m) # no weight
  }else if(is.null(ws)&(!is.null(ws.fun))){
    nu <- 10e-5
    pis.new <- pis
    pis.new[which(pis.new<nu)] <- nu # stabilization
    pis.new[which(pis.new>1-nu)] <- 1-nu # stabilization
    ws <- ws.fun(pis.new)
  }
  
  nu <- 1e-5
  ws[ws<nu] <- nu # Make sure no NaN generate in the process
  
  #=== Perform Non-parametric Emprical Bayesian
  mm <- REBayes::GLmix(x = eta[ind])
  normalized.prob <- mm$y / sum(mm$y)
  
  if(!is.null(mua)){
    mm.pre <- table(mua)
    mm$x <- as.numeric(names(mm.pre))
    mm$y <- mm.pre/sum(mm.pre)
    normalized.prob <- mm$y / sum(mm$y)
  }
  #mm$x <- rep(0,length(mm$x))
  #=== Estimate false discovery
  #=== tm,ta is the threshold with respect to the weighted p-values! Not from tm,ta.
  #=== group.value and group.ind is to accelarate the calculation
  fdr.est <- function (tm, ta, ws, NP.max =0, etype="FDR",
                       q = 1,# Target fdr value
                       Grp_Val_s = NULL,
                       Grp_Ind_s = NULL,
                       Grp_Val_s_up = NULL,
                       #ws.up_est = NULL,
                       is.exact.group = F,
                       tau.tm = 1,
                       tau.ta = 1,
                       ud_est.FD.val = 0,
                       up_est.FD.val = 0,
                       FD.ctl = NULL,
                       adj.up = 0.9,
                       adj.ud = 1.1,
                       need.exact = FALSE){
    
    p <- length(Tm)
    L <- numeric(p)
    L.up_est <- numeric(p)
    L.ud_est<- numeric(p)
    
    #--- Corresponding Rejection
    NP <- sum((Tm>= qnorm(1-pmin(rep(tau.tm,m),ws*tm))-tol) & 
                (Ta>=qnorm(1-pmin(rep(tau.ta,m),ws*ta))-tol))
    
    #--- The index is for accelarating the calculation
    #--- In this way, we could calculate some false discovery probability in one dimension
    ind.tmp.tm <- ws*tm>=tau.tm
    ind.tmp.ta <- ws*ta>=tau.ta
    ind.tmp.NaN <- (ws*tm<=1e-14)|(ws*ta<=1e-14) # Avoid NaN
    #NP <- sum(Tm >= tm & Ta >= ta)
    #NP <- ifelse(is.na(NP), 0, NP)
    
    if(is.null(FD.ctl)){
      FD.ctl <- NP*q
    }
    
    if (NP == 0){
      FDP <- 0
      FD <- 0
    }else if(NP<NP.max){
      FDP <- NA
      FD <- NA
    }else{
      ind.part0 <- which(ind.tmp.NaN) #This part will induce NaN 
      ind.partm <- which(ind.tmp.tm&ind.tmp.ta&(!ind.tmp.NaN))
      ind.parta <- which(!ind.tmp.tm&(ind.tmp.ta)&(!ind.tmp.NaN))
      ind.part3 <- which((ind.tmp.tm)&(!ind.tmp.ta)&(!ind.tmp.NaN))
      ind.part4 <- which((!ind.tmp.tm)&(!ind.tmp.ta)&(!ind.tmp.NaN))
      
      # Common part
      L[ind.part0] <- ws[ind.part0]*tm # A more conservative estimation
      L[ind.partm] <- 1
      #L[ind.parta] <- tau.tm
      L[ind.parta] <- (1-pnorm(qnorm(1-pmin(rep(tau.tm,length(ind.parta)),ws[ind.parta]*tm))))
      
      L.part0 <- sum(pis[ind.part0]*L[ind.part0])
      L.partm <- sum(pis[ind.partm]*L[ind.partm])
      L.parta <- sum(pis[ind.parta]*L[ind.parta])
      L.part0.ud_est <- L.part0.up_est <- L.part0#
      L.partm.ud_est <- L.partm.up_est <- L.partm#
      L.parta.ud_est <- L.parta.up_est <- L.parta#
      #=== With Group
      # If the group is exact, we only need the first part!
      if(!(is.null(Grp_Val_s))&!(is.null(Grp_Ind_s))){
        for(i in (1:nrow(Grp_Val_s))){
          ind.tmp <- Grp_Ind_s[[i]]
          val.ud_est <- Grp_Val_s[i,]; val.up_est <- Grp_Val_s_up[i,]#Vm,VmVa,ws
          ws.ud_est <- val.ud_est[3]; ws.up_est <-val.up_est[3]
          interct.ind3 <- intersect(ind.tmp,ind.part3)
          if(length(interct.ind3)!=0){
            ind.part3.tmp <- intersect(ind.tmp,ind.part3)
            L.ud_est[interct.ind3] <- 1-sum(normalized.prob * pnorm(qnorm(1-ws.ud_est*ta),
                                                                    mean = mm$x))   
            L.up_est[interct.ind3] <- 1-sum(normalized.prob * pnorm(qnorm(1-min(1,ws.up_est*ta)),
                                                                    mean = mm$x))   
          }
        }
        L.part3.ud_est <- sum(pis[ind.part3]*L.ud_est[ind.part3])
        L.part3.up_est <- sum(pis[ind.part3]*L.up_est[ind.part3])
        
        for(i in (1:nrow(Grp_Val_s))){
          ind.tmp <- Grp_Ind_s[[i]]
          val.ud_est <- Grp_Val_s[i,]; val.up_est <- Grp_Val_s_up[i,]#Vm,VmVa,ws
          Va.ud_est <- val.ud_est[1]; VmVa.cov.ud_est <- val.ud_est[2]; ws.ud_est <- val.ud_est[3]
          Va.up_est <- val.up_est[1]; VmVa.cov.up_est <- val.up_est[2]; ws.up_est <- val.up_est[3]
          #print(ws.tmp*tm) 
          #print(ws.tmp*ta) 
          #print(ws.tmp*ta) 
          interct.ind4 <- intersect(ind.tmp,ind.part4)
          if(length(interct.ind4)!=0){
            #print(paste0(ws.ud_est*tm,ws.ud_est*ta))
            if((ws.ud_est*tm<=1e-14)|(ws.ud_est*ta<=1e-14)){
              L.ud_est[interct.ind4] <- ws.ud_est*tm
            }else{
              L.ud_est[interct.ind4] <- L.cal(qnorm(1-ws.ud_est*tm),qnorm(1-ws.ud_est*ta),
                                              mm,normalized.prob,
                                              Va.ud_est,VmVa.cov.ud_est)
            }
            L.up_est[interct.ind4] <- L.cal(qnorm(1-min(1-1e-14,ws.up_est*tm)),
                                            qnorm(1-min(1-1e-14,ws.up_est*ta)),
                                            mm,normalized.prob,Va.up_est,VmVa.cov.up_est)
          }
        }
        L.part4.ud_est <- sum(pis[ind.part4]*L.ud_est[ind.part4])
        L.part4.up_est <- sum(pis[ind.part4]*L.up_est[ind.part4])
        
        #-- A "Conservative" Estimation for FDP
        FD.ud_est.org <- L.part0+L.partm+L.parta+L.part3.ud_est+L.part4.ud_est  + const
        FD.ud_est <- FD.ud_est.org + ud_est.FD.val 
        FD.up_est.org <- L.part0+L.partm+L.parta+L.part3.up_est+L.part4.up_est + const
        FD.up_est <- FD.up_est.org - up_est.FD.val
        FDP.ud_est <- FD.ud_est/NP; FDP.up_est <- FD.up_est/NP
        #print(c(FDP.ud_est,FDP.up_est))
      }else{
        FDP.ud_est = 0 # To make sure the calculation goes to the second part.
        FDP.up_est = 1
      }
      
      #=== Without Group
      # If the conservative estimation indicate this threshold is possible to be active
      # We further explore it
      # adj.up = 0.9;adj.ud = 1.1
      if(is.exact.group){
        FDP <- FDP.ud_est
        FD <- FD.ud_est
      }else if((FDP.up_est>q*adj.up & FDP.ud_est<=q*adj.ud) | (need.exact == TRUE)){
        #if(FDP.ud_est<=q*adj.ud){
        #
        #if(FDP.ud_est<=q & !is.exact.group){
        #if(T){ # For Test
        #-- Va=1
        for(j in (1:p)[ind.part3]){
          L[j] <- 1-sum(normalized.prob * pnorm(qnorm(1-ws[j]*ta),
                                                mean = mm$x))
        }
        L.part3 <- sum(pis[ind.part3]*L[ind.part3])
        #No longer to be grouped
        for(j in (1:p)[ind.part4]){
          L[j] <- L.cal(qnorm(1-ws[j]*tm),qnorm(1-ws[j]*ta),
                        mm,normalized.prob,
                        Va[j],VmVa.cov[j])
        }
        L.part4 <- sum(pis[ind.part4]*L[ind.part4])
        
        FD <- L.part0+L.partm+L.parta+L.part3+L.part4 + const
        FDP <- FD/NP
        
        ud_est.FD.val <- FD-FD.ud_est.org # Update the approximate diff.
        up_est.FD.val <- FD.up_est.org - FD # Update the approximate diff.
        #print(c(FD.ud_est.org/NP,FDP.ud_est,FDP,FDP.up_est,FD.up_est.org/NP))
        #print(c(FD.ud_est.org,FD.ud_est,FD,FD.up_est,FD.up_est.org))
      }else if(FDP.up_est<=q*adj.up){
        FDP <- FDP.up_est
        FD <- FD.up_est
      }else{
        FDP <- FDP.ud_est
        FD <- FD.ud_est
      }
      #if (etype == 'FDR') {
      # When pi0.est, this is simply summation
      #FD <- sum(pis * L) + const
      #FDP <- FD / NP
      #print(paste(FDP,p*(1-B2)/NP))
      #}
    }
    return(list(FDP = FDP,
                NP = NP,
                FD = FD,
                ud_est.FD.val = ud_est.FD.val,
                up_est.FD.val = up_est.FD.val))
  }
  
  #=== Initialize reordered p-values
  #--- weighted p value for Tm
  pv.tm <- 1-pnorm(Tm)
  #pws.tm <- pmin(rep(1,m),pv.tm/ws) #pv.tm/ws
  pws.tm <- pv.tm/ws
  #--- weighted p value for Ta
  pv.ta <- 1-pnorm(Ta)
  #pws.ta<-  pmin(rep(1,m),pv.ta/ws)#pv.ta/ws
  pws.ta<-  pv.ta/ws
  
  #=== Initialize the searching grid
  #--- The cutoff is based on the weighted p-value
  if(is.null(cutoff)){
    pws.tm.cutoff <- pws.tm[pv.tm<=tau.tm & pv.ta<=tau.ta] 
    pws.ta.cutoff <- pws.ta[pv.tm<=tau.tm & pv.ta<=tau.ta] 
    cutoff.rec <- cutoff.pw.gen.rec(pws.tm.cutoff,pws.ta.cutoff,
                                    pws.tm.star,pws.ta.star)
    cutoff <- cutoff.rec$cutoff
    ind.pws.tm <- cutoff.rec$ind.tm
  }
  
  #=== Initilize the group
  #print("init")
  res.hg <- Hard_Group(Org_Value = cbind(Va,VmVa.cov,ws),n.group.max = n.group.max)
  Grp_Val_s <- res.hg$Grp_Val_summarize
  Grp_Val_s_up <- res.hg$Grp_Val_up_summarize
  Grp_Ind_s <- res.hg$Grp_Ind_summarize
  #print("end")
  #ws.up_est <- res.hg$Grp_Value
  is.exact.group <- res.hg$is.exact.group
  #print(is.exact.group)
  #res.hg_up <- Hard_Group_up(Org_Value = cbind(Va,VmVa.cov,ws),n.group.max = n.group.max)
  #Grp_Val_s_up <- res.hg_up$Grp_Val_summarize
  #Grp_Ind_s_up <- res.hg_up$Grp_Ind_summarize
  #is.exact.group_up <- res.hg_up$is.exact.group
  
  #=== Searching for the optimal threshold
  FDP <- NULL
  NP <- NULL
  #tmta <- NULL
  tm.cand.set <- NULL
  ta.cand.set <- NULL
  NP.max <- max.rej
  
  index.pws.tm<- max.rej+1
  
  #-- Adjust for the weighted p-value
  m_prime <- length(ind.pws.tm)-1
  FD.record <- matrix(numeric(0),ncol=2)
  while(index.pws.tm <= m_prime){
    #print(index.pws.tm)
    i.up <- ind.pws.tm[index.pws.tm+1]-1
    cur.ind.pws.tm <- ind.pws.tm[index.pws.tm]
    i.down <- cur.ind.pws.tm
    
    tm.up <- cutoff[i.up,1]
    ta.up <- cutoff[i.up,2]
    
    # We only consider the rej > current max rej
    tm.down <- cutoff[i.down,1]
    ta.down <- cutoff[i.down,2]
    
    NP.down <- sum((Tm>= qnorm(1-pmin(rep(tau.tm,m),ws*tm.down))-tol) & 
                     (Ta>=qnorm(1-pmin(rep(tau.ta,m),ws*ta.down))-tol))
    i.down <- i.down + max(0, NP.max - NP.down)
    ud_est.FD.val <- 0 # Initialize the estimated error
    up_est.FD.val <- 0
    while(i.up >= i.down){
      
      #=== Search Down
      # If not stand, we change i.down to accelerate
      tm.down <- cutoff[i.down,1]
      ta.down <- cutoff[i.down,2]
      obj <- fdr.est(tm = tm.down,ta = ta.down,ws = ws,
                     NP.max = NP.max, etype="FDR",q = q, 
                     Grp_Val_s = Grp_Val_s, Grp_Ind_s = Grp_Ind_s,
                     Grp_Val_s_up = Grp_Val_s_up, 
                     is.exact.group = is.exact.group,
                     tau.tm = tau.tm,
                     tau.ta = tau.ta,
                     adj.up = adj.up,adj.ud = adj.ud,
                     ud_est.FD.val = ud_est.FD.val,
                     up_est.FD.val = up_est.FD.val)
      ud_est.FD.val <- obj$ud_est.FD.val
      up_est.FD.val <- obj$up_est.FD.val
      FD.record <- rbind(FD.record,
                         c(ud_est.FD.val, up_est.FD.val))
      #  print(paste(i.up,i.down,obj$FDP))
      if(is.na(obj$FDP)){
        # In case of ties
        i.down <- i.down + max(0, NP.max - obj$NP)
      } else if(obj$FDP<=q){
        # Save to find minimum FDP
        NP <- c(NP, obj$NP)
        FDP <- c(FDP, obj$FDP)
        tm.cand.set <- c(tm.cand.set, tm.down)
        ta.cand.set <- c(ta.cand.set, ta.down)
        #tmta <- c(tmta, paste(tm.down, ta.down))
        NP.max <- obj$NP
        tm0 <- tm.down
        ta0 <- ta.down
        
        i.down <- i.down +1
      } else{
        # obj$FDP>q
        # Actually, we could alse use this strategy on the x-axis to accelarate
        # We find a way to accelerate the calculation
        FD.down <- obj$FD
        NP.down <- obj$NP
        min.REJ <- ceiling(FD.down/q)
        # Actually, min.REJ>NP.down
        # Current rejecting number is not enough
        #i.down <- i.down + max(0, min.REJ-NP.down)
        i.down.tmp <- i.down + max(min.REJ - NP.down,1)
        i.down <- ifelse((ta.down!=1)&(ta.up==1)&(i.down.tmp>i.up),
                         i.up,i.down.tmp)
      }
    }
    index.pws.tm <- index.pws.tm + 1
  }
  
  if(sum(FD.record[,1]<0 |  FD.record[,2]<0)>0){
    warning_mess <- paste0("seed:",seed,
                           " The surrogate FDP(ud) and FDP(up) are not exact upper and under estimator of FDP with at least ",
                           "Misspecified Under:",sum(FD.record[,1]<0),
                           " and Misspecified Upper:",sum(FD.record[,2]<0)
    )
    warning(warning_mess)
  }
  
  if(is.null(FDP)){
    selected <- NULL
    tm0 <- pws.tm.star # Conservative, to make the program run temproally.
    # For most methods, guaranteed by the one dimensional approach, we don't reach this step;
    # For LAWS, things becomes different because of the numerical digits
    ta0 <- max(1/ws)+tol
  }else{
    NP.cand <- NP[NP==max(NP)]
    FDP.cand <- FDP[NP==max(NP)]
    #tmta.cand <- tmta[NP==max(NP)]
    tm.cand <- tm.cand.set[NP==max(NP)]
    ta.cand <- ta.cand.set[NP==max(NP)]
    tm0 <- tm.cand[which.min(FDP.cand)]
    ta0 <- ta.cand[which.min(FDP.cand)]
    #tmta <- unique(tmta.cand[which.min(FDP.cand)])
    #tm0 <- as.numeric(unlist(strsplit(tmta,split=' '))[1])
    #ta0 <- as.numeric(unlist(strsplit(tmta,split=' '))[2])
  }
  #print(paste0("seed",seed,"start"))
  ind.tmp <- which(ws*tm0>1)
  #indd <- 3;tm0 <- tm.cand[indd];tm0 <- tm.cand[indd]
  # Report Selected point
  pos <- (Tm>= qnorm(1-pmin(rep(tau.tm,m),ws*tm0))-tol) & 
    (Ta>=qnorm(1-pmin(rep(tau.ta,m),ws*ta0))-tol)
  #print(paste0("seed",seed,"end"))
  selected <- which(pos)
  
  final.obj <- fdr.est(tm = tm0,ta = ta0,ws = ws, etype="FDR",q = q, 
                       Grp_Val_s = Grp_Val_s, Grp_Ind_s = Grp_Ind_s,
                       Grp_Val_s_up = Grp_Val_s_up, need.exact=T)
  final.fdr <- final.obj$FDP
  # Report Selected
  #fdp(selected,mu);Pow(selected,mu)
  return(list(selected = selected,
              tm0 = tm0,
              ta0 = ta0,
              final.fdr=final.fdr))
}


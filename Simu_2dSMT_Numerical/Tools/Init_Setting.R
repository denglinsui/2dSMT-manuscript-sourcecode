#==========================================
#==== Generate mu & covariance matrix
#==== One dimensional case
#'@param mu_type: Sparse Medium Dense
#'@param Cov_type: Weak Medium Strong
#'@param magnitude: The magnitude of signal 
#'@return Return mu and Sigma.eps
#'
#'

Init_Setting_1D <- function(mu_type = "Medium",
                            Cov_type = "Medium",
                            mu_gen_machine = "uc.spline",
                            magnitude = 1,
                            point = point,
                            single_bump_prop = 1/20,
                            Dist.p = NULL){
  m <- nrow(point)
  #=== Generate Dist.p
  if(is.null(Dist.p)){
    Dist.p <- as.matrix(dist(point))
  }
  #===========================================================================
  #=== Generate mu
  #single_bump_prop <- 1/10
  if(mu_gen_machine=="uc.spline" | mu_gen_machine == "uc.unif"){
    ind.mu <- switch(mu_type, 
                     Sparse = 1,
                     Medium = 2,
                     Dense = 4)
    num_bump = 1:ind.mu
    mu.res <- mu.fix.gen.1D(point, mu_gen_machine, magnitude,
                            num_bump=num_bump[ind.mu],single_bump_prop=single_bump_prop)
    mu <- mu.res$mu
  }else if(mu_gen_machine=="mvnorm"){
    mu.mean <- switch(mu_type, 
                      Sparse = -2.5,
                      Medium = -2,
                      Dense = -1)
    Sigma.mu <- K.mu.d(Dist.p, k =1, rho.mu = 0.3, sig2.mu=3)
    #set.seed(10)
    mu.res <- mu.fix.gen.1D(point, mu_gen_machine,
                            magnitude = magnitude, mu.mean = mu.mean,
                            Sigma.mu = Sigma.mu)
    mu <- mu.res$mu
  }
  #===========================================================================
  
  #===========================================================================
  #=== Generate Covariance Matrix
  r <- c(0.1,0.3,0.5)
  k <- c(1,1,2)
  rho.eps <- c(0.1,0.1,0.3)
  
  # For testing
  r <- c(0.5,0.9,0.9)
  k <- c(1,1,2)
  rho.eps <- c(0.05,0.1,0.2)
  
  r <- c(0.5,0.8,0.6)
  k <- c(1,1,2)
  rho.eps <- c(0.05,0.1,0.2)
  
  ind.Cov <- switch(Cov_type, 
                    Weak = 1,
                    Medium = 2,
                    Strong = 3)
  if(is.null(Dist.p)){
    Sigma.eps.p <- Sigma.eps(m, r[ind.Cov], k[ind.Cov], rho.eps[ind.Cov],point)
  }else{
    Sigma.eps.p <- K.eps.d(Dist.p, r[ind.Cov], k[ind.Cov], rho.eps[ind.Cov])
  }
  #===========================================================================
  
  return(list(mu = mu,
              Sigma.eps.p = Sigma.eps.p))
}


Init_Setting_1D_spa_grp <- function(
                            Cov_type = "Medium",
                            magnitude = 1,
                            point = point,
                            Dist.p = NULL){
  m <- nrow(point)
  #=== Generate Dist.p
  if(is.null(Dist.p)){
    Dist.p <- as.matrix(dist(point))
  }
  #===========================================================================
  #=== Generate mu
  #single_bump_prop <- 1/10
    single_bump_prop.seq = c(0,0,1/10,1/5)
    num_bump = 4
    
    #=== Choose the boundary
    center.ind <- floor(seq(0,1,length.out=num_bump+2)[-c(1,num_bump+2)]*m)
    each_er.seq <- floor(m*single_bump_prop.seq/2)
    
    #=== Init point
    point.basis <- point %>%  as.data.frame() %>% mutate(mu =0)
    colnames(point.basis) <- c("x","mu")
    
    #=== Start Generating
    for(ind in 1:num_bump){
      cen.ind = center.ind[ind]
      each_er = each_er.seq[ind]
      ind.low <- cen.ind-each_er; ind.up <- cen.ind+each_er
      p.low <- point[ind.low,1];p.up <- point[ind.up,1]
      
      #=== Init Basis
      n.break <- 3
      n.order <- 4
      n.end <- n.order + n.break -2
      if(ind.low!= ind.up){
      sp.basis <- create.bspline.basis(c(p.low, p.up), dropind=c(1,n.end),
                                       breaks = seq(p.low, p.up, length.out = n.break), 
                                       norder = n.order)
      
      point.basis[point.basis$x>=p.low & point.basis$x<=p.up,] <- 
        point.basis %>%
        filter(x>=p.low& x<=p.up) %>% 
        mutate(mu = mu + 1 * rowSums(eval.basis(x,basisobj = sp.basis)))
      }
    }
    
    mu <- magnitude * point.basis$mu
    pis <- as.numeric(mu<=0)
    ## assign group
    grp.start <- c(1,floor((center.ind[1:(num_bump-1)]+ center.ind[2:(num_bump)])/2))
    grp.end <- c(floor((center.ind[1:(num_bump-1)]+ center.ind[2:(num_bump)])/2)+1,m)
    grp <- rep(1,m)
    for(ind in 1:num_bump){
      grp[grp.start[ind]:grp.end[ind]] <- ind
    }
    
  
  #===========================================================================
  
  #===========================================================================
  #=== Generate Covariance Matrix
  r <- c(0.1,0.3,0.5)
  k <- c(1,1,2)
  rho.eps <- c(0.1,0.1,0.3)
  
  # For testing
  r <- c(0.5,0.9,0.9)
  k <- c(1,1,2)
  rho.eps <- c(0.05,0.1,0.2)
  
  r <- c(0.5,0.8,0.6)
  k <- c(1,1,2)
  rho.eps <- c(0.05,0.1,0.2)
  
  ind.Cov <- switch(Cov_type, 
                    Weak = 1,
                    Medium = 2,
                    Strong = 3)
  if(is.null(Dist.p)){
    Sigma.eps.p <- Sigma.eps(m, r[ind.Cov], k[ind.Cov], rho.eps[ind.Cov],point)
  }else{
    Sigma.eps.p <- K.eps.d(Dist.p, r[ind.Cov], k[ind.Cov], rho.eps[ind.Cov])
  }
  #===========================================================================
  
  return(list(mu = mu,
              Sigma.eps.p = Sigma.eps.p,
              grp = grp))
}

Init_Setting_Cov <- function(Cov_type = "Medium",
                             Dist.p = NULL){
  #===========================================================================
  #=== Generate Covariance Matrix
  r <- c(0,0.3,0.3)
  k <- c(1,1,2)
  rho.eps <- c(0.1,0.1,0.3)
  ind.Cov <- switch(Cov_type, 
                    Weak = 1,
                    Medium = 2,
                    Strong = 3)
  
  if(is.null(Dist.p)){
    Sigma.eps.p <- Sigma.eps(m, r[ind.Cov], k[ind.Cov], rho.eps[ind.Cov],point)
  }else{
    Sigma.eps.p <- K.eps.d(Dist.p, r[ind.Cov], k[ind.Cov], rho.eps[ind.Cov])
  }
  #===========================================================================
  
  return(Sigma.eps.p)
}


#==========================================
#==== Generate mu & covariance matrix
#==== Two dimensional case
#'@param mu_type: Sparse Medium Dense
#'@param Cov_type: Weak Medium Strong
#'@param magnitude: The magnitude of signal 
#'@return Return mu and Sigma.eps
#'
#'

Init_Setting_2D <- function(mu_type = "Medium",
                            Cov_type = "Medium",
                            mu_gen_machine = NULL,
                            magnitude = 1,
                            point = point,
                            single_bump_prop = 1/20,
                            Dist.p = NULL){
  m <- nrow(point)
  
  #===========================================================================
  #=== Generate mu
  #single_bump_prop <- 1/10
  methods <- c("uc.square","uc.circle","uc.mix")
  ind.mu <- switch(mu_type, 
                   Sparse = 1,
                   Medium = 2,
                   Dense = 3)
  num_bump = 1:ind.mu
  mu.res <- mu.fix.gen.2D(point, method=methods[ind.mu], magnitude)
  mu <- mu.res$mu
  #===========================================================================
  
  #===========================================================================
  #=== Generate Covariance Matrix
  r <- c(0.1,0.3,0.5)
  k <- c(1,1,2)
  rho.eps <- c(0.1,0.1,0.3)
  
  # For testing
  r <- c(0.5,0.9,0.9)
  k <- c(1,1,2)
  rho.eps <- c(0.05,0.1,0.2)
  
  r <- c(0.5,0.8,0.6)
  k <- c(1,1,2)
  rho.eps <- c(0.05,0.1,0.2)
  
  ind.Cov <- switch(Cov_type, 
                    Weak = 1,
                    Medium = 2,
                    Strong = 3)
  if(is.null(Dist.p)){
    Sigma.eps.p <- Sigma.eps(m, r[ind.Cov], k[ind.Cov], rho.eps[ind.Cov],point)
  }else{
    Sigma.eps.p <- K.eps.d(Dist.p, r[ind.Cov], k[ind.Cov], rho.eps[ind.Cov])
  }
  #===========================================================================
  
  return(list(mu = mu,
              Sigma.eps.p = Sigma.eps.p))
}

Init_Setting_Cov <- function(Cov_type = "Medium",
                             Dist.p = NULL){
  #===========================================================================
  #=== Generate Covariance Matrix
  r <- c(0,0.3,0.3)
  k <- c(1,1,2)
  rho.eps <- c(0.1,0.1,0.3)
  ind.Cov <- switch(Cov_type, 
                    Weak = 1,
                    Medium = 2,
                    Strong = 3)
  
  if(is.null(Dist.p)){
    Sigma.eps.p <- Sigma.eps(m, r[ind.Cov], k[ind.Cov], rho.eps[ind.Cov],point)
  }else{
    Sigma.eps.p <- K.eps.d(Dist.p, r[ind.Cov], k[ind.Cov], rho.eps[ind.Cov])
  }
  #===========================================================================
  
  return(Sigma.eps.p)
}
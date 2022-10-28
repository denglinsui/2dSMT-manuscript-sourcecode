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
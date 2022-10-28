library(fda)
library(dplyr)

#== Estimate fdr & power
fdp = function(selected,mu) sum(mu[selected] <= 0) / max(1, length(selected))
Pow = function(selected,mu){
  if(sum(mu>0)==0){
    return(0)
  }else{
    return(sum(mu[selected] > 0) / sum(mu > 0))
  }
}
# Kernel Function
K.mu <- function(s, t, rho.mu){
  sig.mu^2*exp(-(norm(as.matrix(s - t), type = "F")/rho.mu)^k)
}

K.mu.d <- function(d, k, rho.mu,sig2.mu){
  sig2.mu*exp(-(d/rho.mu)^k)
}

K.eps <- function(s, t, r, k, rho.eps){
  (1-r)*all(s==t) + r*exp(-(norm(as.matrix(s - t), type = "F")/rho.eps)^k)
}

K.eps.d <- function(d, r, k, rho.eps){
  (1-r)*(d==0) + r*exp(-(d/rho.eps)^k)
}

#=== Init Sigma
Dist <- function(m){
  Dist.up <- sapply(1:m, function(i){
    c(sapply(1:i, function(j) {
      norm(as.matrix(point[i,] - point[j,]))
    }),rep(0,m-i))
  })
  Dist <- Dist.up + t(Dist.up) - diag(diag(Dist.up))
}

Sigma.mu <- function(m, rho.mu){
  Sigma.mu.up <- sapply(1:m, function(i){
    c(sapply(1:i, function(j) {
      K.mu(point[i,], point[j,], rho.mu)
    }),rep(0,m-i))
  })
  Sigma.mu <- Sigma.mu.up + t(Sigma.mu.up) - diag(diag(Sigma.mu.up))
}

Sigma.eps <- function(m, r, k, rho.eps, point){
  Sigma.eps.up <- sapply(1:m, function(i){
    c(sapply(1:i, function(j) {
      K.eps(point[i,], point[j,], r, k, rho.eps)
    }),rep(0,m-i))
  })
  Sigma.eps <- Sigma.eps.up + t(Sigma.eps.up) - diag(diag(Sigma.eps.up))
}



user_mvrnorm <- function (n = 1, mu, Sigma, tol = 1e-06, empirical = FALSE, EISPACK = FALSE) 
{
  p <- length(mu)
  if (!all(dim(Sigma) == c(p, p))) 
    stop("incompatible arguments")
  if (EISPACK) 
    stop("'EISPACK' is no longer supported by R", domain = NA)
  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L]))) 
    stop("'Sigma' is not positive definite")
  X <- matrix(rnorm(p * n), n)
  if (empirical) {
    X <- scale(X, TRUE, FALSE)
    X <- X %*% svd(X, nu = 0)$v
    X <- scale(X, FALSE, TRUE)
  }
  len.vec <- length(eS$values)
  eS.vec <- matrix(sign(eS$vectors[1,]),len.vec,len.vec,byrow = T) * eS$vectors
  X <- drop(mu) + eS.vec %*% diag(sqrt(pmax(ev, 0)), p) %*% 
    t(X)
  nm <- names(mu)
  if (is.null(nm) && !is.null(dn <- dimnames(Sigma))) 
    nm <- dn[[1L]]
  dimnames(X) <- list(nm, NULL)
  if (n == 1) 
    drop(X)
  else t(X)
}

## Init mu
if(FALSE){
mu.gen <- function(point, method, magnitude,mu.mean,Sigma.mu){
  m <- nrow(point)
  if(method == "uc.unif"){
    location <- (point[,1]-1/2)^2+(point[,2]-1/2)^2<(1/4)^2
    mu <- runif(nrow(point),1,magnitude)*location
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
    
    mu <- runif(1,1,magnitude) * point.basis$mu
    pis <- as.numeric(mu<=0)
    
  }
  
  if(method == "mvnorm"){
    mu <- MASS::mvrnorm(n = 1, mu = rep(mu.mean,m), Sigma = magnitude*Sigma.mu)
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
}
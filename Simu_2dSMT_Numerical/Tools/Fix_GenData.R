#=================================
## 1 dimension
## Init mu (Fix magnitude)
#=================================
mu.fix.gen.1D <- function(point, method, magnitude,mu.mean,Sigma.mu,
                          num_bump=2,single_bump_prop=1/6 #For uc.spline
                          ){
  m <- length(point)
  
  if(method == "uc.unif"){
    if(num_bump*single_bump_prop>1) warning("The condition for Splines is not satisfied.")
    #=== Choose the boundary
    center.ind <- floor(seq(0,1,length.out=num_bump+2)[-c(1,num_bump+2)]*m)
    each_er <- floor(m*single_bump_prop/2)
    
    #=== Init point
    point.basis <- point %>%  as.data.frame() %>% mutate(mu =0)
    colnames(point.basis) <- c("x","mu")
    
    #=== Start Generating
    for(cen.ind in center.ind){
      ind.low <- cen.ind-each_er; ind.up <- cen.ind+each_er
      p.low <- point[ind.low,1];p.up <- point[ind.up,1]
      
      #=== Init Basis
      n.break <- 3
      n.order <- 4
      n.end <- n.order + n.break -2
      
      sp.basis <- create.bspline.basis(c(p.low, p.up), dropind=c(1,n.end),
                                       breaks = seq(p.low, p.up, length.out = n.break), 
                                       norder = n.order)
      
      point.basis[point.basis$x>=p.low & point.basis$x<=p.up,] <- 
        point.basis %>%
        filter(x>=p.low& x<=p.up) %>% 
        mutate(mu = mu + 1 * rowSums(eval.basis(x,basisobj = sp.basis)))
      
    }
    
    pis <- 1 * point.basis$mu
    pis[pis<=0.01] <- 0.01 #normalize
    mu <- magnitude * rbinom(n=length(pis),size=1,prob=pis)
  }
  
  
  if(method == "uc.spline"){
    if(num_bump*single_bump_prop>1) warning("The condition for Splines is not satisfied.")
    #=== Choose the boundary
    center.ind <- floor(seq(0,1,length.out=num_bump+2)[-c(1,num_bump+2)]*m)
    each_er <- floor(m*single_bump_prop/2)
    
    #=== Init point
    point.basis <- point %>%  as.data.frame() %>% mutate(mu =0)
    colnames(point.basis) <- c("x","mu")
    
    #=== Start Generating
    for(cen.ind in center.ind){
      ind.low <- cen.ind-each_er; ind.up <- cen.ind+each_er
      p.low <- point[ind.low,1];p.up <- point[ind.up,1]
      
      #=== Init Basis
      n.break <- 3
      n.order <- 4
      n.end <- n.order + n.break -2
      
      sp.basis <- create.bspline.basis(c(p.low, p.up), dropind=c(1,n.end),
                                         breaks = seq(p.low, p.up, length.out = n.break), 
                                         norder = n.order)
      
      point.basis[point.basis$x>=p.low & point.basis$x<=p.up,] <- 
        point.basis %>%
        filter(x>=p.low& x<=p.up) %>% 
        mutate(mu = mu + 1 * rowSums(eval.basis(x,basisobj = sp.basis)))
      
    }
    
    mu <- magnitude * point.basis$mu
    pis <- as.numeric(mu<=0)
    
  }
  
  if(method == "mvnorm"){
    mu <- magnitude * MASS::mvrnorm(n = 1, mu = rep(mu.mean,m), Sigma = Sigma.mu)
    pis <- as.numeric(mu<=0)
  }
  
  if(method == "mixture"){
    m <- nrow(point)
    pis <- numeric(m)
    location <- as.numeric(abs(point[,1]-5)<1) + as.numeric(abs(point[,1]-18)<2)
    pis[location] <- 0.9
    pis[-location] <- 0.1
    
    location <- rbinom(m,1,pis)
    mu <- location*magnitude
    pis <- 1-pis
  }
  return(list(mu=mu,
              pis=pis))
}

#=================================
## 2 dimension
## Init mu (Fix magnitude)
#=================================

mu.fix.gen.2D <- function(point, method, magnitude,mu.mean,Sigma.mu){
  m <- nrow(point)
  loc.mid <- (min(point[,1])+max(point[,1]))/2
  range.mid <- (max(point[,1])-min(point[,1]))
  if(method == "uc.circle"|method == "uc.mix"){
    location <- (point[,1]-loc.mid)^2+(point[,2]-loc.mid)^2<(1/4*range.mid)^2
    pis <- numeric(m)
    pis[location] <- 0.9
    pis[!location] <- 0.01
    mu <- magnitude * rbinom(n=length(pis),size=1,prob=pis)
    pis <- as.numeric(mu<=0)
    
  }
  
  if(method == "uc.square"|method == "uc.mix"){
    #=== Init Basis
    dis.low <-(3/4 +sqrt(2)/16-1/8) * range.mid
    dis.up <- (3/4+sqrt(2)/16+1/8) * range.mid
    n.break <- 3
    n.order <- 4
    n.end <- n.order + n.break -2
    sp.basis <- create.bspline.basis(c(dis.low,dis.up), dropind=c(1,n.end),
                                     breaks = seq(dis.low, dis.up, length.out = n.break), 
                                     norder = n.order)
    
    point.basis <- point %>%  as.data.frame() %>% mutate(mu.x = 0) %>% mutate(mu.y = 0)
    colnames(point.basis) <- c("x","y","mu.x","mu.y")
    
    #=== Basis product
    point.basis[point.basis$x>=dis.low & point.basis$x<=dis.up,] <- 
      point.basis %>%
      filter(x>=dis.low & x<=dis.up) %>% 
      mutate(mu.x = sqrt(n.end)/sqrt(n.end-2) * rowSums(eval.basis(x,basisobj = sp.basis)))
    
    point.basis[point.basis$y>=dis.low & point.basis$y<dis.up,] <- 
      point.basis %>%
      filter(y>=dis.low & y<=dis.up) %>% 
      mutate(mu.y = sqrt(n.end)/sqrt(n.end-2) * rowSums(eval.basis(y,basisobj = sp.basis)))
    #sqrt(n.end)/sqrt(n.end-2) is for normalization
    point.basis <- point.basis %>% mutate(mu = mu.x*mu.y) %>% mutate(mu = mu)
    if(method == "uc.mix"){
      mu <- mu + magnitude * point.basis$mu
    }else{
      mu <- magnitude * point.basis$mu
    }
    pis <- as.numeric(mu<=0)
  }
  
  if(method == "mvnorm"){
    mu <- magnitude*MASS::mvrnorm(n = 1, mu = rep(mu.mean,m), Sigma = Sigma.mu)
    pis <- as.numeric(mu<=0)
  }
  
  if(method == "mixture"){
    m <- nrow(point)
    pis <- numeric(m)
    location <- (point[,1]-loc.mid)^2+(point[,2]-loc.mid)^2<(1/4 * range.mid)^2
    pis[location] <- 0.9
    pis[-location] <- 0.1
    
    location <- rbinom(m,1,pis)
    mu <- location*magnitude
    pis <- 1-pis
  }
  return(list(mu=mu,
              pis=pis))
}
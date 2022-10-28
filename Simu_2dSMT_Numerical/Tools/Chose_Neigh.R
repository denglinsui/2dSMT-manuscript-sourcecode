#' We are going to test whether the fast searching algorithm is efficient.
#' 

rm(list=ls())
set.seed(10)
options(digits=5)
require(qvalue)
require(ggplot2)
require(np)
library("doMC")
setwd('~/project/Multiple-Testing-Replication/TDFDR_X/Spatial')
source('Spatial_Detection.R')
source('Tools/Basic.R')
source('Tools/laws_used.R')

#Initialize Parameter
m <- 900 # point size
m.row <- 30
n <- 1 # sample size
mu.mean <- -1
sig.mu <- 2
rho.mu <- 0.2
rho.eps <- 0.05
r <- 0.9
k <- 1
magnitude <- 3
parallel <- F

q <- 0.1
h.bd <- 3
data.type <- "mvnorm"
# Init point
point <- matrix(c(rep(seq(0,1,length.out=m.row),each=m.row),
                  rep(seq(0,1,length.out=m.row),times=m.row)), 
                m, 2, byrow = F) 

#Initial Dist and Sigma
Dist.p <- Dist(m)
Sigma.mu.p <- Sigma.mu(m)
Sigma.eps.p <- diag(rep(1,m))#Sigma.eps(m)

#=== Init Background
point <- matrix(c(rep(seq(0,1,length.out=m.row),each=m.row),
                  rep(seq(0,1,length.out=m.row),times=m.row)), 
                m, 2, byrow = F) 

#Initial Dist, Sigma and X
mu.res <-  mu.gen(point, "uc.unif", magnitude)
mu <- mu.res$mu
pis.hat3 <- mu.res$pis

pluszero <- mu>0
position <- pluszero

X <- MASS::mvrnorm(n = n, mu = mu, 
                   Sigma = Sigma.eps.p)
X <- as.numeric(X)

## For fixed point
par(mfrow=c(4,1))
par(mar=c(2,2,2,2))
plus.ind <- which(pluszero)
i.seq <- c(1,270,plus.ind[c(1,80)])
for(i in i.seq){
  p.dist <- Dist.p[i,]
  p.or <- order(p.dist,decreasing=F)

  magnitude <- cumsum((X[p.or]^2-1))[2:m]/sqrt(c(1:(m-1)))
  plot(magnitude,type="l")
  if(i %in% plus.ind){
    title(main=paste("non null(",round(point[i,1],3),
                     ",",round(point[i,2],3),")")
          )
  }else{
    title(main=paste("null(",round(point[i,1],3),
                     ",",round(point[i,2],3),")")
    )
  }
}

plot(X[p.or]^2-1,type="l")


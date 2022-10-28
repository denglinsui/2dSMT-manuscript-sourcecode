#' We are going to test whether the fast searching algorithm is efficient.
#' 

rm(list=ls())
set.seed(10)
options(digits=5)
require(qvalue)
require(ggplot2)
library("doMC")
library(plotly)
library(ggpubr)
setwd('~/project/Multiple-Testing-Replication/TDFDR_X/Spatial')
source('~/project/Multiple-Testing-Replication/TDFDR_X/Spatial/Spatial_Detection.R')
source('Spatial_Detection.R')
source('Tools/Basic.R')

#Initialize Parameter
m <- 900 # point size
m.row <- 30
n <- 1 # sample size
mu.mean <- -1
sig.mu <- 2
rho.mu <- 0.05
rho.eps <- 0.05
r <- 0.9
k <- 1
parallel <- F
magnitude <- 3

q <- 0.1
#K <-  function(s, t) {exp(- 16 * (s - t) ^ 2)}

# Init point
point <- matrix(c(rep(seq(0,1,length.out=m.row),each=m.row),
                  rep(seq(0,1,length.out=m.row),times=m.row)), 
                m, 2, byrow = F) 

#Initial Dist and Sigma
Dist <- Dist(m)
Sigma.mu <- Sigma.mu(m)
Sigma.eps <- Sigma.eps(m)


# Illustrate the false discovery and undiscovery
# Algorithm Implementation
mu <- mu.gen(point, "uc.spline", magnitude)
pluszero <- mu>0
X <- MASS::mvrnorm(n = n, mu = mu, 
                   Sigma = Sigma.eps)
X <- matrix(X, nrow = n)
  
  #==== Perform Algorithm
  T1 <- apply(X,2,function(x){sum(x)/sqrt(n)})
  
  #=== Without spatial Info: BH
  p.value <- 1 - pnorm(T1)
  result <- qvalue(p.value,pi0 = 1)
  selected.1 <- which(result$qvalues<=q)
  selected <- selected.1
  
  #==== Without spatial Info: One D
  res.1D <- OneD_Decect(T1, q)
  t1.min <- res.1D$t1.min 
  max.rej <- res.1D$max.rej
  
  #== Detect Neibor & Cal T2
  hh <- 5
  Neigh_Detect_res <- Neigh_Detect(hh = hh,
                                   X = X, 
                                   Dist = Dist.p,
                                   Sigma.eps = Sigma.eps.p,
                                   detect.m = "top.k")
  T2 <- Neigh_Detect_res$T2
  V2 <- Neigh_Detect_res$V2
  V1V2.cov <- Neigh_Detect_res$V1V2.cov
  ind <- Neigh_Detect_res$ind
  
  #== Round numerical values to ensures the correctness of taking maximum
  dig <- 5
  T1 <- round(T1,dig)
  T2 <- round(T2,dig)
  V2 <- round(V2,dig)
  V1V2.cov <- round(V1V2.cov,dig)
  #== Run Spatial Selection
  res.fast <- Spatial_Detect_exact_fast(T1, T2, V2, V1V2.cov, ind,
                                        q, t1.min, max.rej)
  selected.fast <- res.fast$selected
  res.exact <- Spatial_Detect_exact(T1, T2, V2, V1V2.cov, ind,
                                    q, max.rej)
  selected.exact <- res.exact$selected
  t1 <- res.exact$t10
  t2 <- res.exact$t20
  
  
  #========
  #== Plot
  #========
  pos.1 <- rep(F,m); pos.exact <- rep(F,m)
  pos.1[selected.1]<-T
  pos.exact[selected.exact]<-T

  selected[pos.1&pos.exact]<-"both"
  selected[!pos.1&pos.exact]<-"Spa"
  selected[pos.1&!pos.exact]<-"BH"
  selected[!pos.1&!pos.exact]<-"Null"
  data <- data.frame(x= point[,1], y= point[,2], X=X,
                     mu = mu, pluszero = pluszero, selected = selected)
  pic <- ggplot(data = data,mapping = aes(x = x, y = y, 
                                          shape = pluszero,
                                          color = selected))+
    geom_point()+
    scale_color_manual(breaks = c("both", "Null", "Spa","BH"),
                       values=c("blue", "gray", "black","green"))
  pic
  pic1 <- ggplot(data=data,aes(x=x,y=y,color=pluszero))+
    geom_point()
  pic1
  
  pic2 <- ggplot(data=data,aes(x=x,y=y,color=X))+
    geom_point()
  pic2
  
  fdp(selected.1,mu);Pow(selected.1,mu) 
  fdp(selected.exact,mu);Pow(selected.exact,mu) 
  fdp(selected.fast,mu);Pow(selected.fast,mu) 


#2D Rejection Region
  result <- qvalue(p.value,pi0 = 1)
  selected.1 <- which(result$qvalues<=q)
  t1.min <- min(T1[selected.1])
  t1 <- min(T1[selected.exact])
  t2 <- min(T2[selected.exact])
  data <- data.frame(T1 = T1, T2 = T2,
                     mu = mu, pluszero = pluszero)
  pic2d <- 
    ggplot(data = data,mapping = aes(x = T1, y = T2, color=pluszero))+
    geom_point(alpha=0.5)+
    geom_vline(mapping = aes(xintercept=t1),color="red")+
    geom_hline(mapping = aes(yintercept=t2),color="red")+
    geom_vline(mapping = aes(xintercept=t1.min),linetype = "dashed")
  pic2d


if(F){##Save all picture in one page
  obj <- ggarrange(pic, pic2d,
                   labels = c("A", "B"),
                   ncol = 2, nrow = 1)
  pdf("mu mvnorm(mag1_rho,mu0.2).pdf", width = 12, height = 4)
  annotate_figure(obj,
                  top = text_grob("mu mvnorm(mag1_rho,mu0.2).pdf", size = 14))
  #print(obj)
  dev.off()
}
  
  data<-data.frame(x=point[,1], y=point[,2],
                   pv = pv.t1, pws = pws,
                   pis.hat = pis.hat,
                   X = X, pluszero = pluszero,
                   p.est = p.est, identify = identify)
p <- ggplot(data, aes(x=pv.t1,y=pws,color=pluszero))+
  geom_point()
p

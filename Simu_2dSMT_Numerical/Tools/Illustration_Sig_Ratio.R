#' We are going to find a way of choosing the optimal radius for k for the neighbor
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
n <- 1 
r <- 0.9
k <- 1
rho.eps <- 0.05

#=== Init Background
point <- matrix(c(rep(seq(0,1,length.out=m.row),each=m.row),
                  rep(seq(0,1,length.out=m.row),times=m.row)), 
                m, 2, byrow = F) 

location <- (point[,1]-1/2)^2+(point[,2]-1/2)^2<(1/4)^2
mu <- runif(nrow(point),0,4)*location

pluszero <- mu>0
position <- pluszero

#Initial Dist, Sigma and X
Dist <- Dist(m)
Sigma.eps <- Sigma.eps(m)
X <- MASS::mvrnorm(n = n, mu = mu, 
                   Sigma = Sigma.eps)
X <- matrix(X, nrow = n)

#==== Save Picture
pdf("SelectNeigh.pdf", width = 10, height = 8)

# Radius.h Based
# Test the relationship between 
# the average number of neighbors and the definition of neighbor
p.len <- dim(Dist)[1]
h.len <- 100
h.test <- seq(0,1,length.out = h.len)
#position <- selected.1
sig.ratio <- sapply(h.test, 
                    function(h){
                      postive.sig <- sum((rowSums(Dist[position,position]<=h)-1)/
                                           ( rowSums(Dist[position,]<=h)-1))
                      negative.sig <- sum((rowSums(Dist[!position,!position]<=h)-1)/
                                            (rowSums(Dist[!position,]<=h)-1))
                      sig.ratio <- postive.sig+negative.sig
                      return(sig.ratio/p.len)
                    })

data <- data.frame(h.test = h.test,
                   sig.ratio = sig.ratio)
ggplot(data,aes(x=h.test,y=sig.ratio))+
  geom_point()+
  labs(title = "Signal Ratio")


# Non-Parametric Bayes Accuracy Based
# Test the relationship between 
# the efficient size of non-parametric Bayes and the definition of neighbor

## Top.k based
h.test <- 1:50
emprical.n.top <- sapply(h.test, 
                     function(hh){
                       S.neigh <- lapply(1:m,
                                         function(i){
                                           ind1 <- order(Dist[i,])[2:(hh+1)]
                                           ind1}) 
                       
                       ind <- Neigh_Partition(S.neigh)
                       return(length(ind))
                     })
qplot(x=h.test, y=emprical.n.top, main = "Emprical n(Rad.h)")

## Rad.h Based
h.test <- seq(0,1,length.out = h.len)
emprical.n.rad <- sapply(h.test, 
                     function(hh){
                       S.neigh <- lapply(1:m,
                                         function(i){
                                           ind1 <- which(Dist[i,]<=hh)
                                           ind <- ind1[ind1!=i]
                                           ind})
                       ind <- Neigh_Partition(S.neigh)
                       return(length(ind))
                     })
qplot(x=h.test, y=emprical.n.rad, main = "Emprical n(Rad.h)")

## Combine sig.ratio and rad.h
library(reshape)
data <- data.frame(h.test = h.test,
                   Sig_divide_n = sig.ratio/emprical.n.rad,
                   Sig_divide_sq.n = sig.ratio/sqrt(emprical.n.rad))
data.plot <- melt(data, id="h.test")

ggplot(data.plot,aes(x=h.test,y=value))+
  geom_point()+
  facet_grid("~variable", scales = "free")

## Try Fan's method
h.test <- 1:50
Significance.test <- sapply(h.test, 
                     function(hh){
                       Neigh_Detect_res <- Neigh_Detect(hh = hh,
                                                        X = X, 
                                                        Dist = Dist,
                                                        Sigma.eps = Sigma.eps,
                                                        detect.m = "top.k")
                       
                       #ind: selected index
                       ind <- Neigh_Detect_res$ind
                       
                       T2.select <- Neigh_Detect_res$T2[ind]
                       V2.select <- Neigh_Detect_res$V2[ind]
                       
                       Significance <- sum((T2.select/V2.select)^2-1)/sqrt(length(ind))
                       
                       return(Significance)
                     })

data.sig <- data.frame(h.test = h.test,
                   Significance.test = Significance.test,
                   emprical.n.top = emprical.n.top)

ggplot(data.sig)+
  geom_point(aes(x=h.test, y=Significance.test, color = "black")) +
  geom_point(aes(x=h.test, y=emprical.n.top, color = "red"))

dev.off()

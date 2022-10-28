#' Simulation for one case

## Start of problem independent section
options(digits=5)
require(qvalue)
require(CompRandFld)
require(np)
library("doMC")
library(fields)
setwd('~/project/2DSpatial')
#setwd("D:/RUC/project/multiple\ testing/2DSpatial")
source('Spatial_Detection_ma.R')
source('Tools/Basic.R')
source('Tools/laws_used.R')
source('Tools/Fix_GenData.R')
source('Tools/Simu_Step_ma.R')
source('Tools/Adapt_pis.R')
source('Tools/Hard_Group.R')
source('Tools/Init_Setting.R')
source('Tools/All_q_est_functions.R')
source('Tools/accumulation_test_functions.R')
source('Tools/DirectOptim/FitComposite2.R')

library(CAMT)
library(IHW)
library(adaptMT)
library(mgcv)
library(FDRreg)
library(dbh)
#source("https://www.stat.uchicago.edu/~rina/sabha/All_q_est_functions.R")
#source('https://www.stat.uchicago.edu/~rina/accumulationtests/accumulation_test_functions.R')
set.seed(10)

source('bash/read_bash.R')
#Initialize Parameter (For Testing)
q <- 0.1
n <- 1 # sample size
data.type <- "SquareSpline"

#Initialize point
m <- 900 # point size
m_x <- 30
m_y <- 30
point_x <- seq(0,5,length.out=m_x)
point_y <- seq(0,5,length.out=m_x)
point <- cbind(rep(point_x,times=length(point_y)), 
               rep(point_y,each=length(point_x))) 

#Initial Dist and Sigma
#Dist.p <- Dist(m)
Dist.p <- as.matrix(dist(point))
#Setting
h <-c(1,2,3,4,7,10,13,16) #c(1,5,9)
h.len <- length(h)
detect.m <- "top.k"
cor.type <- "cor"
k <- 1
rho.eps <- 0.1
cores <- 33
reptime <- 50
const <- q
registerDoMC(cores = cores)

fdp_res <- NULL
pow_res <- NULL
#R.total <- rep(R, each = reptime * length(Magnitude))
#Magnitude.total <- rep(rep(Magnitude, each = reptime), 
#                       times = length(R))

#for(r in R){
#  print(r)
#  Sigma.eps.p <- Sigma.eps(m, r, k, rho.eps)
#  for(magnitude in Magnitude){
I_S <- Init_Setting_2D(mu_type = mu_type, 
                       Cov_type = Cov_type,
                       magnitude = magnitude,
                       #mu_gen_machine = "uc.unif",
                       point = point,
                       Dist.p = Dist.p)
mu <- I_S$mu
Sigma.eps.p <- I_S$Sigma.eps.p
n <- 3
#estcov <- F

#foldername = paste0("Result/Simulation/2d_smoothing_k2(size",n," ",
save.folder.name <- "Simulation2D"
foldername = paste0("Result/",save.folder.name)
if(!dir.exists(foldername)){
  dir.create(foldername)
}
save.image(file="Result/tmp_pre.RData")
print("Start Simu......")
rr <- foreach(jj = 1:reptime,
              .combine = cbind) %dopar% one_step_2D(h, 
                                                    detect.m = detect.m, 
                                                    seed = jj,
                                                    mu = mu,
                                                    const = const,
                                                    magnitude = magnitude,
                                                    Cov_type = Cov_type,
                                                    mu_type = mu_type)
save.image(file="Result/tmp.RData")
seed <- numeric(reptime)
for(i in 1:reptime){
  fdp_res <- rbind(fdp_res, rr[[3*i-2]])
  pow_res <- rbind(pow_res, rr[[3*i-1]])
  seed <- rr[[3*i]]
}

pre.name <- c("BH","LAWS","SABHA",
              "AdaMT","CAMT","FDRreg(T)",
              #"FDRreg(E)",
              "dBH(True)","dBH(m_obs)","dBH(s_obs)","IHW","IHW(NULL)",
              "1D","1D.laws","1D.sabha", 
              "1D.pis2","1D.ihw","1D.ihw.null")
inside.name <- c("2D ","2D.laws ","2D.sabha ",
                 "2D.pis2 ","2D.ihw ","2D.ihw.null ")
fdp_pow_print <- rbind(apply(fdp_res,2,mean),apply(pow_res,2,mean))
if(is.null(h)){
  colnames(fdp_pow_print) <- pre.name
}else{
  if(F){
    colnames(fdp_pow_print) <- c(pre.name,
                                 paste("1D",inside.name),
                                 rep(paste(inside.name,
                                       rep(round(h,2),
                                           each=h.len)),3))
  }
  colnames(fdp_pow_print) <- colnames(fdp_res)
}
rownames(fdp_pow_print) <- c("FDP","POWER")
print(fdp_pow_print)
#  }
#}


foldername = paste0(foldername,"/2d_smoothing_k1(size",n,"const",const," est_cov",estcov)
if(!dir.exists(foldername)){
  dir.create(foldername)
}

filename = paste0(foldername, "/","mag_",magnitude," mu_",mu_type," Cov_",Cov_type," ",
                  gsub(":", "-", Sys.time()),".RData",sep="")
save.image(file = filename)



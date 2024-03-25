#' Simulation for one case

## Start of problem independent section
options(digits=5)

setwd("../")
source("Simulation/Simu_Init.R")

source('bash/read_bash.R')
#Initialize Parameter (For Testing)
set.seed(10) ## Fix unif
n <- 1 # sample size
data.type <- "SquareSpline"


detect.m <- "top.k"
cor.type <- "cor"
registerDoMC(cores = cores)

#Setting
h=4
h.len <- length(h)

for(m in seq(100,2000,100)){
  
  point <- matrix(seq(0,30/900*m,length.out=m), 
                  m, 1, byrow = F) 
  Dist.p <- as.matrix(dist(point))
#R.total <- rep(R, each = reptime * length(Magnitude))
#Magnitude.total <- rep(rep(Magnitude, each = reptime), 
#                       times = length(R))

#for(r in R){
#  print(r)
#  Sigma.eps.p <- Sigma.eps(m, r, k, rho.eps)
#  for(magnitude in Magnitude){
I_S <- Init_Setting_1D(mu_type = "Medium", 
                       Cov_type = "Medium",
                       magnitude = magnitude,
                       #mu_gen_machine = "uc.unif",
                       #mu_gen_machine = "mvnorm",
                       point = point,
                       Dist.p = Dist.p)
mu <- I_S$mu
Sigma.eps.p <- I_S$Sigma.eps.p
n <- 1
estcov <- F

save.folder.name <- paste0("Simulation1D_CalTime")
#foldername = paste0("Result/Simulation/2d_smoothing_k2(size",n," ",
#save.folder.name <- "Simulation1D_unif"
#save.folder.name <- "Simulation1D_mv"
foldername = paste0("Result/",save.folder.name)
if(!dir.exists(foldername)){
  dir.create(foldername)
  dir.create(paste0(foldername,"/Tmp"))
}
save.image(file="Result/tmp_pre.RData")
print("Start Simu......")
    rr <- foreach(jj = 1:reptime,
                  .combine = cbind) %dopar% one_step_1D.CalTime(h, 
                                                     detect.m = detect.m, 
                                                     seed = jj,
                                                     mu = mu,
                                                     const = const,
                                                     magnitude = magnitude,
                                                     Cov_type = Cov_type,
                                                     EmpMethod=EmpMethod,
                                                     mu_type = mu_type,
                                                     tau.tm = tau.tm,
                                                     save.folder.name = save.folder.name)
    save.image(file="Result/tmp.RData")
    seed <- numeric(reptime)
    


foldername = paste0(foldername,"/2d_smoothing_k1(size",n,"const",const," est_cov",estcov)
if(!dir.exists(foldername)){
  dir.create(foldername)
}

filename = paste0(foldername, "/","mag_",magnitude," m_",m," mu_",mu_type," Cov_",Cov_type," ",
                  gsub(":", "-", Sys.time()),".RData",sep="")
save.image(file = filename)
}


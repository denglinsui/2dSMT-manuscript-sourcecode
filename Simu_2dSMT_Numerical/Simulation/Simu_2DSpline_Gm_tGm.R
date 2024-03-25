#' Simulation for one case

## Start of problem independent section
options(digits=5)

setwd("../")
source("Simulation/Simu_Init.R")
source("Tools/DirectOptim/FitComposite2.R")
source('bash/read_bash.R')
#Initialize Parameter (For Testing)
set.seed(10)

q <- 0.1
n <- 1 # sample size
data.type <- "SquareSpline"

registerDoMC(cores = 24)
reptime <- 200
Err_res <- NULL
#Initialize point
for(m in seq(10,50,by=5)^2){
#m <- 900 # point size
m_each <- sqrt(m)
m_x <- m_each
m_y <- m_each

point_x <- seq(0,5/30*m_x,length.out=m_x)
point_y <- seq(0,5/30*m_y,length.out=m_x)
point <- cbind(rep(point_x,times=length(point_y)), 
               rep(point_y,each=length(point_x))) 

#Initial Dist and Sigma
#Dist.p <- Dist(m)
Dist.p <- as.matrix(dist(point))
#Setting
h <-c(4) #c(1,5,9)
h.len <- length(h)
detect.m <- "top.k"
cor.type <- "cor"
#R.total <- rep(R, each = reptime * length(Magnitude))
#Magnitude.total <- rep(rep(Magnitude, each = reptime), 
#                       times = length(R))

#for(r in R){
#  print(r)
#  Sigma.eps.p <- Sigma.eps(m, r, k, rho.eps)
#  for(magnitude in Magnitude){

#estcov <- F

#foldername = paste0("Result/Simulation/2d_smoothing_k2(size",n," ",
save.folder.name <- "Simulation2D_GMLE"
foldername = paste0("Result/",save.folder.name)
if(!dir.exists(foldername)){
  dir.create(foldername)
  dir.create(paste0(foldername,"/Tmp"))
}
save.image(file="Result/tmp_pre.RData")
print("Start Simu......")
rr <- foreach(jj = 1:reptime,
              .combine = cbind) %dopar% one_step_2D_Gm_tGm_comp(h, 
                                                       detect.m = detect.m, 
                                                       seed = jj,
                                                       #mu = mu,
                                                       const = const,
                                                       magnitude = magnitude,
                                                       Cov_type = Cov_type,
                                                       EmpMethod=EmpMethod,
                                                       tau.tm = tau.tm,
                                                       mu_type = mu_type)
save.image(file="Result/tmp.RData")

Err_res <- rbind(Err_res,data.frame(t(rr),m = m))


filename = paste0(foldername, "/Tmp/","mag_",magnitude," mu_",mu_type," Cov_",
                  Cov_type," m", m," ",
                  gsub(":", "-", Sys.time()),".RData",sep="")
save.image(file = filename)

}
foldername = paste0(foldername,"/GMLE_Gm_tGm")
if(!dir.exists(foldername)){
  dir.create(foldername)
}

filename = paste0(foldername, "/","mag_",magnitude," mu_",mu_type," Cov_",Cov_type," ",
                  gsub(":", "-", Sys.time()),".RData",sep="")
save.image(file = filename)



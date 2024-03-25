#' Simulation for one case

## Start of problem independent section
options(digits=5)

setwd("../")
source("Simulation/Simu_Init.R")

source('bash/read_bash.R')
#Initialize Parameter (For Testing)
set.seed(10) ## Fix unif
n <- 1 # sample size
data.type <- "SquareSpline_Spa_Grp"

#Initialize point
m <- 900 # point size

point <- matrix(seq(0,30,length.out=m), 
                m, 1, byrow = F) 

h <-4#c(1,2,3,4,7,10,13,16) #c(1,5,9)

##=======New Setting for m=2000======
#m <- 2000
if(m == 2000){
  point <- matrix(seq(0,60,length.out=m), 
                  m, 1, byrow = F) 
  h <- 4
}

#Initial Dist and Sigma
#Dist.p <- Dist(m)
Dist.p <- as.matrix(dist(point))
#Setting
h.len <- length(h)
detect.m <- "top.k"
cor.type <- "cor"
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
I_S <- Init_Setting_1D_spa_grp(
                       Cov_type = Cov_type,
                       magnitude = magnitude,
                       #mu_gen_machine = "uc.unif",
                       #mu_gen_machine = "mvnorm",
                       point = point,
                       Dist.p = Dist.p)

mu <- I_S$mu
Sigma.eps.p <- I_S$Sigma.eps.p
grp <- I_S$grp
n <- 1
estcov <- F

save.folder.name <- paste0("Simulation1D_spa_grp","_m",m)
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
              .combine = cbind) %dopar% one_step_1D_spa_grp(h, 
                                                            grp = grp,
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
for(i in 1:reptime){
  fdp_res <- rbind(fdp_res, rr[[3*i-2]])
  pow_res <- rbind(pow_res, rr[[3*i-1]])
  seed <- rr[[3*i]]
}


pre.name <- c("1D",
              "1D.pis2","1D.pis2.Grp.SA","1D.pis2.Grp.LAWS")
inside.name <- c("2D ","2D.pis2 ",
                 "2D.pis2.Grp.SA ","2D.pis2.Grp.LAWS ")
fdp_pow_print <- rbind(apply(fdp_res,2,mean),apply(pow_res,2,mean))
if(is.null(h)){
  colnames(fdp_pow_print) <- pre.name
}else{
  if(F){
    colnames(fdp_pow_print) <- c(pre.name,
                                 paste("1D",inside.name),
                                 paste(inside.name,
                                       rep(round(h,2),
                                           each=h.len)))
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



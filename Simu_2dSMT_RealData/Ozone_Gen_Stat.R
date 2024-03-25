## Before run the simulation for Ozone data, please set the dictionary to 2dSMT 
##' Ozone Analysis
#-- Manipulate data
library(reshape2)
library(data.table)
library(dplyr)
library(data.table)
library(backports)
library(CompRandFld)
#-- For detect algorithms
library(qvalue)
load("Simu_2dSMT_RealData/Data/Ozone.RData")

#=== Regression
# Obtain beta_hat/sd(beta_hat)
#'@param o.level The level of ozone averaged by year
#'@param year The year
#'@return the test statistics for the decreasing coeff.
test.stat <- function(o.level,year){
  # measure decreasing
  beta0 <- -0.1; y <- o.level*1000 # ppm -> ppb
  x <- year
  
  # fit
  fit <- lm(y~x)
  coeff.tmp <- summary(fit)$coefficients
  
  # obtain test.stat
  beta.hat <- coeff.tmp[,"Estimate"][2]
  beta.hat.sd <- coeff.tmp[,"Std. Error"][2]
  test.stat <- (beta.hat-beta0)/beta.hat.sd
  return(test.stat)
}


mu.hat <- function(o.level,year){
  # measure decreasing
  beta0 <- -0.1; y <- o.level*1000 # ppm -> ppb
  x <- year
  
  # fit
  fit <- lm(y~x)
  coeff.tmp <- summary(fit)$coefficients
  
  # obtain test.stat
  mu.hat <- coeff.tmp[,"Estimate"][1]
  return(mu.hat)
}

beta.hat <- function(o.level,year){
  # measure decreasing
  beta0 <- -0.1;y <- o.level*1000 # ppm -> ppb
  x <- year
  
  # fit
  fit <- lm(y~x)
  coeff.tmp <- summary(fit)$coefficients
  
  # obtain test.stat
  beta.hat <- coeff.tmp[,"Estimate"][2]
  return(beta.hat)
}

hat.sd <- function(o.level,year){
  # measure decreasing
  beta0 <- -0.1;y <- o.level*1000 # ppm -> ppb
  x <- year
  
  # fit
  fit <- lm(y~x)
  hat.sd <- summary(fit)$sigma
  return(hat.sd)
}

beta.hat.sd <- function(o.level,year){
  # measure decreasing
  beta0 <- -0.1;y <- o.level*1000 # ppm -> ppb
  x <- year
  
  # fit
  fit <- lm(y~x)
  coeff.tmp <- summary(fit)$coefficients
  
  # obtain test.stat
  beta.hat.sd <- coeff.tmp[,"Std. Error"][2]
  return(beta.hat.sd)
}

res.stat <- function(o.level,year){
  # measure decreasing
  beta0 <- -0.1; y <- o.level*1000 # ppm -> ppb
  x <- year
  
  # fit
  fit <- lm(y~x)
  
  # obtain res.stat
  res.stat <- fit$residuals
  return(res.stat)
}

Res_Stat <- Annual_Summary[,
                           .(res.stat = res.stat(Avr,Year),
                             hat.sd = hat.sd(Avr,Year),
                             Year = Year,
                             State.Code=State.Code),
                           by = .(Latitude,Longitude)]

Res_Stat <- Res_Stat[,.(Year,Latitude,Longitude,hat.sd,State.Code,
                        res.stat =res.stat/hat.sd)]

Res_mat <- acast(Res_Stat, Year~Latitude+Longitude, value.var="res.stat")
Res_sd <- acast(Res_Stat, Year~Latitude+Longitude, value.var="hat.sd")[1,]

#=== Estimate covariance
corrmodel <- "exponential"
year.num <- nrow(Res_mat)
loc.num <- ncol(Res_mat)

save.image("Tmp.RData")
est.ind <- 1:loc.num
point.for.cov <- matrix(as.numeric(as.matrix(list2DF(sapply(colnames(Res_mat),
                                                            function(x){strsplit(x,split="_")})))),byrow=T,ncol=2)
est.ind <- 1:loc.num
rep.ind <- 1:year.num
rem <- c(10,0.5,0.5)
names(rem) <- c("scale","sill","nugget")
fit <- FitComposite(Res_mat[rep.ind,est.ind], coordx=point.for.cov[est.ind,],
                    corrmodel=corrmodel, likelihood="Full",type="Standard", 
                    #corrmodel=corrmodel, likelihood='Marginal',type='Pairwise',#optimizer = "CG",  
                    fixed = list(mean=0), start = as.list(rem),
                    distance = "Geod",replicates = length(rep.ind))

param <- as.list(fit$param)

#--- Estimating the covariance matrix
cov.est <- Covmatrix(point.for.cov[,1], point.for.cov[,2], 
                     corrmodel=corrmodel,
                     distance = "Geod",
                     param=param)
Corr.est <- cov.est$covmatrix 


#--- Estimating the covariance matrix with gaussian as correlation model
fit.gau <- FitComposite(Res_mat[rep.ind,est.ind], coordx=point.for.cov[est.ind,],
                        corrmodel="gauss", likelihood="Full",type="Standard", 
                        #corrmodel=corrmodel, likelihood='Marginal',type='Pairwise',#optimizer = "CG",  
                        fixed = list(mean=0), start = as.list(rem),
                        distance = "Geod",replicates = length(rep.ind))

param.gau <- as.list(fit.gau$param)

cov.est.gau <- Covmatrix(point.for.cov[,1], point.for.cov[,2], 
                         corrmodel="gauss",
                         distance = "Geod",
                         param=param.gau)
Corr.est.gau <- cov.est.gau$covmatrix 
#Sigma.eps.est.diff[1:5,1:5]

#--- Estimating the empirical covariance matrix
Corr.est.emp <- ((t(Res_mat) %*% Res_mat)/(year.num-1-1))

#--- Obtain the covariance matrix for beta.hat
TT <- cbind(rep(1,12),2010:2021)
inv_tTT <- solve(t(TT)%*%TT)[2,2]
Corr.beta.est <- Corr.est*inv_tTT
Sigma.beta.est <- Corr.beta.est* (Res_sd%*%t(Res_sd))
#=== Finish estimating the covariance matrix for beta.hat

#---- Obtain Test Statitics
Beta_Stat <- Annual_Summary[,.(beta_hat = beta.hat(Avr,Year),mu_hat=mu.hat(Avr,Year),State.Code=unique(State.Code)),
                            by = .(Latitude,Longitude)]
Sd_Stat <- as.data.table(point.for.cov); colnames(Sd_Stat) <- c("Latitude","Longitude")
Sd_Stat$beta_sd <- sqrt(diag(Sigma.beta.est))

beta0 <- -0.5
T2_Stat <- merge(Beta_Stat,Sd_Stat) 
T2_Stat <- T2_Stat[,.(Latitude,Longitude,beta_hat,beta_sd,mu_hat,State.Code,T2 = (beta_hat-beta0)/beta_sd)]
#Corr.est is the covariance matrix for T2
save.image("Simu_2dSMT_RealData/Data/Ozone_Tmp.RData")

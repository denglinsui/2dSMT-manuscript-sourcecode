# Source necessary files
require(qvalue)
if(system.file(package='CompRandFld')!=""){
  require(CompRandFld)
}
if(system.file(package='fields')!=""){
  require(fields)
}
require(np)
library("doMC")

#setwd('~/project/2dSMT/Simu_2dSMT_Numerical')
library(TwoDSMT)
library(REBayes)
library(distrEx)
#setwd("D:/RUC/project/multiple\ testing/2DSpatial")
#source('Spatial_Detection_ma.R')
source('Tools/Basic.R')
source('Tools/laws_used.R')
source('Tools/Fix_GenData.R')
source('Tools/Simu_Step_ma.R')
source('Tools/Simu_Step_ma_Gm.R')
source('Tools/Simu_Step_ma_Gm_tGm.R')
source('Tools/Simu_Step_ma_spa_grp.R')
source('Tools/Simu_Step_ma_CalTime.R')
source('Tools/Simu_Step_ma_num_neigh.R')
source('Tools/Simu_Step_ma_uneven.R')
#source('Tools/Simu_Step_ma_NPEB.R')
source('Tools/Adapt_pis.R')
#source('Tools/Hard_Group.R')
source('Tools/Init_Setting.R')
source('Tools/All_q_est_functions.R')
source('Tools/accumulation_test_functions.R')

#library(CAMT)
library(IHW)
library(adaptMT)
library(mgcv)
#library(FDRreg)
library(dbh)
library(pbivnorm)


reptime <- 100
q <- 0.1
const <- q
cores <- 30
tau.tm <- q
EmpMethod <- "NPEB"#"DeCor"

#source("https://www.stat.uchicago.edu/~rina/sabha/All_q_est_functions.R")
#source('https://www.stat.uchicago.edu/~rina/accumulationtests/accumulation_test_functions.R')


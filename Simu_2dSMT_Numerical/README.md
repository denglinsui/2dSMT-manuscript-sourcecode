# 2DSMT (Simulation)

This project incorporates the simulation codes for "Powerful Spatial Multiple Testing via Borrowing Neighboring Information".

## Packages

To successfully reproduce the results, please check whether the following R packages are available:

```R
library(latex2exp)
library(ggplot2)
library(ggpubr)
library(RBayes)
library(Rmosek)
library(np)
library(CompRandFld)
library(doMC)
library(argparse)

# Competing Methods
library(qvalue)
library(CAMT)
library(IHW)
library(adaptMT)
library(mgcv)
library(FDRreg)
library(dbh)
```

To install `Rmosek`, please visit https://docs.mosek.com/latest/rmosek/install-interface.html. 

## Simulation

* To reproduce the simulation results in Section 5, input the following codes in the terminal:

```bash
nohup bash ../run_1DSpline.sh &
nohup bash ../run_1Dunif.sh &
nohup bash ../run_1Dmvnorm.sh &
nohup bash ../run_2DSpline.sh &
```

* To create Figures 2 to 9 in Section 5 of the main paper, run `Print_2D.Rmd` and `Print_Result_single.Rmd` (Please also modify the root). 
* To create  Figures S.4-S.9 in Section S.I.2 of the supplementary materials, run `Print_Result_vary_h_from_single.Rmd`  (Please also modify the root). 

## Other Figures

* To create Figures 1 in the main paper, run `Tools/Selected_Cutoff_Simple.R`.  
* To create Figures S.1 to S.3 in the supplementary materials, run  `Tools/SimulationSetting.R`.
* To create Figures S.13 in the supplementary materials, run  `Tools/Reduced_Set.R`.
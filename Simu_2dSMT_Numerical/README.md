# 2DSMT (Simulation)

This project incorporates the simulation codes for "Powerful Spatial Multiple Testing via Borrowing Neighboring Information".

## Packages

To successfully reproduce the results, please check whether the following R packages are available:

```R
library(latex2exp)
library(ggplot2)
library(ggpubr)
library(REBayes)
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

* ## Simulations in Sections 5, S.I.2 and S.I.3

  * To reproduce the simulation results in Section 5, input the following codes in the terminal:

  ```bash
  nohup bash run_1DSpline.sh &
  nohup bash run_1Dunif.sh &
  nohup bash run_1Dmvnorm.sh &
  nohup bash run_2DSpline.sh &
  ```

  * To create Figures 2 to 5 in Section 5 of the main paper, run `Print_2D.Rmd` and `Print_Result_single.Rmd` in the folder `Print_Result` (Please also modify the root). 

  * To create Figures S.4 to S.6 in Section S.I.2 of supplement:

    * uncomment `#m <- 2000` in `Simu_1DSpline.R`, `Simulation_1Dmvnorm.R`, and `Simulation_1dunif.R` in the folder `Simulation`.  
    * Rerun the bash files `run_1DSpline.sh`, `run_1Dunif.sh` and `run_1Dmvnorm.sh`.
    * Run `Print_Result_single.Rmd` in the folder `Print_Result` (Please also modify the root). 

  * To create  Figures S.7-S.12 in Section S.I.3 of the supplement, run `Print_Result_vary_h_from_single.Rmd` in the folder `Print_Result` (Please also modify the root). 

    

* ## Other simulations in the supplement

  * To reproduce the simulation results in Section S.I.4

    * Run  the following codes in the terminal: 

      ```bash
      nohup bash run_1DSpline_neigh.sh &
      nohup bash run_1DUnif_neigh.sh &
      nohup bash run_1Dmv_neigh.sh &
      ```

      

    * To create Figures S.13 and S.14, run `Print_Result_neigh.Rmd`  in the folder `Print_Result`.

  * To reproduce the simulation results in Section S.I.5

    * Run `nohup bash run_1DSpline_spa_grp.sh &` in the terminal;
    * To create Figure S.16, run `Print_Result_spa_grp.Rmd`  in the folder `Print_Result`.


  * To reproduce the simulation results in Section S.VII
    * Run `nohup bash run_1DSpline_CalTime.sh &` in the terminal;
    * To create Figure S.24, run `Print_Result_CalTime.Rmd`  in the folder `Print_Result`.

* ## Other Figures

  * To create Figures 1 in the main paper, run `Tools/Selected_Cutoff_Simple.R`.  
  * To create Figures S.1, S.2, S.3 and  S.15 in the supplement, run  `Tools/SimulationSetting.R`.
  * To create Figures S.20 to S.22 in the supplement, run  `Tools/IlluPowerImp.R`.
  * To create Figures S.23 in the supplement, run  `Tools/Reduced_Set.R`.

* ## More Simulations

  * You can use`run_2D_GMLE.sh` and  `Print_Gm_tGm.R` in `Print_Result` to investigate the convergence of GMLE.
  * You can use`run_1DSpline_uneven.sh.sh` and  `Print_Result_single_uneven.Rmd` in `Print_Result` to study the performance of 2d-SMT when the locations are unevenly distributed and the neighorhood is chosen by setting a threshold for the distance.
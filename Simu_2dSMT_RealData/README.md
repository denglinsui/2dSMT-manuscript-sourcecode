# 2DSMT (Real Data Analysis)

This project incorporates the codes for real data analysis in "Powerful Spatial Multiple Testing via Borrowing Neighboring Information".

## Packages

To successfully reproduce the results, please check whether the following R packages are available:

```R
library(IHW)
library(TwoDSpatial)
library(qvalue)
library(reshape2)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggmap)
library(mapproj)
library(data.table)
library(backports)
library(CompRandFld)
library(cowplot)
```

To install `Rmosek`, please visit https://docs.mosek.com/latest/rmosek/install-interface.html. 

To install `TwoDSpatial`, run the following code in the terminal:

```bash
R CMD INSTALL TwoDSpatial_0.0.0.9000.tar.gz
```



## Data Preparation

* Download data from http://www.epa.gov/airexplorer/index.htm and save in `Data/Ozone` (You should obtain the files from  `annual_conc_by_monitor_2000.csv` to `annual_conc_by_monitor_2021.csv` ).
* Run `LoadOzone.R` to prepare data for Ozone.
* Run `Ozone_Gen_Stat.R` to obtain test statistis for identifying the annual increases of ozone.
* Run `Ozone_Visual.R` to reprocedure Figure S.10 in the supplementary material.

## Real Data Analysis

* Run `Ozone_2D.R` to conduct hypothesis testing ` (Also load `Tools/All_q_est_functions.R` in simulation codes.).
* Run `OzoneCONO.R` to reproduce Figure 10 in the main paper.
* Run `Ozone_2D_Extract_Info.R` to obatin the information of Tables S.1 and S.2 in the supplementary material.
* Run `OzoneSimu.R` to reproduce Table 1 in the main paper.


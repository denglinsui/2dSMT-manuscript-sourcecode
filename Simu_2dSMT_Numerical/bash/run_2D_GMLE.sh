#!/bin/bash

declare Mu_t=("Dense")
declare Cov_t=("Strong")
for mag in $(seq 0.5 0.5 1.5)
do
for (( i=0; i<=0; i++ ))
do
for (( j=0; j<=0; j++ ))
do
echo "1DSpline_h:","mag", $mag, "mu_t", ${Mu_t[$i]}, "cov_t", ${Cov_t[$j]}

nohup Rscript ../Simulation/Simu_2DSpline_Gm_tGm.R --mu_type ${Mu_t[$i]} --Cov_type ${Cov_t[$j]} --magnitude $mag --est_cov True >../logerr/Spline2DGMLE$mag${Mu_t[$i]}${Cov_t[$j]}.out 2>../logerr/Spline2DGMLE$mag${Mu_t[$i]}${Cov_t[$j]}.err

done
done
done


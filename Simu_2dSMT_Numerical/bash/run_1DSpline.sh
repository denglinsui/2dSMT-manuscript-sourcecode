#!/bin/bash
#for magnitude in $(seq 0.1 0.5 10)

declare Mu_t=("Sparse" "Medium" "Dense")
declare Cov_t=("Weak" "Medium" "Strong")
for mag in $(seq 2 1 4)
do
for (( i=0; i<=2; i++ ))
do
for (( j=0; j<=2; j++ ))
do
echo "1DSpline:","mag", $mag, "mu_t", ${Mu_t[$i]}, "cov_t", ${Cov_t[$j]}
#nohup Rscript ../Simulation/Simu_1DSpline_h.R --mu_type ${Mu_t[$i]} --Cov_type ${Cov_t[$j]} --magnitude $mag >../logerr/Spline1D_h_t.out 2>../logerr/Spline1D_h_t.err
#nohup Rscript ../Simulation/Simu_1DSpline.R --mu_type ${Mu_t[$i]} --Cov_type ${Cov_t[$j]} --magnitude $mag --est_cov False >../logerr/Spline1D$mag${Mu_t[$i]}${Cov_t[$j]}.out 2>../logerr/Spline1D$mag${Mu_t[$i]}${Cov_t[$j]}.err
nohup Rscript ../Simulation/Simu_1DSpline.R --mu_type ${Mu_t[$i]} --Cov_type ${Cov_t[$j]} --magnitude $mag --est_cov False >../logerr/Spline1D$mag${Mu_t[$i]}${Cov_t[$j]}.out 2>../logerr/Spline1D$mag${Mu_t[$i]}${Cov_t[$j]}.err
done
done
done

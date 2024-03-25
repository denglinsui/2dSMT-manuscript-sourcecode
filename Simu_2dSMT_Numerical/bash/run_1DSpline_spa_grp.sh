#!/bin/bash

declare Mu_t=("Dense")
declare Cov_t=("Weak" "Medium" "Strong")
for mag in $(seq 2 1 4)
do
for (( i=0; i<=0; i++ ))
do
for (( j=0; j<=2; j++ ))
do
echo "1DSpline_spa_grp:","mag", $mag, "mu_t", ${Mu_t[$i]}, "cov_t", ${Cov_t[$j]}
#nohup Rscript ../Simulation/Simu_1DSpline_h.R --mu_type ${Mu_t[$i]} --Cov_type ${Cov_t[$j]} --magnitude $mag >../logerr/Spline1D_h_t.out 2>../logerr/Spline1D_h_t.err
#nohup Rscript ../Simulation/Simu_1DSpline.R --mu_type ${Mu_t[$i]} --Cov_type ${Cov_t[$j]} --magnitude $mag --est_cov False >../logerr/Spline1D$mag${Mu_t[$i]}${Cov_t[$j]}.out 2>../logerr/Spline1D$mag${Mu_t[$i]}${Cov_t[$j]}.err
nohup Rscript ../Simulation/Simu_1DSpline_spa_grp.R --mu_type ${Mu_t[$i]} --Cov_type ${Cov_t[$j]} --magnitude $mag --est_cov False >../logerr/Spline1D_spa_grp$mag${Mu_t[$i]}${Cov_t[$j]}.out 2>../logerr/Spline1D_spa_grp$mag${Mu_t[$i]}${Cov_t[$j]}.err
#echo "1DSpline_h:","mag", $magnitude, "r", $r, "estcov", "False"
#nohup Rscript ../Simulation/Simu_1DSpline_h.R --Mag $magnitude --R $r --estcov False >../logerr/Spline1D_h_f.out 2>../logerr/Spline1D_h_f.err
done
done
done


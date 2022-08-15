# Longitudinal-trait-GWAs
This repository have an example of the parameter file, and linux script for running  random regression models to estimate SNP effects of longitudinal traits. 
For running the RRM we first run renumf90 using the renum.txt parameter file. After, we run gibbs2f90 to estimate the SNP effects using the renf90.par file generated from renumf90. 
gibb_rrm2.sh is the bash script for running gibbs2f90
postBLUPf90_GWASlongitunidanltrait.R is an R script for taking the solutions from running gibbs2f90 and estimate the SNPs effects by each day after planting. 

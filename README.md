# prom1
This repository contains MATLAB and Python scripts for the reference model of Korenaga et al. (JGR, 2021) and related files. 

**MATLAB**:
"plot_refmodel.m" & "plot_refmodel2.m" are sample scripts to show how to use a function calc_refmodel() contained
in "calc_reflmodel.m". "plot_refmodel2.m" also plots the age-depth and age-heat flow data of normal seafloor, to be
compared with the reference model; it produces something similar to Figures 7a and 8 of Korenaga et al. (JGR, 2021). 
"plot_refmodel.m" and "plot_refmodel2.m" will produce a plot like "plot1.png" and "plot2.png", respectively. 
"plot_refmodel.m" also saves reference model predictions into two ASCII files, "ref_d_q.dat" and "ref_T.dat"; the former
contains depth and heat flow as a function of age, the latter contains thermal structure as a function of age and depth. 

**Python**:
"refmodel.py", "plot_refmodel.py", and "plot_refmodel2.py" correspond to, respectively, "calc_refmodel.m", "plot_refmodel.m", and "plot_refmodel2.m"
described above. 

**Data files**:
"normal_sa_depth_hist.dat" contains age-depth relative freqency, and 
"hf_quartile_filHFnormal.dat" contains age-heat flow data in terms of quartiles. 
lister_nagihara.dat contains the data of Lister et al. (GJI, 1990) and Nagihara et al. (EPSL, 1996)
that are located on the normal seafloor. These data are compiled as described in Korenaga et al. (JGR, 2021). 

**Sample output files**:
"ref_d_q.dat" contains age (Ma), seafloor depth (m), and heat flow (mW/m^2), for 0-180 Ma with an increment of 1 Ma. 
"ref_T.dat" contains age (Ma), depth (km), and temperature (K) for 0-180 Ma with an increment of 1 Ma and 0-300 km with an increment of 1 km. 
These are created by running "plot_refmodel.m". 

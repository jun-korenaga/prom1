# prom1
This repository contains MATLAB and Python scripts for the reference model of Korenaga et al. (JGR, 2021). 

MATLAB:
plot_refmodel.m & plot_refmodel2.m are sample scripts to show how to use a function calc_refmodel() contained
in calc_reflmodel.m. plot_refmodel2.m also plots the age-depth and age-heat flow data of normal seafloor, to be
compared with the reference model; it produces something similar to Figures 7a and 8 of Korenaga et al. (JGR, 2021). 

Python:
refmodel.py, plot_refmodel.py, and plot_refmodel2.py correspond to calc_refmodel.m, plot_refmodel.m, and plot_refmodel2.m
described above. 

Data files:
normal_sa_depth_hist.dat contains age-depth relative freqency, and 
hf_quartile_filHFnormal.dat contains age-heat flow data in terms of quartiles. 
lister_nagihara.dat contains the data of Lister et al. (GJI, 1990) and Nagihara et al. (EPSL, 1996)
that are located on the normal seafloor.

These data are compiled as described in Korenaga et al. (JGR, 2021). 

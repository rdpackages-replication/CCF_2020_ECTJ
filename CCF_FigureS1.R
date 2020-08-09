################################################################################################################################################################
### CALONICO-CATTANEO-FARRELL (2019): Optimal Bandwidth Choice for Robust Bias Corrected Inference in Regression-Discontinuity Designs.
### Figure S1
################################################################################################################################################################
rm(list=ls(all=TRUE))
library(ggplot2)
setwd("...")
source("functions/CCF_functions.R")

for (p in 0:3) rho.optim(p=p, kernel="tria", plot=TRUE)
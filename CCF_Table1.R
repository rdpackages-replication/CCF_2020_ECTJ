################################################################################################################################################################
### CALONICO-CATTANEO-FARRELL (2019): Optimal Bandwidth Choice for Robust Bias Corrected Inference in Regression-Discontinuity Designs.
### Table 1
################################################################################################################################################################
rm(list=ls(all=TRUE))
library(Hmisc)
setwd("...")
source("functions/CCF_functions.R")

table1 = matrix(NA,5,4)
table1[,1] = 0:4
for (p in 0:4) table1[p+1,2] = round(rho.optim(p=p,kernel="tria",plot=FALSE),3)
for (p in 0:4) table1[p+1,3] = round(rho.optim(p=p,kernel="epan",plot=FALSE),3)
for (p in 0:4) table1[p+1,4] = round(rho.optim(p=p,kernel="unif",plot=FALSE),3)
colnames(table1) = c("","Triangular","Epanechnikov","Uniform")
table_tex = latex(table1 , file = paste("CCF_Table1",".txt",sep=""), 
                    n.cgroup = c(1,3) , cgroup = c("$p$","Kernel"),
                    landscape=FALSE, center='none', col.just=rep('c',4), title='', table.env=FALSE)
  


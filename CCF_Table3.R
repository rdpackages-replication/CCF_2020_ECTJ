################################################################################################################################################################
### CALONICO-CATTANEO-FARRELL (2019): Optimal Bandwidth Choice for Robust Bias Corrected Inference in Regression-Discontinuity Designs.
### Table 3
################################################################################################################################################################
rm(list=ls(all=TRUE))
setwd("...")
library(Hmisc)
library(Rcpp)
source("functions/CCF_rdrobust.R")
source("functions/CCF_rdbwselect.R")
source("functions/CCF_functions.R")
sourceCpp("functions/CCF_functions.cpp")

data = read.csv("headstart.csv")

# 1960 Census covariates
X60 = cbind(data$census1960_pop,data$census1960_pctsch1417,data$census1960_pctsch534, data$census1960_pctsch25plus,data$census1960_pop1417,data$census1960_pop534,data$census1960_pop25plus,data$census1960_pcturban,data$census1960_pctblack)

# Outcome, normalized running variable and treatment
cutoff = 59.1984
Y = data$mort_age59_related_postHS
R = data$povrate60
X = data$povrate60 - cutoff

rho.star = 0.850

out.mse1 = rdrobust(Y, X, bwselect = 'mserd')
out.mse2 = rdrobust(Y, X, bwselect = 'mserd', rho=rho.star)
out.mse3 = rdrobust(Y, X, bwselect = 'mserd', rho=1)

out.mse4 = rdrobust(Y, X, bwselect = 'mserd',               covs = X60)
out.mse5 = rdrobust(Y, X, bwselect = 'mserd', rho=rho.star, covs = X60)
out.mse6 = rdrobust(Y, X, bwselect = 'mserd', rho=1,        covs = X60)

out.cer1 = rdrobust(Y, X, bwselect = 'cerdpi')
out.cer2 = rdrobust(Y, X, bwselect = 'cerdpi', rho=rho.star)
out.cer3 = rdrobust(Y, X, bwselect = 'cerdpi', rho=1)

out.cer4 = rdrobust(Y, X, bwselect = 'cerrd')
out.cer5 = rdrobust(Y, X, bwselect = 'cerrd', rho=rho.star)
out.cer6 = rdrobust(Y, X, bwselect = 'cerrd', rho=1)

out.cer7 = rdrobust(Y, X, bwselect = 'cerrd', covs = X60)
out.cer8 = rdrobust(Y, X, bwselect = 'cerrd', rho=rho.star, covs = X60)
out.cer9 = rdrobust(Y, X, bwselect = 'cerrd', rho=1, covs = X60)

table.HS=matrix(NA,5,6)
table.HS[1,1] = round(out.mse1$coef[1],3)
table.HS[1,2] = round(out.mse1$bws[1],3)
table.HS[1,3] = round(out.mse1$bws[1]/out.mse1$bws[2],3)
table.HS[1,4] = paste("[",round(out.mse1$ci[3,1],2)," , ",round(out.mse1$ci[3,2],2),"]",sep="")
table.HS[1,5] = paste("[",round(out.mse2$ci[3,1],2)," , ",round(out.mse2$ci[3,2],2),"]",sep="")
table.HS[1,6] = paste("[",round(out.mse3$ci[3,1],2)," , ",round(out.mse3$ci[3,2],2),"]",sep="")

table.HS[2,1] = round(out.mse4$coef[1],3)
table.HS[2,2] = round(out.mse4$bws[1],3)
table.HS[2,3] = round(out.mse4$bws[1]/out.mse1$bws[2],3)
table.HS[2,4] = paste("[",round(out.mse4$ci[3,1],2)," , ",round(out.mse4$ci[3,2],2),"]",sep="")
table.HS[2,5] = paste("[",round(out.mse5$ci[3,1],2)," , ",round(out.mse5$ci[3,2],2),"]",sep="")
table.HS[2,6] = paste("[",round(out.mse6$ci[3,1],2)," , ",round(out.mse6$ci[3,2],2),"]",sep="")

table.HS[3,1] = round(out.cer1$coef[1],3)
table.HS[3,2] = round(out.cer1$bws[1],3)
table.HS[3,3] = round(out.cer1$bws[1]/out.mse1$bws[2],3)
table.HS[3,4] = paste("[",round(out.cer1$ci[3,1],2)," , ",round(out.cer1$ci[3,2],2),"]",sep="")
table.HS[3,5] = paste("[",round(out.cer2$ci[3,1],2)," , ",round(out.cer2$ci[3,2],2),"]",sep="")
table.HS[3,6] = paste("[",round(out.cer3$ci[3,1],2)," , ",round(out.cer3$ci[3,2],2),"]",sep="")

table.HS[4,1] = round(out.cer4$coef[1],3)
table.HS[4,2] = round(out.cer4$bws[1],3)
table.HS[4,3] = round(out.cer4$bws[1]/out.mse1$bws[2],3)
table.HS[4,4] = paste("[",round(out.cer4$ci[3,1],2)," , ",round(out.cer4$ci[3,2],2),"]",sep="")
table.HS[4,5] = paste("[",round(out.cer5$ci[3,1],2)," , ",round(out.cer5$ci[3,2],2),"]",sep="")
table.HS[4,6] = paste("[",round(out.cer6$ci[3,1],2)," , ",round(out.cer6$ci[3,2],2),"]",sep="")

table.HS[5,1] = round(out.cer7$coef[1],3)
table.HS[5,2] = round(out.cer7$bws[1],3)
table.HS[5,3] = round(out.cer7$bws[1]/out.mse1$bws[2],3)
table.HS[5,4] = paste("[",round(out.cer7$ci[3,1],2)," , ",round(out.cer7$ci[3,2],2),"]",sep="")
table.HS[5,5] = paste("[",round(out.cer8$ci[3,1],2)," , ",round(out.cer8$ci[3,2],2),"]",sep="")
table.HS[5,6] = paste("[",round(out.cer9$ci[3,1],2)," , ",round(out.cer9$ci[3,2],2),"]",sep="")

rownames(table.HS) = c("$\\widehat{h}_{\\tt MSE}$", "$\\widetilde{h}_{\\tt MSE}$",  "$\\widehat{h}_{\\tt RBC}$", "$\\widehat{h}^{\\tt rot}_{\\tt RBC}$", "$\\widetilde{h}^{\\tt rot}_{\\tt RBC}$" )
colnames(table.HS) = c("","$\\widehat{h}$", "$\\widehat{\\rho}$", "$\\widehat{\\rho}$", "$\\rho^*$", "$\\rho=1$" )

table1_tex = latex(table.HS, file = paste("CCF_Table3",".txt",sep=""),
                   landscape=FALSE, outer.size='scriptsize',  center='none', title='', table.env=FALSE,
                   n.rgroup=c(2,3),
                   col.just=rep('c',6), n.cgroup=c(1,2,3), cgroup=c("Point Estimate","Bandwidth", "RBC Confidence Intervals")) 







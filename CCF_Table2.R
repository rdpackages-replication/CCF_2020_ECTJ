################################################################################################################################################################
### CALONICO-CATTANEO-FARRELL (2019): Optimal Bandwidth Choice for Robust Bias Corrected Inference in Regression-Discontinuity Designs.
### Table 2
################################################################################################################################################################
rm(list=ls(all=TRUE))
setwd("...")
library(Hmisc)
library(Rcpp)
source("functions/CCF_rdrobust.R")
source("functions/CCF_rdbwselect.R")
source("functions/CCF_functions.R")
sourceCpp("functions/CCF_functions.cpp")

n = 500
p = 1; q = 2; deriv = 0
kernel = "tri"
sim = 5000
level=5
quant = qnorm((100-level/2)/100)
vce="hc3"
c = 0
rho.star = 0.850
w = 0.5
g.seq = seq(1/(5+2*p), 1/(3+p), length.out=10)
nseq = 10


g0.l.pop = g0.l(c);g0.r.pop = g0.r(c)
g1.l.pop = g1.l(c);g1.r.pop = g1.r(c)
g2.l.pop = g2.l(c);g2.r.pop = g2.r(c)
g3.l.pop = g3.l(c);g3.r.pop = g3.r(c)
g4.l.pop = g4.l(c);g4.r.pop = g4.r(c)
tau.pop = g0.r.pop - g0.l.pop + (gamma.r*gz0.r(c) - gamma.l*gz0.l(c))

f0.pop = fx(c)
s2.pop = sigma.y^2
C1.h = C1.fun(p,   v=0, K=k.fun);  C2.h = C2.fun(p,   v=0, K=k.fun)
C1.b = C1.fun(q,   v=2, K=k.fun);  C2.b = C2.fun(q,   v=2, K=k.fun)
C1.q = C1.fun(q+1, v=3, K=k.fun);  C2.q = C2.fun(q+1, v=3, K=k.fun)

V.IK.pop   = (s2.pop+s2.pop)/f0.pop
Vm0.IK.pop = C2.h*V.IK.pop; Vm2.IK.pop = C2.b*V.IK.pop; Vm3.IK.pop = C2.q*V.IK.pop

B.b.pop = (g3.r.pop-(-1)^(p+q+2)*g3.l.pop)/factorial(3)
D.b.pop = D.fun(q, v=2, B=B.b.pop, K=k.fun)
N.b.pop = N.fun(q, v=2, V=Vm2.IK.pop)
b.mse.pop = h.fun(N.b.pop, D.b.pop, R=0, q, v=2, n)

B.h.pop = (g2.r.pop-(-1)^(p+1)*g2.l.pop)/factorial(2);
D.h.pop = D.fun(p, v=0, B=B.h.pop, K=k.fun);
N.h.pop = N.fun(p, v=0, V=Vm0.IK.pop);

h.mse.pop = h.fun(N.h.pop, D.h.pop,    R=0, p, v=0, n)
h.cer.rot.pop = h.mse.pop*n^(-(p/((3+p)*(3+2*p))))
h.trd.pop = (h.cer.pop*n^(1/(p+3)))*n^-(quantile(g.seq,w))

nbws = 10
h.seq = b.seq = tau.us = tau.bc = se.us = se.rb = matrix(NA,sim,3*nbws)

set.seed(666)

showwhen = 1; showevery=100
for (i in 1:sim) {
  if (i==showwhen) {cat(paste("\nSimulations Completed:",i-1,"of",sim,"- N:",n,"-", Sys.time())); showwhen=showwhen+showevery}
  
  X = x.gen(n) 
  
  u.y  = u.y.rnd(n)
  u.z  = u.z.rnd(n)
  
  z.l = Z.l.fun(X,u.z)
  z.r = Z.r.fun(X,u.z)
  
  y.l = Y.l.fun(X,u.y) + gamma.l*z.l
  y.r = Y.r.fun(X,u.y) + gamma.r*z.r
  
  y = c(y.l[X<c], y.r[X>=c])
  x = c(  X[X<c],   X[X>c])
  z = c(z.l[X<c], z.r[X>=c])
  
  rdbw.out  = rdbwselect(y=y, x=x, c=c, kernel=kernel, vce=vce, scaleregul=1)
  h.mse.hat = rdbw.out$bws[1,1]
  b.mse.hat = rdbw.out$bws[1,3]
  
  eN.l <- sum(x-c>=-h.mse.hat & x<c)
  eN.r <- sum(x-c<= h.mse.hat & x>=c)
  
  h.cer.rot = h.mse.hat*n^(-(p/((3+p)*(3+2*p))))
  
  h.cer.out = rdbwselect(y=y, x=x, c=c, kernel=kernel, vce=vce, bwselect="cerdpi", tau=0.5, scaleregul=1, level=level)
  h.cer.dpi = h.cer.out$bws[1]
  b.cer.dpi = h.cer.out$bws[3]
  
  q1.rbc = h.cer.out$q.terms[1]
  q2.rbc = h.cer.out$q.terms[2]
  q3.rbc = h.cer.out$q.terms[3]
  eta.bc = h.cer.out$eta.bc
  sigma.rbc  = h.cer.out$sigma.rbc
  Q2     = abs(eta.bc^2*q2.rbc)
  H.fun  = function(W) { ( ((1-W)/W) * ((1+2*deriv)/(5+2*p)) * (4*sigma.rbc^2*quant^2)/Q2) ^ (1/(6+2*p+2*deriv)) }
  H.star = H.fun(w)
  
  h.trd.hat =  H.star*n^-quantile(g.seq,w)

  rdbwcov.out  = rdbwselect(y=y, x=x, c=c, kernel=kernel, vce=vce, scaleregul=1, covs=z)
  h.mse.hat.cov = rdbwcov.out$bws[1,1]
  b.mse.hat.cov = rdbwcov.out$bws[1,3]
  h.cer.rot.cov = h.mse.hat.cov*n^(-(p/((3+p)*(3+2*p))))
  
  h.list      = c(h.mse.pop, h.mse.hat, h.mse.hat.cov, h.cer.pop, h.cer.dpi, h.cer.rot.pop, h.cer.rot, h.cer.rot.cov, h.trd.pop, h.trd.hat)
  b.list      = c(b.mse.pop, b.mse.hat, b.mse.hat.cov, b.mse.pop, b.mse.hat, b.mse.pop,     b.mse.hat, b.mse.hat,     b.mse.pop, b.mse.hat)
  b.list.rho1 = h.list/rho.star
  b.list.rho2 = h.list
  
  h.seq[i,] = c(h.list, h.list,      h.list)
  b.seq[i,] = c(b.list, b.list.rho1, b.list.rho2)
  
  
  for (j in 1:(3*nbws)) {
    eN.l <- sum(x-c>=-h.seq[i,j] & x<c)
    eN.r <- sum(x-c<= h.seq[i,j] & x>=c)
    if (eN.l>5 & eN.r>5) {
      rd.out      = rdrobust(y=y, x=x, c=c, h=h.seq[i,j], b=b.seq[i,j], p=p, kernel=kernel, vce=vce, all = TRUE)
      tau.us[i,j] = rd.out$coef[1]
      tau.bc[i,j] = rd.out$coef[2]
      se.us[i,j]  = rd.out$se[1]
      se.rb[i,j]  = rd.out$se[3]
    }
  }
}

T.us  = (tau.us - tau.pop) / se.us
T.bc  = (tau.bc - tau.pop) / se.us
T.rb  = (tau.bc - tau.pop) / se.rb
ec.us = colMeans(1*(abs(T.us)<=quant),  na.rm = T)
ec.bc = colMeans(1*(abs(T.bc)<=quant),  na.rm = T)
ec.rb = colMeans(1*(abs(T.rb)<=quant),  na.rm = T)
il.us = colMeans(2*quant*se.us,         na.rm = T)
il.rb = colMeans(2*quant*se.rb,         na.rm = T)
h.mean= colMeans(h.seq, na.rm = T)
b.mean= colMeans(b.seq, na.rm = T)

ind.h1 = 1:nbws
ind.h2 = ind.h1 + nbws
ind.h3 = ind.h2 + nbws

table.ec = formatC(100*cbind(ec.us[ind.h1],rep(NA,10), ec.rb[ind.h1], ec.rb[ind.h2], ec.rb[ind.h3]), format = "f", digits = 1)
table.il = formatC(    cbind(il.us[ind.h1],rep(NA,10), il.rb[ind.h1], il.rb[ind.h2], il.rb[ind.h3]), format = "f", digits = 2)
table.ec[,2] = rep("",10)
table.il[,2] = rep("",10)
table.rd = cbind(table.ec, table.il)
h.list = c(h.mean[ind.h1])
b.list = c(b.mean[ind.h1])
table.bw = formatC(cbind(h.list,h.list/b.list),        format = "f", digits = 3)

out.rd = data.frame(table.bw, table.rd, row.names = NULL, stringsAsFactors = FALSE)
colnames(out.rd)= c("$h$", "$\\widehat{\\rho}$", "US", "RBC:", "$\\widehat{\\rho}$", "$\\rho^*$", "$\\rho=1$", "US", "RBC:", "$\\widehat{\\rho}$", "$\\rho^*$", "$\\rho=1$")
row.names = c("$h_{\\tt MSE}$",  "$\\widehat{h}_{\\tt MSE}$", "$\\widetilde{h}_{\\tt MSE}$" , 
              "$h_{\\tt RBC}$",  "$\\widehat{h}_{\\tt RBC}$", "$h^{\\tt rot}_{\\tt RBC}$","$\\widehat{h}^{\\tt rot}_{\\tt RBC}$", "$\\widetilde{h}^{\\tt rot}_{\\tt RBC}$",
              "$h_{\\tt TO}$", "$\\widehat{h}_{\\tt TO}$")

table_tex = latex(out.rd , file = paste("CCF_Table2",".txt",sep=""), 
                  cgroup = c("Bandwidth","Empirical Coverage", "Interval Length"), n.cgroup = c(2,5,5),   
                  n.rgroup = c(3,5,2), rowname=row.names, 
                  landscape=FALSE, center='none', col.just=rep('c',12), title='', table.env=FALSE)



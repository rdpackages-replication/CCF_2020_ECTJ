################################################################################################################################################################
### CALONICO-CATTANEO-FARRELL (2019): Optimal Bandwidth Choice for Robust Bias Corrected Inference in Regression-Discontinuity Designs.
### Auxiliary Functions
################################################################################################################################################################

c.r = c(0.26, 18.49, -54.81, 74.30, -45.02, 9.83)
c.l = c(3.71,  2.30,   3.28,  1.45,   0.23, 0.03)
cz.r = c(.49062457, 1.0600066-.45019195, 5.7429879-5.5137502, 17.142579-20.605129, 19.752974-13.317815, 7.4750031-10.955869)
cz.l = c(.49062457, 1.0600066,           5.7429879,           17.142579,           19.752974,           7.4750031)

gamma.r = gamma.l = 0
h.cer.pop = 0.1452419

sigma.y = 6.136/10
sigma.z = 0.2053

x.gen  = function(n)   { 2*rbeta(n,shape1=2,shape2=4)-1 } 
fx     = function(x) return(dbeta(x=(x+1)/2, shape1=2,shape2=4)/2)
  
g0.l = function(x) { c.l[1] + c.l[2]*x + c.l[3]*x^2 + c.l[4]*x^3   + c.l[5]*x^4     +   c.l[6]*x^5 }
g0.r = function(x) { c.r[1] + c.r[2]*x + c.r[3]*x^2 + c.r[4]*x^3   + c.r[5]*x^4     +   c.r[6]*x^5 }
g1.l = function(x) {          c.l[2] + 2*c.l[3]*x + 3*c.l[4]*x^2 + 4*c.l[5]*x^3   +   5*c.l[6]*x^4 }
g1.r = function(x) {          c.r[2] + 2*c.r[3]*x + 3*c.r[4]*x^2 + 4*c.r[5]*x^3   +   5*c.r[6]*x^4 }
g2.l = function(x) {                   2*c.l[3] + 2*3*c.l[4]*x + 3*4*c.l[5]*x^2 +   4*5*c.l[6]*x^3 }
g2.r = function(x) {                   2*c.r[3] + 2*3*c.r[4]*x + 3*4*c.r[5]*x^2 +   4*5*c.r[6]*x^3 }
g3.l = function(x) {                              2*3*c.l[4] + 2*3*4*c.l[5]*x +   3*4*5*c.l[6]*x^2 }
g3.r = function(x) {                              2*3*c.r[4] + 2*3*4*c.r[5]*x +   3*4*5*c.r[6]*x^2 }
g4.l = function(x) {                                           2*3*4*c.l[5]   + 2*3*4*5*c.l[6]*x   }
g4.r = function(x) {                                           2*3*4*c.r[5]   + 2*3*4*5*c.r[6]*x   }

gz0.l = function(x) { cz.l[1] + cz.l[2]*x + cz.l[3]*x^2 + cz.l[4]*x^3   + cz.l[5]*x^4  +   cz.l[6]*x^5 }
gz0.r = function(x) { cz.r[1] + cz.r[2]*x + cz.r[3]*x^2 + cz.r[4]*x^3   + cz.r[5]*x^4  +   cz.r[6]*x^5 }
  
u.y.rnd = function(n) { sigma.y*rnorm(n) }
u.z.rnd = function(n) { sigma.z*rnorm(n) }
Y.l.fun = function(a,u) { g0.l(a)  +  u }; Y.r.fun = function(a,u) { g0.r(a)  +  u }
Z.l.fun = function(a,u) { gz0.l(a) +  u }; Z.r.fun = function(a,u) { gz0.r(a) +  u }
  
k.fun = function(u){
    if (kernel=="epanechnikov" | kernel=="epa") {
      w = (0.75*(1-u^2)*(abs(u)<=1))
    }
    else if (kernel=="uniform" | kernel=="uni") {
      w = (0.5*(abs(u)<=1))
    }
    else {
      w = ((1-abs(u))*(abs(u)<=1))
    }
    return(w)  
}
  
m1 = function(i,j,k) integrate(function(x) x^i*x^j*k(x),0,Inf)$value
m2 = function(i,j,k) integrate(function(x) x^i*x^j*(k(x))^2,0,Inf)$value
m3 = function(i,j,k,a) integrate(function(x) x^i*(a*x)^j*k(x)*k(a*x),0,Inf)$value
GAMMA = function(p,k) {out=matrix(NA,p+1,p+1); for (i in 0:p) for (j in 0:p) out[i+1,j+1]=m1(i,j,k); return(out)}
NU    = function(p,k) {out=matrix(NA,p+1,1); for (i in 0:p) out[i+1,1]=m1(i,p+1,k); return(out)}
PSI   = function(p,k) {out=matrix(NA,p+1,p+1); for (i in 0:p) for (j in 0:p) out[i+1,j+1]=m2(i,j,k); return(out)}
B.lp  = function(p,k) {out=solve(GAMMA(p,k))%*%NU(p,k); out[1]}
C1.fun = function(p0,v,K) {
    S.inv = solve(GAMMA(p0,K)) 
    C1 = (S.inv%*%NU(p0,K))[v+1]
    return(C1)
}
C2.fun = function(p0,v,K) {
    S.inv = solve(GAMMA(p0,K)) 
    C2 = (S.inv%*%PSI(p0,K)%*%S.inv)[v+1,v+1]
    return(C2)
}
  
h.fun = function(V,B,R,p0,v,N) {(V/(N*(B+R)))^(1/(2*p0+3))}
N.fun = function(p0,v,V)       {(2*v+1)*V}
D.fun = function(p0,v,B,K)     {2*(p0+1-v)*(C1.fun(p0,v,K)*B)^2}

  
rho.optim = function(p,kernel,plot=FALSE) {
    deriv = 0
    p.comp = p+1
    ### Kernel Functions
    if (kernel=="norm") {
      k.f = function(x) (2*pi)^(-1/2)*exp(-x^2/2) 
      kernel.type="Normal"
    }
    if (kernel=="unif") {
      k.f = function(x) (abs(x)<=1)/2
      kernel.type="Uniform"
    }
    
    if (kernel=="epan") {
      k.f = function(x) 3*(1-x^2)*(abs(x)<=1)/4
      kernel.type="Epanechnikov"
    }
    if (kernel=="tria") {
      k.f = function(x) (1-abs(x))*(abs(x)<=1)
      kernel.type="Triangular"
    }
    
    k.unif = function(x) (abs(x)<=1)/2
    k.tria = function(x) (1-abs(x))*(abs(x)<=1)
    
    ### Matrix Functions
    integral.low=0

    m1 = function(i,j,k)   integrate(function(x) x^i*x^j*k(x),            integral.low,1)$value
    m2 = function(i,j,k)   integrate(function(x) x^i*x^j*(k(x))^2,        integral.low,1)$value
    m3 = function(i,j,k,a) integrate(function(x) x^i*(a*x)^j*k(x)*k(a*x), integral.low,1)$value
    
    GAMMA = function(p,k) {out=matrix(NA,p+1,p+1); for (i in 0:p) for (j in 0:p) out[i+1,j+1]=m1(i,j,k);   return(out)}
    NU    = function(p,k) {out=matrix(NA,p+1,1);   for (i in 0:p)                out[i+1,1]  =m1(i,p+1,k); return(out)}
    PSI   = function(p,k) {out=matrix(NA,p+1,p+1); for (i in 0:p) for (j in 0:p) out[i+1,j+1]=m2(i,j,k);   return(out)}
    
    B.lp = function(p,k) {out=solve(GAMMA(p,k))%*%NU(p,k); out[deriv+1]}
    
    ### Cheng-Fan-Marron Functions
    K.equi = function(p,t,k){
      v=deriv
      S.inv = solve(GAMMA(p,k));
      out=0; for (j in 0:p) out = out + S.inv[v+1,j+1]*t^j
      #return(factorial(v)*out*k(t))
      return(out*k(t))
    }
    
    ### CCT Functions
    K.equi.CCT = function(p,q,t,k,a){
      S.inv = solve(GAMMA(q,k));
      out=0; for (j in 0:q) out = out + S.inv[p+2,j+1] * (a*t)^j
      return(K.equi(p,t,k) - a^(p+2) * out * B.lp(p,k) * k(a*t))
    }
    
    K.equi.plot = function(t){
      out=K.equi(p.comp,t,k.unif)
      out[abs(t)>1]=NA
      return(out)
    }
    
    #### Optimal \rho
    dif1 = function(p,q,t,a) (abs(K.equi(p=p.comp,t=t,k=k.unif) -     K.equi.CCT(p=p,q=q,t=t,k=k.f,a=a)))^2
    integ = function(a,p,q) integrate(dif1, lower=integral.low, upper=1, p=p, q=q, a=a)$value
    
    q=p+1
    rho.star   = optimize(integ, interval=c(0,2), p=p,   q=q)$minimum

    if (plot==TRUE) {
        leg1 = paste("Uniform Equivalent Kernel",sep="")
        leg2 = paste(kernel.type," BC Equivalent Kernel",sep="")
        
        ggplot(data.frame(x = c(0, 1.5)), aes(x)) +
          stat_function(fun = K.equi.plot, geom = "line", aes(colour = leg1), size=1.2) +
          stat_function(fun = K.equi.CCT, args=list(p=p, q=q,k=k.f, a=rho.star), geom = "line", 
                        linetype = "twodash", aes(colour = leg2), size=1.2) +
          scale_colour_manual("", values = c("red","black")) +
          theme(legend.position="none", panel.background = element_blank(), 
                text = element_text(size=14), axis.line = element_line(color = "black", linetype = "solid")) +
          xlab("") + ylab("") +
          geom_segment(aes(x =  1, y = 0, xend = 1.5,  yend = 0), size=1.2) 
        
        ggsave(paste("CCF_k",kernel,"_p",p,".pdf", sep=""))  
    }
    
    return(rho.star)
  }
  
  
rdrobust_bw = function(Y, X, T, Z, C, W, c, o, nu, o_B, h_V, h_B, scale, vce, nnmatch, kernel, dups, dupsid){
  dT = dZ = dC = eC = 0
  w = rdrobust_kweight(X, c, h_V, kernel)
  dW = length(W)
  if (dW>1) {
    w = W*w
  }
  
  ind_V = w> 0; eY = Y[ind_V];eX = X[ind_V];eW = w[ind_V]
  n_V = sum(ind_V)
  D_V = eY
  R_V = matrix(NA,n_V,o+1)
  for (j in 1:(o+1)) R_V[,j] = (eX-c)^(j-1)
  invG_V = qrXXinv(R_V*sqrt(eW))
  e_v = matrix(0,(o+1),1); e_v[nu+1]=1
  s = 1
  eT=eC=eZ=NULL
  if (!is.null(T)) {
    dT = 1
    eT = T[ind_V]
    D_V = cbind(D_V,eT)
  }
  if (!is.null(Z)) {
    dZ = ncol(Z)
    eZ = Z[ind_V,,drop=FALSE]
    D_V = cbind(D_V,eZ)
    U = crossprod(R_V*eW,D_V)
    ZWD  = crossprod(eZ*eW,D_V)
    colsZ = (2+dT):max(c(2+dT+dZ-1,(2+dT)))
    UiGU =  crossprod(matrix(U[,colsZ],nrow=o+1),invG_V%*%U) 
    ZWZ = ZWD[,colsZ] - UiGU[,colsZ] 
    ZWY = ZWD[,1:(1+dT)] - UiGU[,1:(1+dT)] 
    gamma = chol2inv(chol(ZWZ))%*%ZWY
    s = c(1 , -gamma[,1])
  }
  if (!is.null(C)) {
    dC = 1
    eC =  C[ind_V] 
  }
  beta_V = invG_V%*%crossprod(R_V*eW,D_V)	
  if (is.null(Z) & !is.null(T)) {	
    tau_Y = factorial(nu)*beta_V[nu+1,1]
    tau_T = factorial(nu)*beta_V[nu+1,2]
    s = c(1/tau_T , -(tau_Y/tau_T^2))
  }
  if (!is.null(Z) & !is.null(T)) {	
    s_T = c(1 , -gamma[,2])
    tau_Y = factorial(nu)*t(s)%*%  c(beta_V[nu+1,1],beta_V[nu+1,colsZ])
    tau_T = factorial(nu)*t(s_T)%*%c(beta_V[nu+1,2],beta_V[nu+1,colsZ])
    s = c(1/tau_T , -(tau_Y/tau_T^2) , -(1/tau_T)*gamma[,1] + (tau_Y/tau_T^2)*gamma[,2])
  }	
  dups_V=dupsid_V=predicts_V=0
  
  if (vce=="nn") {
    dups_V   = dups[ind_V]
    dupsid_V = dupsid[ind_V]
  }
  
  if (vce=="hc0" | vce=="hc1" | vce=="hc2" | vce=="hc3") {
    predicts_V=R_V%*%beta_V
    if (vce=="hc2" | vce=="hc3") {
      hii=matrix(NA,n_V,1)	
      for (i in 1:n_V) {
        hii[i] = R_V[i,]%*%invG_V%*%(R_V*eW)[i,]
      }
    }
  }	
  res_V = rdrobust_res(eX, eY, eT, eZ, predicts_V, hii, vce, nnmatch, dups_V, dupsid_V, o+1)
  V_V = (invG_V%*%rdrobust_vce(dT+dZ, s, R_V*eW, res_V, eC)%*%invG_V)[nu+1,nu+1]
  v = crossprod(R_V*eW,((eX-c)/h_V)^(o+1))
  Hp = 0
  for (j in 1:(o+1)) Hp[j] = h_V^((j-1))
  BConst = (Hp*(invG_V%*%v))[nu+1]
  
  w = rdrobust_kweight(X, c, h_B, kernel)
  if (dW>1) {
    w = W*w
  }
  ind = w> 0 
  n_B = sum(ind)
  eY = Y[ind];eX = X[ind];eW = w[ind]
  D_B = eY
  R_B = matrix(NA,n_B,o_B+1)
  for (j in 1:(o_B+1)) R_B[,j] = (eX-c)^(j-1)
  invG_B = qrXXinv(R_B*sqrt(eW))
  eT=eC=eZ=NULL
  if (!is.null(T)) {
    eT = T[ind]
    D_B = cbind(D_B,eT)
  }
  if (!is.null(Z)) {
    eZ = Z[ind,,drop=FALSE]
    D_B = cbind(D_B,eZ)
  }
  if (!is.null(C)) {
    eC=C[ind]
  }	
  beta_B = invG_B%*%crossprod(R_B*eW,D_B)	
  BWreg=0
  if (scale>0) {
    e_B = matrix(0,(o_B+1),1); e_B[o+2]=1
    dups_B=dupsid_B=hii=predicts_B=0
    if (vce=="nn") {
      dups_B   = dups[ind]
      dupsid_B = dupsid[ind]
    }
    if (vce=="hc0" | vce=="hc1" | vce=="hc2" | vce=="hc3") {
      predicts_B=R_B%*%beta_B
      if (vce=="hc2" | vce=="hc3") {
        hii=matrix(NA,n_B,1)	
        for (i in 1:n_B) {
          hii[i] = R_B[i,]%*%invG_B%*%(R_B*eW)[i,]
        }
      }
    }	
    res_B = rdrobust_res(eX, eY, eT, eZ, predicts_B, hii, vce, nnmatch, dups_B, dupsid_B,o_B+1)
    V_B = (invG_B%*%rdrobust_vce(dT+dZ, s, R_B*eW, res_B, eC)%*%invG_B)[o+2,o+2]
    BWreg = 3*BConst^2*V_B
  }
  B =  sqrt(2*(o+1-nu))*BConst%*%(t(s)%*%(beta_B[o+2,]))
  V = (2*nu+1)*h_V^(2*nu+1)*V_V
  R = scale*(2*(o+1-nu))*BWreg
  rate = 1/(2*o+3)
  output = list(V=V,B=B,R=R,rate=rate)
  return(output)
}


rdrobust_kweight = function(X, c,  h,  kernel){
  u = (X-c)/h
  if (kernel=="epanechnikov" | kernel=="epa") {
    w = (0.75*(1-u^2)*(abs(u)<=1))/h
  }
  else if (kernel=="uniform" | kernel=="uni") {
    w = (0.5*(abs(u)<=1))/h
  }
  else {
    w = ((1-abs(u))*(abs(u)<=1))/h
  }
  return(w)	
}

rdrobust_res = function(X, y, T, Z, m, hii, vce, matches, dups, dupsid, d) {
  n = length(y)
  dT=dZ=0
  if (!is.null(T)) dT = 1
  if (!is.null(Z)) dZ = ncol(Z)
  res = matrix(NA,n,1+dT+dZ)  	
  
  if (vce=="nn") {
    for (pos in 1:n) {
      rpos = dups[pos] - dupsid[pos]
      lpos = dupsid[pos] - 1
      while (lpos+rpos < min(c(matches,n-1))) {
        if (pos-lpos-1 <= 0) rpos = rpos + dups[pos+rpos+1]
        else if (pos+rpos+1>n) lpos = lpos + dups[pos-lpos-1]
        else if ((X[pos]-X[pos-lpos-1]) > (X[pos+rpos+1]-X[pos])) rpos = rpos + dups[pos+rpos+1]
        else if ((X[pos]-X[pos-lpos-1]) < (X[pos+rpos+1]-X[pos])) lpos = lpos + dups[pos-lpos-1]
        else {
          rpos = rpos + dups[pos+rpos+1]
          lpos = lpos + dups[pos-lpos-1]
        }
      }
      ind_J = max(c(0,(pos-lpos))):min(c(n,(pos+rpos)))
      y_J   = sum(y[ind_J])-y[pos]
      Ji = length(ind_J)-1
      res[pos,1] = sqrt(Ji/(Ji+1))*(y[pos] - y_J/Ji)
      if (!is.null(T)) {
        T_J = sum(T[ind_J])-T[pos]
        res[pos,2] = sqrt(Ji/(Ji+1))*(T[pos] - T_J/Ji)
      }
      if (!is.null(Z)) {
        for (i in 1:dZ) {
          Z_J = sum(Z[ind_J,i])-Z[pos,i]
          res[pos,1+dT+i] = sqrt(Ji/(Ji+1))*(Z[pos,i] - Z_J/Ji)
        }
      }
    }		
  }
  else {
    if (vce=="hc0") w = 1
    else if (vce=="hc1") w = sqrt(n/(n-d))
    else if (vce=="hc2") w = sqrt(1/(1-hii))
    else                 w =      1/(1-hii)
    res[,1] = w*(y-m[,1])
    if (dT==1) res[,2] = w*(T-m[,2])
    if (dZ>0) {
      for (i in 1:dZ) {
        res[,1+dT+i] = w*(Z[,i]-m[,1+dT+i])
      }
    }
  }
  return(res)
}

qrXXinv = function(x, ...) {
  #tcrossprod(solve(qr.R(qr(x, tol = 1e-10)), tol = 1e-10))
  #tcrossprod(solve(qr.R(qr(x))))
  chol2inv(chol(crossprod(x)))
}


rdrobust_vce = function(d, s, RX, res, C) {	
  k = ncol(as.matrix(RX))
  M = matrix(0,k,k)
  n  = length(C)
  if (is.null(C)) {
    w = 1
    if (d==0){
      M  = crossprod(c(res)*RX)
    }
    else {
      for (i in 1:(1+d)) {
        SS = res[,i]*res
        for (j in 1:(1+d)) {
          M = M + crossprod(RX*(s[i]*s[j])*SS[,j],RX)
        }
      }
    }
  }
  else {	
    clusters = unique(C)
    g     = length(clusters)
    w=((n-1)/(n-k))*(g/(g-1))
    if (d==0){
      for (i in 1:g) {
        ind=C==clusters[i]
        Xi = RX[ind,,drop=FALSE]
        ri = res[ind,,drop=FALSE]
        M = M + crossprod(t(crossprod(Xi,ri)),t(crossprod(Xi,ri)))
      }
    }
    else {
      for (i in 1:g) {
        ind=C==clusters[i]
        Xi = RX[ind,,drop=FALSE]
        ri = res[ind,,drop=FALSE]
        for (l in 1:(1+d)) {	
          for (j in 1:(1+d)) {
            M = M + crossprod(t(crossprod(Xi,s[l]*ri[,l])),t(crossprod(Xi,s[j]*ri[,j])))
          }	
        }					
      }
    }
  }
  return(w*M)		
}


rdbwselect.ce.dpi = function(y, x, h, b, c, p, q, deriv, rho, kernel, vce, nnmatch, tau, quant){
  
  rho <- 1
  bwregul <- 0
  
  N = length(x)
  range = max(x)-min(x)
  
  dups <- dupsid <- hii <- predicts <- NULL
  
  if (vce=="nn") {
    order.x <- order(x)
    x <- x[order.x]
    y <- y[order.x]
  }
  
  if (!is.null(rho)){
    b <- h/rho
  } else {
    rho <- h/b
  }
  
  x.l = x[x<c];   
  x.l.min = min(x.l)
  x.l.max = max(x.l)
  range.l = abs(c-x.l.min)
  x.r = x[x>=c]
  x.r.min = min(x.r)
  x.r.max = max(x.r)
  range.r = abs(c-x.r.max)
  
  y.l = y[x<c];    y.r = y[x>=c]
  N.l = length(x.l);   N.r = length(x.r)
  x.min=min(x);  x.max=max(x)
  N = N.r + N.l
  
  X.h.l <- (x.l-c)/h;  X.h.r <- (x.r-c)/h
  X.b.l <- (x.l-c)/b;  X.b.r <- (x.r-c)/b
  K.h.l <- W.fun(X.h.l, kernel);  K.h.r <- W.fun(X.h.r, kernel)
  L.b.l <- W.fun(X.b.l, kernel);  L.b.r <- W.fun(X.b.r, kernel)
  ind.h.l <- K.h.l>0; ind.h.r <- K.h.r>0    
  ind.b.l <- L.b.l>0; ind.b.r <- L.b.r>0
  ind.l <- ind.h.l; ind.r <- ind.h.r
  if (h>b) {
    ind.l = ind.h.l
    ind.r = ind.h.r
  }
  eN.l   <- sum(ind.l); eN.r   <- sum(ind.r)
  eY.l   <-   y.l[ind.l]; eY.r   <-   y.r[ind.r]  
  eX.l   <-   x.l[ind.l]; eX.r   <-   x.r[ind.r]
  eX.h.l <- X.h.l[ind.l]; eX.h.r <- X.h.r[ind.r]  
  eX.b.l <- X.b.l[ind.l]; eX.b.r <- X.b.r[ind.r]
  eK.h.l <- K.h.l[ind.l]; eK.h.r <- K.h.r[ind.r]  
  eL.b.l <- L.b.l[ind.l]; eL.b.r <- L.b.r[ind.r]
  
  W.p.l <- eK.h.l/h; W.p.r <- eK.h.r/h  
  W.q.l <- eL.b.l/b; W.q.r <- eL.b.r/b
  R.p.2.l <- matrix(NA,eN.l,(p+3)); R.p.2.r <- matrix(NA,eN.r,(p+3))
  for (j in 1:(p+3))  {
    R.p.2.l[,j] <- eX.h.l^(j-1)
    R.p.2.r[,j] <- eX.h.r^(j-1)
  }
  R.p.1.l <- R.p.2.l[,1:(p+2)]; R.p.1.r <- R.p.2.r[,1:(p+2)]
  R.p.l   <- R.p.2.l[,1:(p+1)]; R.p.r   <- R.p.2.r[,1:(p+1)]
  R.q.l   <- matrix(NA,eN.l,(q+1)); R.q.r   <- matrix(NA,eN.r,(q+1))
  for (j in 1:(q+1))  {
    R.q.l[,j] <- eX.b.l^(j-1)
    R.q.r[,j] <- eX.b.r^(j-1)
  }
  
  L.p.1.l <- crossprod(R.p.l*W.p.l, eX.h.l^(p+1))/eN.l
  L.p.2.l <- crossprod(R.p.l*W.p.l, eX.h.l^(p+2))/eN.l
  L.p.3.l <- crossprod(R.p.l*W.p.l, eX.h.l^(p+3))/eN.l
  L.q.1.l <- crossprod(R.q.l*W.q.l, eX.b.l^(q+1))/eN.l
  L.q.2.l <- crossprod(R.q.l*W.q.l, eX.b.l^(q+2))/eN.l
  L.q.3.l <- crossprod(R.q.l*W.q.l, eX.b.l^(q+3))/eN.l
  invG.p.l <- eN.l*qrXXinv(R.p.l*sqrt(W.p.l))
  invG.q.l <- eN.l*qrXXinv(R.q.l*sqrt(W.q.l))
  
  L.p.1.r <- crossprod(R.p.r*W.p.r, eX.h.r^(p+1))/eN.r
  L.p.2.r <- crossprod(R.p.r*W.p.r, eX.h.r^(p+2))/eN.r
  L.p.3.r <- crossprod(R.p.r*W.p.r, eX.h.r^(p+3))/eN.r
  L.q.1.r <- crossprod(R.q.r*W.q.r, eX.b.r^(q+1))/eN.r
  L.q.2.r <- crossprod(R.q.r*W.q.r, eX.b.r^(q+2))/eN.r
  L.q.3.r <- crossprod(R.q.r*W.q.r, eX.b.r^(q+3))/eN.r
  invG.p.r <- eN.r*qrXXinv(R.p.r*sqrt(W.p.r))
  invG.q.r <- eN.r*qrXXinv(R.q.r*sqrt(W.q.r))
  
  edups.l = edupsid.l = matrix(0,eN.l,1)
  edups.r = edupsid.r = matrix(0,eN.r,1)
  
  if (vce=="nn") {
    for (i in 1:eN.l) edups.l[i]=sum(eX.l==eX.l[i])
    for (i in 1:eN.r) edups.r[i]=sum(eX.r==eX.r[i])
    i=1
    while (i<=eN.l) {
      edupsid.l[i:(i+edups.l[i]-1)] = 1:edups.l[i]
      i = i+edups.l[i]
    }
    i=1
    while (i<=eN.r) {
      edupsid.r[i:(i+edups.r[i]-1)]=1:edups.r[i]
      i=i+edups.r[i]
    }
  }          
  
  hii.l=hii.r=predicts.p.l=predicts.p.r=predicts.q.l=predicts.q.r=0
  if (vce=="hc0" | vce=="hc1" | vce=="hc2" | vce=="hc3") {
    
    H.q <- 0
    for (j in 1:(q+1)) H.q[j] <- b^(-(j-1))
    
    beta.q.l <- H.q*invG.q.l%*%crossprod(R.q.l*W.q.l,eY.l)/eN.l
    beta.q.r <- H.q*invG.q.r%*%crossprod(R.q.r*W.q.r,eY.r)/eN.r
    
    r.q.l = matrix(NA,eN.l,q+1)
    r.q.r = matrix(NA,eN.r,q+1)
    
    predicts.q.l = predicts.q.r = 0
    for (j in 1:(q+1))  {
      r.q.l[,j] <- (eX.l-c)^(j-1)
      r.q.r[,j] <- (eX.r-c)^(j-1)
    }
    
    for (j in 1:eN.l) predicts.q.l[j] <- r.q.l[j,]%*%beta.q.l
    for (j in 1:eN.r) predicts.q.r[j] <- r.q.r[j,]%*%beta.q.r
    
    if (vce=="hc2" | vce=="hc3") {
      hii.l <- matrix(NA,eN.l,1)
      hii.r <- matrix(NA,eN.r,1)
      for (j in 1:eN.l) hii.l[j] <- (R.p.l[j,]%*%invG.p.l%*%(R.p.l*W.p.l)[j,])/eN.l
      for (j in 1:eN.r) hii.r[j] <- (R.p.r[j,]%*%invG.p.r%*%(R.p.r*W.p.r)[j,])/eN.r
    }
  }
  
  res.q.l <- rdrobust_res(eX.l, eY.l, NULL, NULL, as.matrix(predicts.q.l), hii.l, vce, nnmatch, edups.l, edupsid.l, q+1)
  res.q.r <- rdrobust_res(eX.r, eY.r, NULL, NULL, as.matrix(predicts.q.r), hii.r, vce, nnmatch, edups.r, edupsid.r, q+1)
  
  e.p.1 <- matrix(0,(q+1),1); e.p.1[p+2]=1
  e.0   <- matrix(0,(p+1),1); e.0[1]=1
  
  ### Bias
  #k <- p+3
  #r.k <- matrix(NA,N,k+3)
  #for (j in 1:(k+3))  r.k[,j] <- x^(j-1)
  #gamma <- lm(y~r.k-1)
  #m.p.3 <- gamma$coeff[p+4]*factorial(p+3) + gamma$coeff[p+5]*factorial(p+4)*c + gamma$coeff[p+6]*factorial(p+5)*c^2/2
  #r.k <- r.k[,1:(k+2)]
  #gamma <- lm(y~r.k[,1:(k+2)]-1)
  #m.p.2 <- gamma$coeff[p+3]*factorial(p+2) + gamma$coeff[p+4]*factorial(p+3)*c + gamma$coeff[p+5]*factorial(p+4)*c^2/2
  #m.p.2.l =m.p.2.r=m.p.2
  #m.p.3.l =m.p.3.r=m.p.3
  
  #if (is.na(m.p.3)) m.p.3 <- lprobust(y=y, x=x, h=range, c=c, p=p+4, deriv=(p+3), kernel=kernel, vce=vce)$Estimate[5]
  #if (is.na(m.p.2)) m.p.2.l <- rdrobust(y=y.l, x=x.l, h=range, c=c, p=p+4, deriv=(p+2), kernel=kernel, vce=vce)$Estimate[5]
  m.p.2   <- rdrobust(y=y, x=x, c=c, p=p+3, deriv=(p+2), kernel=kernel, vce=vce)
  m.p.2.r <- factorial(p+2)*m.p.2$beta_p_r[p+3]
  m.p.2.l <- factorial(p+2)*m.p.2$beta_p_l[p+3]
  
  #m.p.2.l = g2.mnus(c)
  #m.p.2.r = g2.plus(c)
  
  #eta.bc.r  <- (m.p.2.r/factorial(p+2))*(t(e.0)%*%invG.p.r%*%L.p.2.r/h - rho^(-2)*b^(q-p-1)*(t(e.0)%*%invG.p.r%*%L.p.1.r%*%t(e.p.1)%*%invG.q.r%*%L.q.1.r/b))
  #eta.bc.l  <- (m.p.2.l/factorial(p+2))*(t(e.0)%*%invG.p.l%*%L.p.2.l/h - rho^(-2)*b^(q-p-1)*(t(e.0)%*%invG.p.l%*%L.p.1.l%*%t(e.p.1)%*%invG.q.l%*%L.q.1.l/b))
  
  eta.bc.r  <- (m.p.2.r/factorial(p+2))*(t(e.0)%*%invG.p.r%*%L.p.2.r - rho^(-1)*(t(e.0)%*%invG.p.r%*%L.p.1.r%*%t(e.p.1)%*%invG.q.r%*%L.q.1.r))
  eta.bc.l  <- (m.p.2.l/factorial(p+2))*(t(e.0)%*%invG.p.l%*%L.p.2.l - rho^(-1)*(t(e.0)%*%invG.p.l%*%L.p.1.l%*%t(e.p.1)%*%invG.q.l%*%L.q.1.l))
  eta.bc   <- c(eta.bc.r - eta.bc.l)
  
  dquant <- dnorm(quant)
  zeros.l = rep(0,eN.l)
  zeros.r = rep(0,eN.r)
  
  eK.h.l = c(eK.h.l, zeros.r)
  eL.b.l = c(eL.b.l, zeros.r)
  eK.h.r = c(zeros.l, eK.h.r)
  eL.b.r = c(zeros.l, eL.b.r)
  
  eX = c(eX.l, eX.r)
  eY = c(eY.l, eY.r)
  
  res.q.l = c(res.q.l, zeros.r)
  res.q.r = c(zeros.l, res.q.r)
  
  #res.q.l = c(resid.pob.l[ind.l], zeros.r)
  #res.q.r = c(zeros.l, resid.pob.r[ind.r])
  
  q.terms <- rdbwce(y_l = eY, y_r=eY, x_l = eX, x_r=eX,  K_l=eK.h.l,  K_r=eK.h.r, L_l=eL.b.l,  L_r=eL.b.r, res_l= res.q.l, res_r= res.q.r, c=c, p=p, q=q, h=h, b=b, deriv=deriv, fact=factorial(deriv), z=quant)
  
  q1.rbc  <- 2*c(q.terms$q1rbc)
  q2.rbc  <- 2*c(q.terms$q2rbc)
  q3.rbc  <- 2*c(q.terms$q3rbc)
  s2.rbc  <- c(q.terms$s2rbc)
  
  #eta.bc.1.l <- (t(e.0)%*%invG.p.l)%*%( L.p.2.l - rho^(-1)*L.p.1.l%*%t(e.p.1)%*%invG.q.l%*%L.q.1.l)*(m.p.2.l/factorial(p+2))
  #eta.bc.1.r <- (t(e.0)%*%invG.p.r)%*%( L.p.2.r - rho^(-1)*L.p.1.r%*%t(e.p.1)%*%invG.q.r%*%L.q.1.r)*(m.p.2.r/factorial(p+2))
  #eta.bc.2.l <- (t(e.0)%*%invG.p.l)%*%( L.p.3.l - rho^(-2)*L.p.1.l%*%t(e.p.1)%*%invG.q.l%*%L.q.2.l)*(m.p.3.l/factorial(p+3))
  #eta.bc.2.r <- (t(e.0)%*%invG.p.r)%*%( L.p.3.r - rho^(-2)*L.p.1.r%*%t(e.p.1)%*%invG.q.r%*%L.q.2.r)*(m.p.3.r/factorial(p+3))
  #eta.bc.1 = eta.bc.1.r - eta.bc.1.l
  #eta.bc.2 = eta.bc.2.r - eta.bc.2.l
  #Reg      <- 3*(t(e.0)%*%invG.p)%*%(L.p.2/h)^2*V.reg
  
  
  CE.fun <- function(H) {H^(-1)*q1.rbc + H^(5+2*p)*(eta.bc^2)*q2.rbc + H^(p+2)*(eta.bc)*q3.rbc}
  H.bc   <-function(H) {check.f(CE.fun(H), tau)}
  h.bc     <- optimize(H.bc , interval=c(.Machine$double.eps, range))
  h.ce.dpi <- as.numeric(h.bc$minimum*N^(-1/(p+3)))
  
  #sigma = sqrt(s2.rbc*N*h)
  sigma = sqrt(s2.rbc)
  
  out <- list(h=h.ce.dpi, sigma=sigma, eta.bc=eta.bc, q.terms=c(q1.rbc,q2.rbc,q3.rbc))
  return(out)
}


rdbwselect.ce.pop = function(y, x, u, h, b, c, p, q, deriv, m2, rho, kernel, tau, quant){
  
  rho <- 1
  bwregul <- 0
  
  N = length(x)
  range = max(x)-min(x)
  
  if (!is.null(rho)){
    b <- h/rho
  } else {
    rho <- h/b
  }
  
  x.l = x[x<c];   
  x.l.min = min(x.l)
  x.l.max = max(x.l)
  range.l = abs(c-x.l.min)
  x.r = x[x>=c]
  x.r.min = min(x.r)
  x.r.max = max(x.r)
  range.r = abs(c-x.r.max)
  
  y.l = y[x<c];    y.r = y[x>=c]
  u.l = u[x<c];    u.r = u[x>=c]
  
  N.l = length(x.l);   N.r = length(x.r)
  x.min=min(x);  x.max=max(x)
  N = N.r + N.l
  
  X.h.l <- (x.l-c)/h;  X.h.r <- (x.r-c)/h
  X.b.l <- (x.l-c)/b;  X.b.r <- (x.r-c)/b
  K.h.l <- W.fun(X.h.l, kernel);  K.h.r <- W.fun(X.h.r, kernel)
  L.b.l <- W.fun(X.b.l, kernel);  L.b.r <- W.fun(X.b.r, kernel)
  ind.h.l <- K.h.l>0; ind.h.r <- K.h.r>0    
  ind.b.l <- L.b.l>0; ind.b.r <- L.b.r>0
  ind.l <- ind.h.l; ind.r <- ind.h.r
  if (h>b) {
    ind.l = ind.h.l
    ind.r = ind.h.r
  }
  eN.l   <- sum(ind.l); eN.r   <- sum(ind.r)
  eY.l   <-   y.l[ind.l]; eY.r   <-   y.r[ind.r]  
  eX.l   <-   x.l[ind.l]; eX.r   <-   x.r[ind.r]
  eX.h.l <- X.h.l[ind.l]; eX.h.r <- X.h.r[ind.r]  
  eX.b.l <- X.b.l[ind.l]; eX.b.r <- X.b.r[ind.r]
  eK.h.l <- K.h.l[ind.l]; eK.h.r <- K.h.r[ind.r]  
  eL.b.l <- L.b.l[ind.l]; eL.b.r <- L.b.r[ind.r]
  
  W.p.l <- eK.h.l/h; W.p.r <- eK.h.r/h  
  W.q.l <- eL.b.l/b; W.q.r <- eL.b.r/b
  R.p.2.l <- matrix(NA,eN.l,(p+3)); R.p.2.r <- matrix(NA,eN.r,(p+3))
  for (j in 1:(p+3))  {
    R.p.2.l[,j] <- eX.h.l^(j-1)
    R.p.2.r[,j] <- eX.h.r^(j-1)
  }
  R.p.1.l <- R.p.2.l[,1:(p+2)]; R.p.1.r <- R.p.2.r[,1:(p+2)]
  R.p.l   <- R.p.2.l[,1:(p+1)]; R.p.r <- R.p.2.r[,1:(p+1)]
  R.q.l   <- matrix(NA,eN.l,(q+1)); R.q.r   <- matrix(NA,eN.r,(q+1))
  for (j in 1:(q+1))  {
    R.q.l[,j] <- eX.b.l^(j-1)
    R.q.r[,j] <- eX.b.r^(j-1)
  }
  
  L.p.1.l <- crossprod(R.p.l*W.p.l, eX.h.l^(p+1))/eN.l
  L.p.2.l <- crossprod(R.p.l*W.p.l, eX.h.l^(p+2))/eN.l
  L.p.3.l <- crossprod(R.p.l*W.p.l, eX.h.l^(p+3))/eN.l
  L.q.1.l <- crossprod(R.q.l*W.q.l, eX.b.l^(q+1))/eN.l
  L.q.2.l <- crossprod(R.q.l*W.q.l, eX.b.l^(q+2))/eN.l
  L.q.3.l <- crossprod(R.q.l*W.q.l, eX.b.l^(q+3))/eN.l
  invG.p.l <- eN.l*qrXXinv(R.p.l*sqrt(W.p.l))
  invG.q.l <- eN.l*qrXXinv(R.q.l*sqrt(W.q.l))
  
  L.p.1.r <- crossprod(R.p.r*W.p.r, eX.h.r^(p+1))/eN.r
  L.p.2.r <- crossprod(R.p.r*W.p.r, eX.h.r^(p+2))/eN.r
  L.p.3.r <- crossprod(R.p.r*W.p.r, eX.h.r^(p+3))/eN.r
  L.q.1.r <- crossprod(R.q.r*W.q.r, eX.b.r^(q+1))/eN.r
  L.q.2.r <- crossprod(R.q.r*W.q.r, eX.b.r^(q+2))/eN.r
  L.q.3.r <- crossprod(R.q.r*W.q.r, eX.b.r^(q+3))/eN.r
  invG.p.r <- eN.r*qrXXinv(R.p.r*sqrt(W.p.r))
  invG.q.r <- eN.r*qrXXinv(R.q.r*sqrt(W.q.r))
  
  e.p.1 <- matrix(0,(q+1),1); e.p.1[p+2]=1
  e.0   <- matrix(0,(p+1),1); e.0[1]=1
  
  ### Bias
  m.p.2   <- rdrobust(y=y, x=x, c=c, p=p+3, deriv=(p+2), kernel=kernel, vce=vce)
  m.p.2.r <- factorial(p+2)*m.p.2$beta_p_r[p+3]
  m.p.2.l <- factorial(p+2)*m.p.2$beta_p_l[p+3]
  
  #m.p.2.l = m2[1]
  #m.p.2.r = m2[2]
  
  #eta.bc.r  <- (m.p.2.r/factorial(p+2))*(t(e.0)%*%invG.p.r%*%L.p.2.r/h - rho^(-2)*b^(q-p-1)*(t(e.0)%*%invG.p.r%*%L.p.1.r%*%t(e.p.1)%*%invG.q.r%*%L.q.1.r/b))
  #eta.bc.l  <- (m.p.2.l/factorial(p+2))*(t(e.0)%*%invG.p.l%*%L.p.2.l/h - rho^(-2)*b^(q-p-1)*(t(e.0)%*%invG.p.l%*%L.p.1.l%*%t(e.p.1)%*%invG.q.l%*%L.q.1.l/b))
  
  eta.bc.r  <- (m.p.2.r/factorial(p+2))*(t(e.0)%*%invG.p.r%*%L.p.2.r - rho^(-1)*(t(e.0)%*%invG.p.r%*%L.p.1.r%*%t(e.p.1)%*%invG.q.r%*%L.q.1.r))
  eta.bc.l  <- (m.p.2.l/factorial(p+2))*(t(e.0)%*%invG.p.l%*%L.p.2.l - rho^(-1)*(t(e.0)%*%invG.p.l%*%L.p.1.l%*%t(e.p.1)%*%invG.q.l%*%L.q.1.l))
  
  eta.bc   <- c(eta.bc.r - eta.bc.l)
  
  dquant <- dnorm(quant)
  zeros.l = rep(0,eN.l)
  zeros.r = rep(0,eN.r)
  
  eK.h.l = c(eK.h.l, zeros.r)
  eL.b.l = c(eL.b.l, zeros.r)
  eK.h.r = c(zeros.l, eK.h.r)
  eL.b.r = c(zeros.l, eL.b.r)
  
  eX = c(eX.l, eX.r)
  eY = c(eY.l, eY.r)
  
  #res.q.l = c(res.q.l, zeros.r)
  #res.q.r = c(zeros.l, res.q.r)
  
  res.q.l = c(u.l[ind.l], zeros.r)
  res.q.r = c(zeros.l, u.r[ind.r])
  
  q.terms <- rdbwce(y_l = eY, y_r=eY, x_l = eX, x_r=eX,  K_l=eK.h.l,  K_r=eK.h.r, L_l=eL.b.l,  L_r=eL.b.r, res_l= res.q.l, res_r= res.q.r, c=c, p=p, q=q, h=h, b=b, deriv=deriv, fact=factorial(deriv), z=quant)
  
  q1.rbc  <- 2*c(q.terms$q1rbc)
  q2.rbc  <- 2*c(q.terms$q2rbc)
  q3.rbc  <- 2*c(q.terms$q3rbc)
  s2.rbc  <- c(q.terms$s2rbc)
  
  CE.fun <- function(H) {H^(-1)*q1.rbc + H^(5+2*p)*(eta.bc^2)*q2.rbc + H^(p+2)*(eta.bc)*q3.rbc}
  H.bc   <-function(H) {check.f(CE.fun(H), tau)}
  h.bc     <- optimize(H.bc , interval=c(.Machine$double.eps, range))
  h.ce.pop <- as.numeric(h.bc$minimum*N^(-1/(p+3)))
  
  out <- list(h=c(h.ce.pop))
  return(out)
}

W.fun = function(u,kernel){
  if (kernel=="epa"|kernel=="epanechnikov") w = 0.75*(1-u^2)*(abs(u)<=1)
  if (kernel=="uni"|kernel=="uniform")      w =          0.5*(abs(u)<=1)
  if (kernel=="tri"|kernel=="triangular")   w =   (1-abs(u))*(abs(u)<=1)
  return(w)
}

check.f = function(u,tau) {
  (tau - 1*(u<0))*u
}


qrreg = function(x,y,w,s2=0,var.comp=TRUE, ...) {
  M.X = sqrt(w)*x
  X.M.X_inv = qrXXinv(M.X) 
  X.M.Y = crossprod(M.X,sqrt(w)*y)
  beta.hat = X.M.X_inv%*%X.M.Y
  Psi.hat=Sigma.hat=0
  if (var.comp==TRUE) {
    Psi.hat = crossprod((w*s2*w)*x,x)
    Sigma.hat = crossprod(Psi.hat%*%X.M.X_inv,X.M.X_inv)
  }
  output = list(X.M.X_inv=X.M.X_inv, X.M.Y=X.M.Y, beta.hat=beta.hat, Psi.hat=Psi.hat, Sigma.hat=Sigma.hat)
  return(output)
}





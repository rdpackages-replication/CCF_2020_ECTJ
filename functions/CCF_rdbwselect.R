rdbwselect = function(y, x, c=NULL, fuzzy = NULL, deriv=NULL, p=NULL, q=NULL, covs = NULL,  
                      kernel="tri", weights=NULL, bwselect="mserd", vce="nn", cluster = NULL, 
                      nnmatch=3,  scaleregul=1, sharpbw=FALSE,  all=NULL, subset = NULL, tau=0.5, level=5){
  
  quant = qnorm((100-level/2)/100)
  
  if (!is.null(subset)) {
    x <- x[subset]
    y <- y[subset]
  }
  
  if (is.null(all)) all<-FALSE
  if (is.null(c)) c <- 0
    
  # p
  if (length(p) == 0) {
    flag_no_p <- TRUE
    p <- 1
  } else if ((length(p) != 1) | !(p[1]%in%0:20)) {
    stop("Polynomial order p incorrectly specified.\n")
  } else {
    flag_no_p <- FALSE
  }
  
  # q
  if (length(q) == 0) {
    flag_no_q <- TRUE
    q <- p + 1
  } else if ((length(q) > 1) | !(q[1]%in%c(0:20)) | (q[1]<p)) {
    stop("Polynomial order (for bias correction) q incorrectly specified.\n")
  } else {
    flag_no_q <- FALSE
  }
  
  # deriv
  if (length(deriv) == 0) {
    flag_no_deriv <- TRUE
    deriv <- 0
  } else if ((length(deriv) > 1) | !(deriv[1]%in%c(0:20)) | (deriv[1]>p)) {
    stop("Derivative order incorrectly specified.\n")
  } else {
    flag_no_deriv <- FALSE
  }
  
  na.ok <- complete.cases(x) & complete.cases(y)
  
  if (!is.null(cluster)){
    if (!is.null(subset))  cluster <- cluster[subset]
    na.ok <- na.ok & complete.cases(cluster)
  } 
  
  if (!is.null(covs)){
    if (!is.null(subset))  covs <- subset(covs,subset)
    na.ok <- na.ok & complete.cases(covs)
  } 
  
  if (!is.null(fuzzy)){
    if (!is.null(subset)) fuzzy <- fuzzy[subset]
    na.ok <- na.ok & complete.cases(fuzzy)
  } 
  
  if (!is.null(weights)){
    if (!is.null(subset)) weights <- weights[subset]
    na.ok <- na.ok & complete.cases(weights)
  } 
  
  x = as.matrix(x[na.ok])
  y = as.matrix(y[na.ok])
  if (!is.null(covs))    covs    = as.matrix(covs)[na.ok, , drop = FALSE]
  if (!is.null(fuzzy))   fuzzy   = as.matrix(fuzzy[na.ok])
  if (!is.null(cluster)) cluster = as.matrix(cluster[na.ok])
  if (!is.null(weights)) weights = as.matrix(weights[na.ok])
  
  
  if (vce=="nn") {
    order_x = order(x)
    x = x[order_x,,drop=FALSE]
    y = y[order_x,,drop=FALSE]
    if (!is.null(covs))    covs    =  as.matrix(covs)[order_x,,drop=FALSE]
    if (!is.null(fuzzy))   fuzzy   =   fuzzy[order_x,,drop=FALSE]
    if (!is.null(cluster)) cluster = cluster[order_x,,drop=FALSE]
    if (!is.null(weights)) weights = weights[order_x,,drop=FALSE]
  }
  
  ### reescaling
  y_sd = sd(y)
	y = y/y_sd
	x_sd = sd(x)
  x = x/x_sd
	c_orig = c
  c = c/x_sd
	
  #x_sd = sd(x)
  #x_iq = quantile(x,.75, type=6) - quantile(x,.25, type=6)
  x_iq = quantile(x,.75,type=2) - quantile(x,.25,type=2)
  
  X_l = x[x<c];   
  x_l_min = min(X_l)
  x_l_max = max(X_l)
  range_l = abs(c-x_l_min)
  X_r = x[x>=c]
  x_r_min = min(X_r)
  x_r_max = max(X_r)
  range_r = abs(c-x_r_max)

  Y_l = y[x<c];    Y_r = y[x>=c]
  N_l = length(X_l);   N_r = length(X_r)
  x_min=min(x);  x_max=max(x)
  N = N_r + N_l

    exit=0
    #################  ERRORS
    if (kernel!="uni" & kernel!="uniform" & kernel!="tri" & kernel!="triangular" & kernel!="epa" & kernel!="epanechnikov" & kernel!="" ){
      print("kernel incorrectly specified")
      exit = 1
    }
    
    if  (bwselect!="cerdpi" & bwselect!="mserd" & bwselect!="msetwo" & bwselect!="msesum" & bwselect!="msecomb1" & bwselect!="msecomb2"  & bwselect!="cerrd" & bwselect!="certwo" & bwselect!="cersum" & bwselect!="cercomb1" & bwselect!="cercomb2" & bwselect!=""){
      print("bwselect incorrectly specified")  
      exit = 1
    }
    
    if (bwselect=="CCT" | bwselect=="IK" | bwselect=="CV" | bwselect=="cct" | bwselect=="ik" | bwselect=="cv"){
      print("bwselect options IK, CCT and CV have been depricated. Please see help for new options")  
      exit = 1
    }
    
    if (vce!="nn" & vce!="" & vce!="hc1" & vce!="hc2" & vce!="hc3" & vce!="hc0"){ 
      print("vce incorrectly specified")
      exit = 1
    }

    if (c<=x_min | c>=x_max){
      print("c should be set within the range of x")
      exit = 1
    }
    
    if (p<0 | q<0 | deriv<0 | nnmatch<=0 ){
      print("p, q, deriv and matches should be positive integers")
      exit = 1
    }
    
    if (p>=q){
      print("q should be set higher than p")
      exit = 1
    }
    
    if (deriv>p){
      print("deriv can only be equal or lower p")
      exit = 1
    }
    
    p_round = round(p)/p;    q_round = round(q)/q;    d_round = round(deriv+1)/(deriv+1);    m_round = round(nnmatch)/nnmatch
        
    if ((p_round!=1 &p>0) | (q_round!=1&q>0) | d_round!=1 | m_round!=1 ){
      print("p,q,deriv and matches should be integer numbers")
      exit = 1
    }
    
    if (exit>0) stop()
    
  
  if (kernel=="epanechnikov" | kernel=="epa") {
    kernel_type = "Epanechnikov"
    C_c=2.34
  }  else if (kernel=="uniform" | kernel=="uni") {
    kernel_type = "Uniform"
    C_c=1.843
  }   else  {
    kernel_type = "Triangular"
    C_c=2.576
  }
  
    vce_type = "NN"
    if (vce=="hc0")     		vce_type = "HC0"
    if (vce=="hc1")      	  vce_type = "HC1"
    if (vce=="hc2")      	  vce_type = "HC2"
    if (vce=="hc3")      	  vce_type = "HC3"
    if (vce=="cluster")  	  vce_type = "Cluster"
    if (vce=="nncluster") 	vce_type = "NNcluster"

    
  #***********************************************************************
  dZ=Z_l=Z_r=T_l=T_r=C_l=C_r=Cind_l=Cind_r=g_l=g_r=NULL

  dups_l = dupsid_l = matrix(0,N_l,1)
  dups_r = dupsid_r = matrix(0,N_r,1)
  
  if (vce=="nn") {
    for (i in 1:N_l) {
      dups_l[i]=sum(X_l==X_l[i])
    }
    for (i in 1:N_r) {
      dups_r[i]=sum(X_r==X_r[i])
    }
    i=1
    while (i<=N_l) {
      dupsid_l[i:(i+dups_l[i]-1)] = 1:dups_l[i]
      i = i+dups_l[i]
    }
    i=1
    while (i<=N_r) {
      dupsid_r[i:(i+dups_r[i]-1)]=1:dups_r[i]
      i=i+dups_r[i]
    }
  } 
  
  if (!is.null(covs)) {
    dZ = ncol(covs)
    Z_l  = covs[x<c,,drop=FALSE];  Z_r  = covs[x>=c,,drop=FALSE]
  }
  perf_comp=FALSE
  if (!is.null(fuzzy)) {
    T_l  = fuzzy[x<c,,drop=FALSE];  T_r  = fuzzy[x>=c,,drop=FALSE]; 
    if (var(T_l)==0 | var(T_r)==0) perf_comp=TRUE
    if (perf_comp==TRUE | sharpbw==TRUE) {
      T_l = T_r = NULL
      }
    }
   
  if (!is.null(cluster)) {
    C_l  = cluster[x<c,,drop=FALSE]; C_r= cluster[x>=c,,drop=FALSE]
    g_l = length(unique(C_l));	g_r = length(unique(C_r))
  }
  fw_l = fw_r = 0 
  if (!is.null(weights)) {
    fw_l=weights[x<c];  fw_r=weights[x>=c]
  }                                                                           
    #***********************************************************************
    c_bw = C_c*min(c(1,x_iq/1.349))*N^(-1/5)
  
    #*** Step 1: d_bw
    C_d_l = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c=c, o=q+1, nu=q+1, o_B=q+2, h_V=c_bw, h_B=range_l, 0, vce, nnmatch, kernel, dups_l, dupsid_l)
    C_d_r = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c=c, o=q+1, nu=q+1, o_B=q+2, h_V=c_bw, h_B=range_r, 0, vce, nnmatch, kernel, dups_r, dupsid_r)
  
    #if (C_d_l$V=="NaN" | C_d_l$B=="NaN" | C_d_l$R=="NaN") C_d_l = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, c=c, o=q+1, nu=q+1, o_B=q+2, h_V=c_bw, h_B=c_bw, 0, vce, nnmatch, kernel, dups_l, dupsid_l)
    #if (C_d_r$V=="NaN" | C_d_r$B=="NaN" | C_d_r$R=="NaN") C_d_r = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, c=c, o=q+1, nu=q+1, o_B=q+2, h_V=c_bw, h_B=c_bw, 0, vce, nnmatch, kernel, dups_r, dupsid_r)
    #if (C_d_l$V==. | C_d_l$B==. | C_d_l$R==.) display("Invertibility problem in the computation of preliminary bandwidth below the threshold")  
    #if (C_d_r$V==. | C_d_r$B==. | C_d_r$R==.) display("Invertibility problem in the computation of preliminary bandwidth above the threshold")  
    #if (C_d_l$V==0 | C_d_l$B==0) display("Not enough variability to compute the preliminary bandwidth below the threshold. Range defined by bandwidth: ")  
    #if (C_d_r$V==0 | C_d_r$B==0) display("Not enough variability to compute the preliminary bandwidth above the threshold. Range defined by bandwidth: ")  
 
    #*** TWO
    if  (bwselect=="msetwo" |  bwselect=="certwo" | bwselect=="msecomb2" | bwselect=="cercomb2"  | all=="TRUE")  {		
      d_bw_l = c((  C_d_l$V              /   C_d_l$B^2             )^C_d_l$rate)
      d_bw_r = c((  C_d_r$V              /   C_d_r$B^2             )^C_d_l$rate)
      C_b_l  = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c=c, o=q, nu=p+1, o_B=q+1, h_V=c_bw, h_B=d_bw_l, scaleregul, vce, nnmatch, kernel, dups_l, dupsid_l)
      b_bw_l = c((  C_b_l$V              /   (C_b_l$B^2 + scaleregul*C_b_l$R)        )^C_b_l$rate)
      C_b_r  = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c=c, o=q, nu=p+1, o_B=q+1, h_V=c_bw, h_B=d_bw_r, scaleregul, vce, nnmatch, kernel, dups_r, dupsid_r)
      b_bw_r = c((  C_b_r$V              /   (C_b_r$B^2 + scaleregul*C_b_r$R)        )^C_b_l$rate)
      C_h_l  = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c=c, o=p, nu=deriv, o_B=q, h_V=c_bw, h_B=b_bw_l, scaleregul, vce, nnmatch, kernel, dups_l, dupsid_l)
      h_bw_l = c((  C_h_l$V              /   (C_h_l$B^2 + scaleregul*C_h_l$R)         )^C_h_l$rate)
      C_h_r  = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c=c, o=p, nu=deriv, o_B=q, h_V=c_bw, h_B=b_bw_r, scaleregul, vce, nnmatch, kernel, dups_r, dupsid_r)
      h_bw_r = c((  C_h_r$V              /   (C_h_r$B^2 + scaleregul*C_h_r$R)         )^C_h_l$rate)
                  
      #if (C_b_l$V==0 | C_b_l$B==0) printf("{err}Not enough variability to compute the bias bandwidth (b) below the threshold. Range defined by bandwidth = %f\n", d_bw_l)  
      #if (C_b_r$V==0 | C_b_r$B==0) printf("{err}Not enough variability to compute the bias bandwidth (b) above the threshold. Range defined by bandwidth = %f\n", d_bw_r)  
      #if (C_h_l$V==0 | C_h_l$B==0) printf("{err}Not enough variability to compute the loc. poly. bandwidth (h) below the threshold. Range defined by bandwidth = %f\n", b_bw_l) 
      #if (C_h_r$V==0 | C_h_r$B==0) printf("{err}Not enough variability to compute the loc. poly. bandwidth (h) above the threshold. Range defined by bandwidth = %f\n", b_bw_r) 
    }
  
#  *** SUM
  if  (bwselect=="msesum" | bwselect=="cersum" |  bwselect=="msecomb1" | bwselect=="msecomb2" |  bwselect=="cercomb1" | bwselect=="cercomb2"  |  all=="TRUE")  {
    d_bw_s = c(( (C_d_l$V + C_d_r$V)  /  (C_d_r$B + C_d_l$B)^2 )^C_d_l$rate)
    C_b_l  = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c=c, o=q, nu=p+1, o_B=q+1, h_V=c_bw, h_B=d_bw_s, scaleregul, vce, nnmatch, kernel, dups_l, dupsid_l)
    C_b_r  = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c=c, o=q, nu=p+1, o_B=q+1, h_V=c_bw, h_B=d_bw_s, scaleregul, vce, nnmatch, kernel, dups_r, dupsid_r)
    b_bw_s = c(( (C_b_l$V + C_b_r$V)  /  ((C_b_r$B + C_b_l$B)^2 + scaleregul*(C_b_r$R+C_b_l$R)) )^C_b_l$rate)
    C_h_l  = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c=c, o=p, nu=deriv, o_B=q, h_V=c_bw, h_B=b_bw_s, scaleregul, vce, nnmatch, kernel, dups_l, dupsid_l)
    C_h_r  = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c=c, o=p, nu=deriv, o_B=q, h_V=c_bw, h_B=b_bw_s, scaleregul, vce, nnmatch, kernel, dups_r, dupsid_r)
    h_bw_s = c(( (C_h_l$V + C_h_r$V)  /  ((C_h_r$B + C_h_l$B)^2 + scaleregul*(C_h_r$R + C_h_l$R)) )^C_h_l$rate)

    #if (C_b_l$V==0 | C_b_l$B==0) printf("{err}Not enough variability to compute the bias bandwidth (b) below the threshold. Range defined by bandwidth = %f\n", d_bw_s)  
    #if (C_b_r$V==0 | C_b_r$B==0) printf("{err}Not enough variability to compute the bias bandwidth (b) above the threshold. Range defined by bandwidth = %f\n", d_bw_s)  
    #if (C_h_l$V==0 | C_h_l$B==0) printf("{err}Not enough variability to compute the loc. poly. bandwidth (h) below the threshold. Range defined by bandwidth = %f\n", b_bw_s) 
    #if (C_h_r$V==0 | C_h_r$B==0) printf("{err}Not enough variability to compute the loc. poly. bandwidth (h) above the threshold. Range defined by bandwidth = %f\n", b_bw_s) 
}

    #                     *** RD
if  (bwselect=="mserd" | bwselect=="cerrd" | bwselect=="cerdpi" | bwselect=="msecomb1" | bwselect=="msecomb2" | bwselect=="cercomb1" | bwselect=="cercomb2" | bwselect=="" | all=="TRUE" ) {
  d_bw_d = c(( (C_d_l$V + C_d_r$V)  /  (C_d_r$B - C_d_l$B)^2 )^C_d_l$rate)
  C_b_l  = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c=c, o=q, nu=p+1, o_B=q+1, h_V=c_bw, h_B=d_bw_d, scaleregul, vce, nnmatch, kernel, dups_l, dupsid_l)
  C_b_r  = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c=c, o=q, nu=p+1, o_B=q+1, h_V=c_bw, h_B=d_bw_d, scaleregul, vce, nnmatch, kernel, dups_r, dupsid_r)
  b_bw_d = c(( (C_b_l$V + C_b_r$V)  /  ((C_b_r$B - C_b_l$B)^2 + scaleregul*(C_b_r$R + C_b_l$R)) )^C_b_l$rate)
  C_h_l  = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c=c, o=p, nu=deriv, o_B=q, h_V=c_bw, h_B=b_bw_d, scaleregul, vce, nnmatch, kernel, dups_l, dupsid_l)
  C_h_r  = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c=c, o=p, nu=deriv, o_B=q, h_V=c_bw, h_B=b_bw_d, scaleregul, vce, nnmatch, kernel, dups_r, dupsid_r)
  h_bw_d = c(( (C_h_l$V + C_h_r$V)  /  ((C_h_r$B - C_h_l$B)^2 + scaleregul*(C_h_r$R + C_h_l$R)) )^C_h_l$rate)
    
  #if (C_b_l$V==0 | C_b_l$B==0) printf("{err}Not enough variability to compute the bias bandwidth (b) below the threshold. Range defined by bandwidth = %f\n", d_bw_d)  
  #if (C_b_r$V==0 | C_b_r$B==0) printf("{err}Not enough variability to compute the bias bandwidth (b) above the threshold. Range defined by bandwidth = %f\n", d_bw_d)  
  #if (C_h_l$V==0 | C_h_l$B==0) printf("{err}Not enough variability to compute the loc. poly. bandwidth (h) below the threshold. Range defined by bandwidth = %f\n", b_bw_d) 
  #if (C_h_r$V==0 | C_h_r$B==0) printf("{err}Not enough variability to compute the loc. poly. bandwidth (h) above the threshold. Range defined by bandwidth = %f\n", b_bw_d) 
}	
                     
#if (C_b_l$V==. | C_b_l$B==. | C_b_l$R==.) printf("{err}Invertibility problem in the computation of bias bandwidth (b) below the threshold") 			
#if (C_b_r$V==. | C_b_r$B==. | C_b_r$R==.) printf("{err}Invertibility problem in the computation of bias bandwidth (b) above the threshold")  
#if (C_h_l$V==. | C_h_l$B==. | C_h_l$R==.) printf("{err}Invertibility problem in the computation of loc. poly. bandwidth (h) below the threshold") 
#if (C_h_r$V==. | C_h_r$B==. | C_h_r$R==.) printf("{err}Invertibility problem in the computation of loc. poly. bandwidth (h) above the threshold") 
         
if  (bwselect=="mserd" | bwselect=="cerrd" | bwselect=="cerdpi" | bwselect=="msecomb1" | bwselect=="msecomb2" | bwselect=="cercomb1" | bwselect=="cercomb2" | bwselect=="" | all=="TRUE" ) {
  h_mserd = x_sd*h_bw_d
  b_mserd = x_sd*b_bw_d
}	
if  (bwselect=="msesum" | bwselect=="cersum" |  bwselect=="msecomb1" | bwselect=="msecomb2" |  bwselect=="cercomb1" | bwselect=="cercomb2"  |  all=="TRUE")  {
  h_msesum = x_sd*h_bw_s
  b_msesum = x_sd*b_bw_s
}
if  (bwselect=="msetwo" |  bwselect=="certwo" | bwselect=="msecomb2" | bwselect=="cercomb2"  | all=="TRUE")  {		
  h_msetwo_l = x_sd*h_bw_l
  h_msetwo_r = x_sd*h_bw_r
  b_msetwo_l = x_sd*b_bw_l
  b_msetwo_r = x_sd*b_bw_r
}
if  (bwselect=="msecomb1" | bwselect=="cercomb1" | all=="TRUE" ) {
  h_msecomb1 = min(c(h_mserd,h_msesum))
  b_msecomb1 = min(c(b_mserd,b_msesum))
}
if  (bwselect=="msecomb2" | bwselect=="cercomb2" |  all=="TRUE" ) {
  h_msecomb2_l = median(c(h_mserd,h_msesum,h_msetwo_l))
  h_msecomb2_r = median(c(h_mserd,h_msesum,h_msetwo_r))
  b_msecomb2_l = median(c(b_mserd,b_msesum,b_msetwo_l))
  b_msecomb2_r = median(c(b_mserd,b_msesum,b_msetwo_r))
}
cer_h = N^(-(p/((3+p)*(3+2*p))))

if (!is.null(cluster)) {
  cer_h = (g_l+g_r)^(-(p/((3+p)*(3+2*p))))
}

#cer_b = N^(-(q/((3+q)*(3+2*q))))
cer_b = 1
	if  (bwselect=="cerrd" | all=="TRUE" ){
		h_cerrd = h_mserd*cer_h
		b_cerrd = b_mserd*cer_b
	}
	if  (bwselect=="cersum" | all=="TRUE" ){
		h_cersum = h_msesum*cer_h
		b_cersum=  b_msesum*cer_b
		}
	if  (bwselect=="certwo" | all=="TRUE" ){
		h_certwo_l   = h_msetwo_l*cer_h
		h_certwo_r   = h_msetwo_r*cer_h
		b_certwo_l   = b_msetwo_l*cer_b
		b_certwo_r   = b_msetwo_r*cer_b
		}
	if  (bwselect=="cercomb1" | all=="TRUE" ){
		h_cercomb1 = h_msecomb1*cer_h
		b_cercomb1 = b_msecomb1*cer_b
		}
	if  (bwselect=="cercomb2" | all=="TRUE" ){
		h_cercomb2_l = h_msecomb2_l*cer_h
		h_cercomb2_r = h_msecomb2_r*cer_h
		b_cercomb2_l = b_msecomb2_l*cer_b
		b_cercomb2_r = b_msecomb2_r*cer_b
	}

sigma.rbc = eta.bc = q.terms = NA

if  (bwselect=="cerdpi" |  all=="TRUE" ) {
  cer.out = rdbwselect.ce.dpi(y*y_sd, x*x_sd, h_mserd, h_mserd, c, p, q, deriv, rho, kernel, vce, nnmatch, tau, quant)
  #cer.out = rdbwselect.ce.dpi(y*y_sd, x*x_sd, h.mse.pob, h.mse.pob, c, p, q, deriv, rho, kernel, vce, nnmatch, tau)
  h_cerdpi = cer.out$h
  b_cerdpi = x_sd*b_bw_d

  #z = qnorm(0.975)
  sigma.rbc = cer.out$sigma
  eta.bc    = cer.out$eta.bc
  q.terms   = cer.out$q.terms
  #H.star = function(W) { (  ((1-W)/W)   * ((1+2*deriv)/(5+2*p)) * (sigma*z)/Q2 )^(2/(9+4*deriv+2*deriv))} 
}	


if (all=="FALSE"){
  bw_list = bwselect
  bws = matrix(NA,1,4)
  colnames(bws)=c("h (left)","h (right)","b (left)","b (right)")
  rownames(bws)=bwselect
  if  (bwselect=="mserd" | bwselect=="") bws[1,] = c(h_mserd,      h_mserd,      b_mserd,      b_mserd)
  if  (bwselect=="msetwo")               bws[1,] = c(h_msetwo_l,   h_msetwo_r,   b_msetwo_l,   b_msetwo_r)
  if  (bwselect=="msesum")               bws[1,] = c(h_msesum,     h_msesum,     b_msesum,     b_msesum)
  if  (bwselect=="msecomb1")             bws[1,] = c(h_msecomb1,   h_msecomb1,   b_msecomb1,   b_msecomb1)
  if  (bwselect=="msecomb2")             bws[1,] = c(h_msecomb2_l, h_msecomb2_r, b_msecomb2_l, b_msecomb2_r) 
  if  (bwselect=="cerrd")                bws[1,] = c(h_cerrd,      h_cerrd,      b_cerrd,      b_cerrd)
  if  (bwselect=="certwo")               bws[1,] = c(h_certwo_l,   h_certwo_r,   b_certwo_l,   b_certwo_r)
  if  (bwselect=="cersum")               bws[1,] = c(h_cersum,     h_cersum,     b_cersum,     b_cersum)
  if  (bwselect=="cercomb1")             bws[1,] = c(h_cercomb1,   h_cercomb1,   b_cercomb1,   b_cercomb1)
  if  (bwselect=="cercomb2")             bws[1,] = c(h_cercomb2_l, h_cercomb2_r, b_cercomb2_l, b_cercomb2_r)
  if  (bwselect=="cerdpi")               bws[1,] = c(h_cerdpi,     h_cerdpi,     b_cerdpi,     b_cerdpi)
}

  if (all=="TRUE"){
    bwselect="All"
    bws = matrix(NA,10,4)
    colnames(bws)=c("h (left)","h (right)","b (left)","b (right)")
    bw_list=c("mserd","msetwo","msesum","msecomb1","msecomb2","cerrd","certwo","cersum","cercomb1","cercomb2") 
    rownames(bws)=c("mserd","msetwo","msesum","msecomb1","msecomb2","cerrd","certwo","cersum","cercomb1","cercomb2") 
    bws[1,] =c(h_mserd,      h_mserd,      b_mserd,      b_mserd)
    bws[2,] =c(h_msetwo_l,   h_msetwo_r,   b_msetwo_l,   b_msetwo_r)
    bws[3,] =c(h_msesum,     h_msesum,     b_msesum,     b_msesum)
    bws[4,] =c(h_msecomb1,   h_msecomb1,   b_msecomb1,   b_msecomb1)
    bws[5,] =c(h_msecomb2_l, h_msecomb2_r, b_msecomb2_l, b_msecomb2_r) 
    bws[6,] =c(h_cerrd,      h_cerrd,      b_cerrd,      b_cerrd)
    bws[7,] =c(h_certwo_l,   h_certwo_r,   b_certwo_l,   b_certwo_r)
    bws[8,] =c(h_cersum,     h_cersum,     b_cersum,     b_cersum)
    bws[9,] =c(h_cercomb1,   h_cercomb1,   b_cercomb1,   b_cercomb1)
    bws[10,]=c(h_cercomb2_l, h_cercomb2_r, b_cercomb2_l, b_cercomb2_r)
  }
  

  out = list(bws=bws, 
             bwselect=bwselect, bw_list=bw_list, kernel=kernel_type, p=p, q=q, c=c, N=c(N_l,N_r), vce=vce_type, sigma.rbc=sigma.rbc, eta.bc=eta.bc, q.terms=q.terms)
  out$call <- match.call()
  class(out) <- "rdbwselect"
  return(out)
}

print.rdbwselect <- function(x,...){
  cat("Call: rdbwselect\n\n")
  cat(paste("Number of Obs.           ",  format(sprintf("%10.0f",x$N[1]+x$N[2], width=10, justify="right")),"\n", sep=""))
  cat(paste("BW type                  ",  format(x$bwselect, width=10, justify="right"),"\n", sep=""))
  cat(paste("Kernel                   ",  format(x$kernel,   width=10, justify="right"),"\n", sep=""))
  cat(paste("VCE method               ",  format(x$vce,      width=10, justify="right"),"\n", sep=""))
  cat("\n")
  cat(paste("Number of Obs.           ",  format(sprintf("%9.0f",x$N[1], width=10, justify="right")),  "   ", format(sprintf("%9.0f",x$N[2],width=10, justify="right")),        "\n", sep=""))
  cat(paste("Order est. (p)           ",  format(sprintf("%9.0f",x$p,    width=10, justify="right")),  "   ", format(sprintf("%9.0f",x$p,  width=10, justify="right")),       "\n", sep=""))
  cat(paste("Order bias  (p)          ",  format(sprintf("%9.0f",x$q,    width=10, justify="right")),  "   ", format(sprintf("%9.0f",x$q,  width=10, justify="right")),       "\n", sep=""))
  cat("\n")
}

summary.rdbwselect <- function(object,...) {
  x    <- object
  args <- list(...)
  
  cat("Call: rdbwselect\n\n")
  
  cat(paste("Number of Obs.           ",  format(sprintf("%10.0f",x$N[1]+x$N[2], width=10, justify="right")),"\n", sep=""))
  cat(paste("BW type                  ",  format(x$bwselect, width=10, justify="right"),"\n", sep=""))
  cat(paste("Kernel                   ",  format(x$kernel,   width=10, justify="right"),"\n", sep=""))
  cat(paste("VCE method               ",  format(x$vce,      width=10, justify="right"),"\n", sep=""))
  cat("\n")
  cat(paste("Number of Obs.           ",  format(sprintf("%9.0f",x$N[1], width=10, justify="right")),  "   ", format(sprintf("%9.0f",x$N[2],width=10, justify="right")),        "\n", sep=""))
  cat(paste("Order est. (p)           ",  format(sprintf("%9.0f",x$p,    width=10, justify="right")),  "   ", format(sprintf("%9.0f",x$p,  width=10, justify="right")),       "\n", sep=""))
  cat(paste("Order bias  (q)          ",  format(sprintf("%9.0f",x$q,    width=10, justify="right")),  "   ", format(sprintf("%9.0f",x$q,  width=10, justify="right")),       "\n", sep=""))
  cat("\n")
  
  
    col1.names = c("","BW est. (h)", "BW bias (b)")
    col2.names = c("","Left of c", "Right of c","Left of c", "Right of c")

  ### print output
  cat(paste(rep("=", 15 + 10*4), collapse="")); cat("\n")

    cat(format(col1.names  , width=14, justify="right"))
    cat("\n")
    cat(format(col2.names            , width=10, justify="right"))
    cat("\n")
    
  cat(paste(rep("=", 15 + 10*4), collapse="")); cat("\n")
  
    for (j in 1:nrow(x$bws)) {
      #cat(format(toString(j), width=4))
      cat(format(x$bw_list[j]           , width=10, justify="right"))
      cat(format(sprintf("%3.3f", x$bws[j,]), width=10, justify="right"))
      cat("\n")
    }
    cat(paste(rep("=", 15 + 10*ncol(x$bws)), collapse="")); cat("\n")   
}


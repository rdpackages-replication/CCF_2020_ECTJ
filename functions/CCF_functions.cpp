#include <RcppArmadillo.h>
#include <iostream> 

using namespace Rcpp;

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

List rdbwce(arma::vec y_l, arma::vec y_r, arma::vec x_l, arma::vec x_r, arma::vec K_l, arma::vec K_r, arma::vec L_l, arma::vec L_r, arma::vec res_l, arma::vec res_r, double c, int p, int q, double h, double b, int deriv, int fact, double z) {

  int N   = y_l.n_rows;
  int N_l = y_l.n_rows;
  int N_r = y_r.n_rows;
  float rho = h / b;

  arma::vec Wp_l = K_l/h;
  arma::vec Wp_r = K_r/h;
  arma::vec Wq_l = L_l/b;
  arma::vec Wq_r = L_r/b;
  arma::vec Xh_l = (x_l-c)/h;
  arma::vec Xh_r = (x_r-c)/h;
  arma::vec Xb_l = (x_l-c)/b;
  arma::vec Xb_r = (x_r-c)/b;
  arma::mat Rq_l(N_l,q+1,arma::fill::zeros);
  arma::mat Rq_r(N_r,q+1,arma::fill::zeros);
  arma::mat Rp_l(N_l,p+1,arma::fill::zeros);
  arma::mat Rp_r(N_r,p+1,arma::fill::zeros);
  for (int i=0; i<q+1 ; i++) Rq_l.col(i) = pow(Xb_l,i);
  for (int i=0; i<q+1 ; i++) Rq_r.col(i) = pow(Xb_r,i);
  for (int i=0; i<p+1 ; i++) Rp_l.col(i) = pow(Xh_l,i);
  for (int i=0; i<p+1 ; i++) Rp_r.col(i) = pow(Xh_r,i);

  arma::mat dWp_l = diagmat(Wp_l);
  arma::mat dWp_r = diagmat(Wp_r);
  arma::mat dWq_l = diagmat(Wq_l);
  arma::mat dWq_r = diagmat(Wq_r);
  arma::mat dK_l  = diagmat(K_l);
  arma::mat dK_r  = diagmat(K_r);
  arma::mat dL_l  = diagmat(L_l);
  arma::mat dL_r  = diagmat(L_r);

  arma::mat Lp1_l = Rp_l.t()*dWp_l*pow(Xh_l,p+1)/N_l;
  arma::mat Lp1_r = Rp_r.t()*dWp_r*pow(Xh_r,p+1)/N_r;
  arma::mat Gp_l  = Rp_l.t()*dWp_l*Rp_l/N_l;
  arma::mat Gp_r  = Rp_r.t()*dWp_r*Rp_r/N_r;
  arma::mat Gq_l  = Rq_l.t()*dWq_l*Rq_l/N_l;
  arma::mat Gq_r  = Rq_r.t()*dWq_r*Rq_r/N_r;

  arma::mat iGp_l = Gp_l.i();
  arma::mat iGp_r = Gp_r.i();
  arma::mat iGq_l = Gq_l.i();
  arma::mat iGq_r = Gq_r.i();

  arma::mat ep1(q+1,1,arma::fill::zeros);
  ep1(p+1) = 1;
  arma::mat e0(p+1,1,arma::fill::zeros);
  e0(deriv) = fact;

  arma::mat lus0_l = e0.t()*iGp_l*(dK_l*Rp_l).t();
  arma::mat lus0_r = e0.t()*iGp_r*(dK_r*Rp_r).t();
  arma::mat lus0 = lus0_r - lus0_l;
  arma::mat lbc0_l = lus0_l - pow(rho,p+1)*(e0.t()*iGp_l)*Lp1_l*ep1.t()*iGq_l*(dL_l*Rq_l).t();
  arma::mat lbc0_r = lus0_r - pow(rho,p+1)*(e0.t()*iGp_r)*Lp1_r*ep1.t()*iGq_r*(dL_r*Rq_r).t();
  arma::mat lbc0 = lbc0_r - lbc0_l;
  arma::vec vx_l = pow(res_l,2);
  arma::vec vx_r = pow(res_r,2);
  arma::vec vx = vx_l + vx_r;
  arma::vec res = res_l + res_r;
  
  double sums2_l=0;
  for (int i = 0 ; i < N_l ; i++) {
    sums2_l += pow(lbc0_l(i),2)*vx_l(i);
  }
  double s2_l=sums2_l/(N_l*h);
  
  double sums2_r=0;
  for (int i = 0 ; i < N_r ; i++) {
    sums2_r += pow(lbc0_r(i),2)*vx_r(i);
  }
  double s2_r=sums2_r/(N_r*h);
  double s2 = s2_l + s2_r;

  arma::mat Krrp_l(p+1,p+1,arma::fill::zeros);
  arma::mat Krrp_r(p+1,p+1,arma::fill::zeros);
  arma::mat Krxip_l(1, p+1,arma::fill::zeros);
  arma::mat Krxip_r(1, p+1,arma::fill::zeros);
  arma::mat Krxp_l(1,  p+1,arma::fill::zeros);
  arma::mat Krxp_r(1,  p+1,arma::fill::zeros);
  arma::mat Lrrq_l(q+1,q+1,arma::fill::zeros);
  arma::mat Lrrq_r(q+1,q+1,arma::fill::zeros);

  for (int i = 0 ; i < N_l ; i++) {
    arma::mat Rpi_l = Rp_l.row(i);
    arma::mat Rpi_r = Rp_r.row(i);
    arma::mat Rqi_l = Rq_l.row(i);
    arma::mat Rqi_r = Rq_r.row(i);
    Krrp_l  += K_l(i)*Rpi_l.t()*Rpi_l;
    Krrp_r  += K_r(i)*Rpi_r.t()*Rpi_r;
    Lrrq_l  += L_l(i)*Rqi_l.t()*Rqi_l;
    Lrrq_r  += L_r(i)*Rqi_r.t()*Rqi_r;
    Krxip_l += K_l(i)*Rpi_l*pow(Xh_l(i),p+1) ;
    Krxip_r += K_r(i)*Rpi_r*pow(Xh_r(i),p+1) ;
    for (int j = 0 ; j < N_l ; j++) {
      if (j != i) Krxp_l += K_l(i)*Rpi_l*pow(Xh_l(j),p+1);
    }
    for (int j = 0 ; j < N_r ; j++) {
      if (j != i) Krxp_r += K_r(i)*Rpi_r*pow(Xh_r(j),p+1);
    }
  }
  
  arma::mat EKrrp_l = Krrp_l/N_l;
  arma::mat EKrrp_r = Krrp_r/N_r;
  arma::mat EKrxp_l = Krxp_l/(N_l*(N_l-1));
  arma::mat EKrxp_r = Krxp_r/(N_r*(N_r-1));
  arma::mat EKrxip_l = Krxip_l/N_l;
  arma::mat EKrxip_r = Krxip_r/N_r;
  arma::mat ELrrq_l = Lrrq_l/N_l;
  arma::mat ELrrq_r = Lrrq_r/N_r;

  arma::mat q1(1,1,arma::fill::zeros);
  arma::mat q2(1,1,arma::fill::zeros);
  arma::mat q3(1,1,arma::fill::zeros);
  arma::mat q4l(1,1,arma::fill::zeros);
  arma::mat q4r(1,1,arma::fill::zeros);
  arma::mat q5al(1,q+1,arma::fill::zeros);
  arma::mat q5ar(1,q+1,arma::fill::zeros);
  arma::mat q5bl(q+1,1,arma::fill::zeros);
  arma::mat q5br(q+1,1,arma::fill::zeros);
  arma::mat q6l(1,1,arma::fill::zeros);
  arma::mat q6r(1,1,arma::fill::zeros);
  arma::mat q7al(1,q+1,arma::fill::zeros);
  arma::mat q7ar(1,q+1,arma::fill::zeros);
  arma::mat q7bl(q+1,q+1,arma::fill::zeros);
  arma::mat q7br(q+1,q+1,arma::fill::zeros);
  arma::mat q7cl(q+1,1,arma::fill::zeros);
  arma::mat q7cr(q+1,1,arma::fill::zeros);
  arma::mat q8(1,1,arma::fill::zeros);
  arma::mat q9(1,1,arma::fill::zeros);
  arma::mat q10(1,1,arma::fill::zeros);
  arma::mat q11(1,1,arma::fill::zeros);
  arma::mat q12(1,1,arma::fill::zeros);
  arma::mat q3a(1,1,arma::fill::zeros);

  for (int i = 0 ; i < N ; i++) {
    arma::mat Rpi_l = Rp_l.row(i);
    arma::mat Rpi_r = Rp_r.row(i);
    arma::mat Rqi_l = Rq_l.row(i);
    arma::mat Rqi_r = Rq_r.row(i);

    q1 += pow(lbc0(i)*res(i),3);

    arma::mat lus1_l = fact * iGp_l.row(deriv) * (EKrrp_l - K_l(i)*Rpi_l.t()*Rpi_l)*iGp_l*K_l(i)*Rpi_l.t();
    arma::mat lus1_r = fact * iGp_r.row(deriv) * (EKrrp_r - K_r(i)*Rpi_r.t()*Rpi_r)*iGp_r*K_r(i)*Rpi_r.t();
    
    arma::mat T1_l   = fact * iGp_l.row(deriv) * ((EKrrp_l - K_l(i)*Rpi_l.t()*Rpi_l)*iGp_l*Lp1_l*ep1.t())    *(iGq_l*L_l(i)*Rqi_l.t());
    arma::mat T2_l   = fact * iGp_l.row(deriv) * ((K_l(i)*Rpi_l*pow(Xh_l(i),p+1) - EKrxip_l).t())*ep1.t()*(iGq_l*L_l(i)*Rqi_l.t());
    arma::mat T3_l   = fact * iGp_l.row(deriv) * ((Lp1_l*ep1.t()*iGq_l)*(ELrrq_l - L_l(i)*Rqi_l.t()*Rqi_l))  *(iGq_l*L_l(i)*Rqi_l.t());
    arma::mat lbc1_l = lus1_l - pow(rho,p+1)*(T1_l + T2_l + T3_l);
    arma::mat T1_r   = fact * iGp_r.row(deriv) * ((EKrrp_r - K_r(i)*Rpi_r.t()*Rpi_r)*iGp_r*Lp1_r*ep1.t())    *(iGq_r*L_r(i)*Rqi_r.t());
    arma::mat T2_r   = fact * iGp_r.row(deriv) * ((K_r(i)*Rpi_r*pow(Xh_r(i),p+1) - EKrxip_r).t())*ep1.t()*(iGq_r*L_r(i)*Rqi_r.t());
    arma::mat T3_r   = fact * iGp_r.row(deriv) * ((Lp1_r*ep1.t()*iGq_r)*(ELrrq_r - L_r(i)*Rqi_r.t()*Rqi_r))  *(iGq_r*L_r(i)*Rqi_r.t());
    arma::mat lbc1_r = lus1_r - pow(rho,p+1)*(T1_r + T2_r + T3_r);
    arma::mat lbc1 =  lbc1_r - lbc1_l;
    
    q2 += lbc1*lbc0(i)*pow(res(i),2);

    q3 += pow(lbc0(i),4)*(pow(res(i),4)-pow(vx(i),2));

    q4l += pow(lbc0(i),2)*(Rqi_l*iGq_l*L_l(i)*Rqi_l.t())*pow(res_l(i),2);
    q4r += pow(lbc0(i),2)*(Rqi_r*iGq_r*L_r(i)*Rqi_r.t())*pow(res_r(i),2);

    q5al += pow(lbc0(i),3)*(Rqi_l*iGq_l)*pow(res_l(i),2);
    q5ar += pow(lbc0(i),3)*(Rqi_r*iGq_r)*pow(res_r(i),2);
    q5bl += L_l(i)*Rqi_l.t()*lbc0(i)*pow(res_l(i),2);
    q5br += L_r(i)*Rqi_r.t()*lbc0(i)*pow(res_r(i),2);

    q7al += lbc0(i)*pow(res_l(i),2)*L_l(i)*Rqi_l*iGq_l;
    q7ar += lbc0(i)*pow(res_r(i),2)*L_r(i)*Rqi_r*iGq_r;
    q7bl += pow(lbc0(i),2)*Rqi_l.t()*Rqi_l*iGq_l;
    q7br += pow(lbc0(i),2)*Rqi_r.t()*Rqi_r*iGq_r;
    q7cl += lbc0(i)*pow(res_l(i),2)*L_l(i)*Rqi_l.t();
    q7cr += lbc0(i)*pow(res_r(i),2)*L_r(i)*Rqi_r.t();

    q8  += pow(lbc0(i)*res(i),4);
    q9  += (pow(lbc0(i),2)*vx(i)-h*s2)*pow(lbc0(i)*res(i),2);

    q12 += pow(pow(lbc0(i),2)*vx(i)-h*s2,2);

    q3a += pow(lbc0[i]*res[i],3);


    for (int j = 0 ; j < N ; j++) {
      if (j != i) {
        arma::mat Rpj_l = Rp_l.row(j);
        arma::mat Rqj_l = Rq_l.row(j);
        arma::mat Rpj_r = Rp_r.row(j);
        arma::mat Rqj_r = Rq_r.row(j);

        arma::mat lus1_l = fact * iGp_l.row(deriv) *  (EKrrp_l - K_l(j)*Rpj_l.t()*Rpj_l)*iGp_l*K_l(i)*Rpi_l.t();
        arma::mat T1_l   = fact * iGp_l.row(deriv) * ((EKrrp_l - K_l(j)*Rpj_l.t()*Rpj_l)*iGp_l*Lp1_l*ep1.t())    *(iGq_l*L_l(i)*Rqi_l.t());
        arma::mat T2_l   = fact * iGp_l.row(deriv) * ((K_l(j)*Rpj_l*pow(Xh_l(i),p+1) - EKrxp_l).t())*ep1.t() *(iGq_l*L_l(i)*Rqi_l.t());
        arma::mat T3_l   = fact * iGp_l.row(deriv) * ((Lp1_l*ep1.t()*iGq_l)*(ELrrq_l - L_l(j)*Rqj_l.t()*Rqj_l))  *(iGq_l*L_l(i)*Rqi_l.t());
        arma::mat lbc1_l = lus1_l - pow(rho,p+1)*(T1_l + T2_l + T3_l);
        
        arma::mat lus1_r = fact * iGp_r.row(deriv) *  (EKrrp_r - K_r(j)*Rpj_r.t()*Rpj_r)*iGp_r*K_r(i)*Rpi_r.t();
        arma::mat T1_r   = fact * iGp_r.row(deriv) * ((EKrrp_r - K_r(j)*Rpj_r.t()*Rpj_r)*iGp_r*Lp1_r*ep1.t())    *(iGq_r*L_r(i)*Rqi_r.t());
        arma::mat T2_r   = fact * iGp_r.row(deriv) * ((K_r(j)*Rpj_r*pow(Xh_r(i),p+1) - EKrxp_r).t())*ep1.t() *(iGq_r*L_r(i)*Rqi_r.t());
        arma::mat T3_r   = fact * iGp_r.row(deriv) * ((Lp1_r*ep1.t()*iGq_r)*(ELrrq_r - L_r(j)*Rqj_r.t()*Rqj_r))  *(iGq_r*L_r(i)*Rqi_r.t());
        arma::mat lbc1_r = lus1_r - pow(rho,p+1)*(T1_r + T2_r + T3_r);

        arma::mat lbc1 =  lbc1_r - lbc1_l;
        
        q10 += lbc1*lbc0(i)*pow(lbc0(j)*res(j),2)*vx(i);
        q11 += lbc1*lbc0(i)*(pow(lbc0(j),2)*vx(j)-h*s2)*pow(res(i),2);

        q6l += pow(lbc0(i),2)*pow(Rqi_l*iGq_l*L_l(j)*Rqj_l.t(),2)*pow(res_l(j),2);
        q6r += pow(lbc0(i),2)*pow(Rqi_r*iGq_r*L_r(j)*Rqj_r.t(),2)*pow(res_r(j),2);
       }
    }
  }

  arma::mat Eq1  = pow(q1/(N*h),2);
  arma::mat Eq2  = q2/(N*h);
  arma::mat Eq3  = q3/(N*h);
  arma::mat Eq4  = q4l/(N_l*h) + q4r/(N_r*h);
  arma::mat Eq5  = (q5al/(N_l*h))*(q5bl/(N_l*h)) + (q5ar/(N_r*h))*(q5br/(N_r*h));
  arma::mat Eq6  = q6l/(N_l*(N_l-1)*pow(h,2)) + q6r/(N_r*(N_r-1)*pow(h,2));
  arma::mat Eq7  = (q7al/(N_l*h))*(q7bl/(N_l*h))*(q7cl/(N_l*h)) + (q7ar/(N_r*h))*(q7br/(N_r*h))*(q7cr/(N_r*h));
  arma::mat Eq8  = q8/(N*h);
  arma::mat Eq9  = q9/(N*h);
  arma::mat Eq10 = q10/(N*(N-1)*pow(h,2));
  arma::mat Eq11 = q11/(N*(N-1)*pow(h,2));
  arma::mat Eq12 = q12/(N*h);

  arma::mat  q1rbc =   
       Eq1*(pow(z,3)/3+7*z/4+s2*z*(pow(z,2)-3)/4)/pow(s2,3)
    +  Eq2*(-z*(pow(z,2)-3)/2)/s2
    +  Eq3*(z*(pow(z,2)-3)/8)/pow(s2,2)
    -  Eq4*(z*(pow(z,2)-1)/2)/s2
    -  Eq5*(z*(pow(z,2)-1))/pow(s2,2)
    +  Eq6*(z*(pow(z,2)-1)/4)/s2
    +  Eq7*(z*(pow(z,2)-1)/2)/pow(s2,2)
    +  Eq8*(-z*(pow(z,2)-3)/24)/pow(s2,2)
    +  Eq9*(z*(pow(z,2)-1)/4)/pow(s2,2)
    +  Eq10*(z*(pow(z,2)-3))/pow(s2,2)
    +  Eq11*(-z)/pow(s2,2)
    +  Eq12*(-z*(pow(z,2)+1)/8)/pow(s2,2);

    double q2rbc = -z/(2*s2);

    arma::mat Eq3a  = q3a/(N*h);
    arma::mat q3rbc = Eq3a/(pow(s2,2))*(pow(z,3)/3);

  return Rcpp::List::create(
    Rcpp::Named("q1rbc") = q1rbc,
    Rcpp::Named("q2rbc") = q2rbc,
    Rcpp::Named("q3rbc") = q3rbc,
    Rcpp::Named("s2rbc") = s2
  ) ;

}


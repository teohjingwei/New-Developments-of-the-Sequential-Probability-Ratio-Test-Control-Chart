#include <Rcpp.h>
#include <iostream>
#include <cmath>
using namespace std;
using namespace Rcpp;

// [[Rcpp::depends(BH)]]
#include <boost/math/distributions/non_central_chi_squared.hpp>

#define M 200

double weight[9] = {0.0812743884, 0.1806481607, 0.2606106964, 0.3123470770, 0.3302393550, 0.3123470770, 0.2606106964, 0.1806481607, 0.0812743884};
double abscissa[9] = {-0.9681602395, -0.8360311073, -0.6133714327, -0.3242534234, 0., 0.3242534234,  0.6133714327,  0.8360311073,  0.9681602395};

/*
#################################
###   Run-length properties   ###
#################################
 */

double arl_sprt(int method, double g, double h, double k, double gamma, double deltamu, double deltasigma, double C[], double R[], double N[])
{ double arl, P0;
  double pnorm(double x), pnchisq(double x, double df, double ncp);
  void Rmat(int method, double g, double h, double k, double gamma, double deltamu, double deltasigma, double R[]), Qvec(int method, double g, double h, double k, double gamma, double deltamu, double deltasigma, double C[]), Nvec(int cases, double C[], double R[], double N[]), Cvec(int method, double g, double h, double k, double gamma, double deltamu, double deltasigma, double C[]);
  
  Rmat(method, g, h, k, gamma, deltamu, deltasigma, R);
  
  Qvec(method, g, h, k, gamma, deltamu, deltasigma, C);
  
  Nvec(2, C, R, N);
  
  Cvec(method, g, h, k, gamma, deltamu, deltasigma, C);
  
  if(method == 1){
    P0 = (g + gamma > 0.0) ? pnchisq((g + gamma)/(deltasigma * deltasigma),1,(deltamu+k)*(deltamu+k)/(deltasigma * deltasigma)) : 0.0;
  }else if(method == 2){
    P0 = (g + gamma > 0.0) ? pnorm((sqrt(g + gamma) - (deltamu + k))/deltasigma) - pnorm((-sqrt(g + gamma) - (deltamu + k))/deltasigma) : 0.0;
  }
  
  for(int i = 0; i < M; ++i)
  { P0 += C[i] * N[i];
  }
  
  arl = 1./(1. - P0);
  
  return arl;
}

double ats_sprt(int method, double g, double h, double k, double gamma, double deltamu, double deltasigma, double d, double C[], double R[], double N[])
{ double ats, P0;
  double pnorm(double x), pnchisq(double x, double df, double ncp);
  void Rmat(int method, double g, double h, double k, double gamma, double deltamu, double deltasigma, double R[]), Qvec(int method, double g, double h, double k, double gamma, double deltamu, double deltasigma, double C[]), Nvec(int cases, double C[], double R[], double N[]), Cvec(int method, double g, double h, double k, double gamma, double deltamu, double deltasigma, double C[]);
  
  Rmat(method, g, h, k, gamma, deltamu, deltasigma, R);
  
  Qvec(method, g, h, k, gamma, deltamu, deltasigma, C);
  
  Nvec(2, C, R, N);
  
  Cvec(method, g, h, k, gamma, deltamu, deltasigma, C);
  
  if(method == 1){
    P0 = (g + gamma > 0.0) ? pnchisq((g + gamma)/(deltasigma * deltasigma),1,(deltamu+k)*(deltamu+k)/(deltasigma * deltasigma)) : 0.0;
  }else if(method == 2){
    P0 = (g + gamma > 0.0) ? pnorm((sqrt(g + gamma) - (deltamu + k))/deltasigma) - pnorm((-sqrt(g + gamma) - (deltamu + k))/deltasigma) : 0.0;
  }
  
  for(int i = 0; i < M; ++i)
  { P0 += C[i] * N[i];
  }
  
  ats = d/(1. - P0);
  
  if((deltamu != 0.0) || (deltasigma != 1.0))
    ats -= d * 0.5; 
  
  return ats;
}

double sdrl_sprt(int method, double g, double h, double k, double gamma, double deltamu, double deltasigma, double C[], double R[], double N[])
{ double sdrl, P0;
  double pnorm(double x), pnchisq(double x, double df, double ncp);
  void Rmat(int method, double g, double h, double k, double gamma, double deltamu, double deltasigma, double R[]), Qvec(int method, double g, double h, double k, double gamma, double deltamu, double deltasigma, double C[]), Nvec(int cases, double C[], double R[], double N[]), Cvec(int method, double g, double h, double k, double gamma, double deltamu, double deltasigma, double C[]);
  
  Rmat(method, g, h, k, gamma, deltamu, deltasigma, R);
  
  Qvec(method, g, h, k, gamma, deltamu, deltasigma, C);
  
  Nvec(2, C, R, N);
  
  Cvec(method, g, h, k, gamma, deltamu, deltasigma, C);
  
  if(method == 1){
    P0 = (g + gamma > 0.0) ? pnchisq((g + gamma)/(deltasigma * deltasigma),1,(deltamu+k)*(deltamu+k)/(deltasigma * deltasigma)) : 0.0;
  }else if(method == 2){
    P0 = (g + gamma > 0.0) ? pnorm((sqrt(g + gamma) - (deltamu + k))/deltasigma) - pnorm((-sqrt(g + gamma) - (deltamu + k))/deltasigma) : 0.0;
  }
  
  for(int i = 0; i < M; ++i)
  { P0 += C[i] * N[i];
  }
  
  sdrl = sqrt(P0/((1.0 - P0)*(1.0 - P0)));
  
  return sdrl;
}

double sdts_sprt(int method, double g, double h, double k, double gamma, double deltamu, double deltasigma, double d, double C[], double R[], double N[])
{ double sdts, P0;
  double pnorm(double x), pnchisq(double x, double df, double ncp);
  void Rmat(int method, double g, double h, double k, double gamma, double deltamu, double deltasigma, double R[]), Qvec(int method, double g, double h, double k, double gamma, double deltamu, double deltasigma, double C[]), Nvec(int cases, double C[], double R[], double N[]), Cvec(int method, double g, double h, double k, double gamma, double deltamu, double deltasigma, double C[]);
  
  Rmat(method, g, h, k, gamma, deltamu, deltasigma, R);
  
  Qvec(method, g, h, k, gamma, deltamu, deltasigma, C);
  
  Nvec(2, C, R, N);
  
  Cvec(method, g, h, k, gamma, deltamu, deltasigma, C);
  
  if(method == 1){
    P0 = (g + gamma > 0.0) ? pnchisq((g + gamma)/(deltasigma * deltasigma),1,(deltamu+k)*(deltamu+k)/(deltasigma * deltasigma)) : 0.0;
  }else if(method == 2){
    P0 = (g + gamma > 0.0) ? pnorm((sqrt(g + gamma) - (deltamu + k))/deltasigma) - pnorm((-sqrt(g + gamma) - (deltamu + k))/deltasigma) : 0.0;
  }
  
  for(int i = 0; i < M; ++i)
  { P0 += C[i] * N[i];
  }
  
  if((deltamu == 0.0) && (deltasigma == 1.0)){
    sdts = d * sqrt(P0/((1.0 - P0)*(1.0 - P0)));
  }else{
    sdts = d * sqrt(1.0/12.0 + P0/((1.0 - P0)*(1.0 - P0)));
  }
  
  return sdts;
}

double asn_sprt(int method, double g, double h, double k, double gamma, double deltamu, double deltasigma, double C[], double R[], double N[])
{ 
  double pnorm(double x), pnchisq(double x, double df, double ncp);
  void Rmat(int method, double g, double h, double k, double gamma, double deltamu, double deltasigma, double R[]), Qvec(int method, double g, double h, double k, double gamma, double deltamu, double deltasigma, double C[]), Nvec(int cases, double C[], double R[], double N[]), Cvec(int method, double g, double h, double k, double gamma, double deltamu, double deltasigma, double C[]);
  
  Rmat(method, g, h, k, gamma, deltamu, deltasigma, R);
  
  Qvec(method, g, h, k, gamma, deltamu, deltasigma, C);
  
  Nvec(1, C, R, N);
  
  Cvec(method, g, h, k, gamma, deltamu, deltasigma, C);
  
  double asn = 1.0;
  for(int i = 0; i < M; ++i)
  { asn += C[i] * N[i];
  }
  
  return asn;
}

double anos_sprt(int method, double g, double h, double k, double gamma, double deltamu, double deltasigma, double C[], double R[], double N[])
{
  return arl_sprt(method, g, h, k, gamma, deltamu, deltasigma, C, R, N) * asn_sprt(method, g, h, k, gamma, deltamu, deltasigma, C, R, N);
}

double aeql_sprt(int method, double g, double h, double k, double gamma, double d, double C[], double R[], double N[], double deltamumin, double deltamumax, double deltasigmamin, double deltasigmamax)
{
  double sum = 0.0;
  for(int j = 0; j < 9; ++j){
    for(int i = 0; i < 9; ++i){
      double deltamuval = 0.5 * ((deltamumax - deltamumin)*abscissa[i] + (deltamumax + deltamumin));
      double deltasigmaval = 0.5 * ((deltasigmamax - deltasigmamin)*abscissa[j] + (deltasigmamax + deltasigmamin));
      sum += weight[j] * weight[i] * (deltamuval * deltamuval + deltasigmaval * deltasigmaval - 1.0) * ats_sprt(method, g, h, k, gamma, deltamuval, deltasigmaval, d, C, R, N);
    }
  }
  sum *= 0.25;
  return sum;
}

double eanos_sprt(int method, double g, double h, double k, double gamma, double C[], double R[], double N[], double deltamumin, double deltamumax, double deltasigmamin, double deltasigmamax)
{
  double sum = 0.0;
  for(int j = 0; j < 9; ++j){
    for(int i = 0; i < 9; ++i){
      double deltamuval = 0.5 * ((deltamumax - deltamumin)*abscissa[i] + (deltamumax + deltamumin));
      double deltasigmaval = 0.5 * ((deltasigmamax - deltasigmamin)*abscissa[j] + (deltasigmamax + deltasigmamin));
      sum += weight[j] * weight[i] * anos_sprt(method, g, h, k, gamma, deltamuval, deltasigmaval, C, R, N);
    }
  }
  sum *= 0.25;
  return sum;
}


/*
##############################
###   Vectors & Matrices   ###
##############################
 */

void Rmat(int method, double g, double h, double k, double gamma, double deltamu, double deltasigma, double R[])
{ int rn, kk, j;
  double start, temp, last, U, ncp, mu, valstart;
  double pnorm(double x), pnchisq(double x, double df, double ncp);
  
  double r = (h - g) / M;
  
  if(method == 1){
    ncp = (deltamu+k)*(deltamu+k)/(deltasigma * deltasigma);
    U = (-0.5*r + gamma)/(deltasigma * deltasigma);
    start = (U > 0.0) ? pnchisq(U,1,ncp) : 0.0;
    
    last = start;
    for (int i = 0; i < M; ++i)
    { U += r/(deltasigma * deltasigma);
      temp = (U > 0.0) ? pnchisq(U,1,ncp) : 0.0;
      R[i] = temp - last;
      last = temp;
    }
    
    last = start;
    U = (-0.5*r + gamma)/(deltasigma * deltasigma);
    for (int i = 1; i < M; ++i)
    { U -= r/(deltasigma * deltasigma);
      temp = (U > 0.0) ? pnchisq(U,1,ncp) : 0.0;
      R[i * M] = last - temp;
      last = temp;
    }
  }else if(method == 2){
    mu = (deltamu + k)/deltasigma;
    U = -0.5*r + gamma;
    start = (U > 0.0) ? pnorm(sqrt(U)/deltasigma - mu) - pnorm(-sqrt(U)/deltasigma - mu) : 0.0;
    
    if(U > 0.0){
      last = start;
      for (int i = 0; i < M; ++i)
      { U += r;
        temp = pnorm(sqrt(U)/deltasigma - mu) - pnorm(-sqrt(U)/deltasigma - mu);
        R[i] = temp - last;
        last = temp;
      }
    }else{
      if((static_cast<double>(M) - 0.5)*r + gamma < 0.0){
        for(int i = 0; i < M; ++i){
          R[i] = 0.0;
        }
      }else{
        for(int i = 0; i < M; ++i){
          U += r;
          R[i] = 0.0;
          if(U > 0){
            last = pnorm(sqrt(U)/deltasigma - mu) - pnorm(-sqrt(U)/deltasigma - mu);
            R[i] = last;
            j = i;
            break;
          }
        }
        for(int i = (j + 1); i < M; ++i){
          U += r;
          temp = pnorm(sqrt(U)/deltasigma - mu) - pnorm(-sqrt(U)/deltasigma - mu);
          R[i] = temp - last;
          last = temp;
        }
      }
    }
    
    U = -0.5*r + gamma;
    
    if(U <= 0.0){
      for (int i = 1; i < M; ++i)
      { R[i * M] = 0.0;
      }
    }else{
      valstart = U;
      last = pnorm(sqrt(U)/deltasigma - mu) - pnorm(-sqrt(U)/deltasigma - mu);
      
      for (int i = 1; i < M; ++i)
      { U -= r;
        
        if(valstart > 0.0){
          if(U <= 0.0){
            R[i * M] = pnorm(sqrt(valstart)/deltasigma - mu) - pnorm(-sqrt(valstart)/deltasigma - mu);
          }else{
            temp = pnorm(sqrt(U)/deltasigma - mu) - pnorm(-sqrt(U)/deltasigma - mu);
            R[i * M] = last - temp;
            last = temp;
          }
        }else{
          R[i * M] = 0.0;
        }
        valstart = U;
      }
    }
  }
  for (int i = 1; i < M; ++i)
  { rn = i * M;
    for (int j = 1; j < M; ++j)
    { kk = rn + j;
      R[kk] = R[kk-M-1];
    }
  }
}

void Qvec(int method, double g, double h, double k, double gamma, double deltamu, double deltasigma, double C[])
{ double U;
  double pnorm(double x), pnchisq(double x, double df, double ncp);
  
  double r = (h - g) / M;
  
  for(int i = 0; i < M; ++i){ 
    if(i == 0){
      U = gamma - 0.5 * r;
    }else{
      U -= r;
    }
    if(method == 1){
      C[i] = (U > 0.0) ? pnchisq(U/(deltasigma * deltasigma),1,(deltamu+k)*(deltamu+k)/(deltasigma * deltasigma)) : 0.0;
    }else if(method == 2){
      C[i] = (U > 0.0) ? pnorm((sqrt(U)-(deltamu + k))/deltasigma) - pnorm((-sqrt(U)-(deltamu + k))/deltasigma) : 0.0;
    }
  }
}

void Nvec(int cases, double C[], double R[], double N[])
{
  double g_equation(double *a, double *b, double *x, int n, int kk);
  
  for(int i = 0; i < M; ++i)
  { int ij = i * M - 1;
    
    for(int j = 0; j < M; ++j)
    { ++ij;
      if(i == j)
        R[ij] = 1.0 - R[ij];
      else 
        R[ij] = -R[ij];
    }
    
    if (cases == 1)
      C[i] = 1.;
  }
  
  g_equation(R, C, N, M, -1);
}

void Cvec(int method, double g, double h, double k, double gamma, double deltamu, double deltasigma, double C[])
{
  int j;
  double temp, start, last, U, ncp, mu;
  double pnorm(double x), pnchisq(double x, double df, double ncp);
  
  double r = (h - g) / M;
  if(method == 1){
    ncp = (deltamu+k)*(deltamu+k)/(deltasigma * deltasigma);
    U = (g + gamma)/(deltasigma * deltasigma);
    last = (U > 0.0) ? pnchisq(U,1,ncp) : 0.0;
    
    for(int i = 0; i < M; ++i){
      U += r/(deltasigma * deltasigma);
      temp = (U > 0.0) ? pnchisq(U,1,ncp) : 0.0;
      C[i] = temp - last;
      last = temp;
    }
  }else if(method == 2){
    
    mu = (deltamu + k)/deltasigma;
    U = g + gamma;
    start = (U > 0.0) ? pnorm(sqrt(U)/deltasigma - mu) - pnorm(-sqrt(U)/deltasigma - mu) : 0.0;
    
    if(U > 0.0){
      last = start;
      for (int i = 0; i < M; ++i)
      { U += r;
        temp = pnorm(sqrt(U)/deltasigma - mu) - pnorm(-sqrt(U)/deltasigma - mu);
        C[i] = temp - last;
        last = temp;
      }
    }else{
      if(g + static_cast<double>(M)*r + gamma > 0.0){
        for(int i = 0; i < M; ++i){
          U += r;
          C[i] = 0.0;
          if(U > 0){
            last = pnorm(sqrt(U)/deltasigma - mu) - pnorm(-sqrt(U)/deltasigma - mu);
            C[i] = last;
            j = i;
            break;
          }
        }
        for(int i = (j + 1); i < M; ++i){
          U += r;
          temp = pnorm(sqrt(U)/deltasigma - mu) - pnorm(-sqrt(U)/deltasigma - mu);
          C[i] = temp - last;
          last = temp;
        }
      }else{
        for(int i = 0; i < M; ++i){
          C[i] = 0.0;
        }
      }
    }
  }
}

/*
##############################################
###   Functions for Scientific Computing   ###
##############################################
 */

/* solve linear simultaneous equations
 <NUMERICAL RECIPES, the art of scientific computing>,
 William H. Press, 1986, page 31 */
double g_equation(double *a, double *b, double *x, int n, int kk)
{ int i, j, k, imax=0;
  double temp, vmax, small;
  
  small = 1.e-20;
  
  if (kk >= 0 && kk < n)
  { for (i = 0; i < n; ++i)
  { j = i * n;
    temp = (*(a+j+kk));
    *(a+j+kk) = (*(a+j+n-1));
    *(a+j+n-1) = temp;
  }
  }
  
  for (i = 0; i < n; ++i)  /* scaling each row */
  { vmax = 0.;
    for (j = 0; j < n; ++j)
    { temp = (*(a+n*i+j));
      temp = (temp >= 0.) ? temp : (-temp);
      if (temp > vmax)
        vmax = temp;
    }
    if (vmax > -small && vmax < small)
    { printf("\n\n ZERO PIVOT IN equation\n");
      vmax = small;
    }
    x[i] = 1./vmax;
    /*x array is used to store scaling factor temporarily*/
  }
  
  for (j = 0; j < n; ++j)  /* loop each colum */
  { /* evaluate upper triangle elmemts */
    for (i = 0; i < j; ++i)
    { temp = (*(a+n*i+j));
      for (k = 0; k < i; ++k)
        temp -= (*(a+n*i+k)) * (*(a+n*k+j));
      
      *(a + i*n + j) = temp;
    }
    
    vmax = 0.;
    
    for (i = j; i < n; ++i)   /*evaluate diagonal and lower angle
     element*/
    { temp = (*(a+n*i+j));
      for (k = 0; k < j; ++k)
        temp -= (*(a+n*i+k)) * (*(a+n*k+j));
      
      *(a + i*n + j) = temp;
      temp = (temp >= 0.) ? temp : (-temp);
      temp *= x[i];
      
      if (temp > vmax)   /* evaluate the merit of pivot */
      { vmax = temp;
        imax = i;
      }
    }
    
    if (j != imax)
    { for (i = 0; i < n; ++i)  /* shift row between j and imax */
    { temp = (*(a+n*imax+i));
      *(a + imax*n + i) = (*(a+n*j+i));
      *(a + j*n + i) = temp;
    }
    x[imax] = x[j];
      temp = b[imax];
      b[imax] = b[j];
      b[j] = temp;
    }
    
    if ((*(a+n*j+j)) > -small && (*(a+n*j+j)) < small)
    { printf("\n\n ZERO PIVOT IN equation\n");
      *(a + j*n + j) = small;
    }
    
    if (j != n-1)     /* complete lower trangle element */
    { temp = 1. / (*(a+n*j+j));
      for (i = j+1; i < n; ++i)
        *(a + i*n + j) *= temp;
    }
  }
  
  for (i = 0; i < n; ++i)    /* evaluate mediate solution vector */
  { temp = b[i];
    for (j = 0; j < i; ++j)
      temp -= (*(a+n*i+j)) * b[j];
    
    b[i] = temp;
  }
  
  if (kk >= 0 && kk < n)
    temp = b[n-1] / (*(a+n*n-1));
  else
  { for (i = n-1; i >=0; --i)  /* evaluate final solution vector */
  { temp = b[i];
    for (j = i+1; j < n; ++j)
      temp -= (*(a+n*i+j)) * b[j];
    
    b[i] = temp / (*(a+n*i+i));
  }
  
  for (i = 0; i < n; ++i) { /* write solution to array x */
    x[i] = b[i];
    
    temp = 0.;
  }
  }
  return (temp);
}

double pnorm(double x)
{
  return 0.5 * erfc(-x * M_SQRT1_2);
}

double pnchisq(double x, double df, double ncp) {
  double y;
  boost::math::non_central_chi_squared_distribution<double> non_central_chi_squared_distribution(df, ncp);
  
  y = cdf(non_central_chi_squared_distribution, x);
  
  return y;
}


// [[Rcpp::export]]
double ARL(int method, double g, double h, double k, double gamma, double deltamu, double deltasigma) {
  double C[M], R[M*M], N[M];
  return arl_sprt(method, g, h, k, gamma, deltamu, deltasigma, C, R, N);
}

// [[Rcpp::export]]
double ATS(int method, double g, double h, double k, double gamma, double deltamu, double deltasigma, double d) {
  double C[M], R[M*M], N[M];
  return ats_sprt(method, g, h, k, gamma, deltamu, deltasigma, d, C, R, N);
}

// [[Rcpp::export]]
double SDRL(int method, double g, double h, double k, double gamma, double deltamu, double deltasigma) {
  double C[M], R[M*M], N[M];
  return sdrl_sprt(method, g, h, k, gamma, deltamu, deltasigma, C, R, N);
}

// [[Rcpp::export]]
double SDTS(int method, double g, double h, double k, double gamma, double deltamu, double deltasigma, double d) {
  double C[M], R[M*M], N[M];
  return sdts_sprt(method, g, h, k, gamma, deltamu, deltasigma, d, C, R, N);
}

// [[Rcpp::export]]
double ASN(int method, double g, double h, double k, double gamma, double deltamu, double deltasigma) {
  double C[M], R[M*M], N[M];
  return asn_sprt(method, g, h, k, gamma, deltamu, deltasigma, C, R, N);
}


// [[Rcpp::export]]
double ANOS(int method, double g, double h, double k, double gamma, double deltamu, double deltasigma) {
  double C[M], R[M*M], N[M];
  return anos_sprt(method, g, h, k, gamma, deltamu, deltasigma, C, R, N);
}

// [[Rcpp::export]]
double AEQL(int method, double g, double h, double k, double gamma, double d, double deltamumin, double deltamumax, double deltasigmamin, double deltasigmamax){
  double C[M], R[M*M], N[M];
  return aeql_sprt(method, g, h, k, gamma, d, C, R, N, deltamumin, deltamumax, deltasigmamin, deltasigmamax);
}

// [[Rcpp::export]]
double EANOS(int method, double g, double h, double k, double gamma, double deltamumin, double deltamumax, double deltasigmamin, double deltasigmamax){
  double C[M], R[M*M], N[M];
  return eanos_sprt(method, g, h, k, gamma, C, R, N, deltamumin, deltamumax, deltasigmamin, deltasigmamax);
}

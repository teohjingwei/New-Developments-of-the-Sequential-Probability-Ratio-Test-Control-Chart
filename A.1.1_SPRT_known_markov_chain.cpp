#include <Rcpp.h>
#include <iostream>
#include <cmath>
using namespace std;
using namespace Rcpp;

#define M 800

double weight[9] = {0.0812743884, 0.1806481607, 0.2606106964, 0.3123470770, 0.3302393550, 0.3123470770, 0.2606106964, 0.1806481607, 0.0812743884};
double abscissa[9] = {-0.9681602395, -0.8360311073, -0.6133714327, -0.3242534234, 0., 0.3242534234,  0.6133714327,  0.8360311073,  0.9681602395};

/*
#################################
###   Run-length properties   ###
#################################
 */

double arl_sprt(double g, double h, double gamma, double delta, double C[], double R[], double N[])
{ double arl;
  double pnorm(double x);
  void Rmat(double g, double h, double gamma, double delta, double R[]), Qvec(double g, double h, double gamma, double delta, double C[]), Nvec(int cases, double C[], double R[], double N[]), Cvec(double g, double h, double gamma, double delta, double C[]);
  
  Rmat(g, h, gamma, delta, R);
  
  Qvec(g, h, gamma, delta, C);
  
  Nvec(2, C, R, N);
  
  Cvec(g, h, gamma, delta, C);
  
  double P0 = pnorm(g + gamma - delta);
  for(int i = 0; i < M; ++i)
  { P0 += C[i] * N[i];
  }
  
  arl = 1./(1. - P0);
  
  return arl;
}

double ats_sprt(double g, double h, double gamma, double delta, double d, double C[], double R[], double N[])
{ double ats;
  double pnorm(double x);
  void Rmat(double g, double h, double gamma, double delta, double R[]), Qvec(double g, double h, double gamma, double delta, double C[]), Nvec(int cases, double C[], double R[], double N[]), Cvec(double g, double h, double gamma, double delta, double C[]);
  
  Rmat(g, h, gamma, delta, R);
  
  Qvec(g, h, gamma, delta, C);
  
  Nvec(2, C, R, N);
  
  Cvec(g, h, gamma, delta, C);
  
  double P0 = pnorm(g + gamma - delta);
  for(int i = 0; i < M; ++i)
  { P0 += C[i] * N[i];
  }
  
  ats = d/(1. - P0);
  
  if(delta > 0.0)
    ats -= d * 0.5; 
  
  return ats;
}

double sdrl_sprt(double g, double h, double gamma, double delta, double C[], double R[], double N[])
{ double sdrl;
  double pnorm(double x);
  void Rmat(double g, double h, double gamma, double delta, double R[]), Qvec(double g, double h, double gamma, double delta, double C[]), Nvec(int cases, double C[], double R[], double N[]), Cvec(double g, double h, double gamma, double delta, double C[]);
  
  Rmat(g, h, gamma, delta, R);
  
  Qvec(g, h, gamma, delta, C);
  
  Nvec(2, C, R, N);
  
  Cvec(g, h, gamma, delta, C);
  
  double P0 = pnorm(g + gamma - delta);
  for(int i = 0; i < M; ++i)
  { P0 += C[i] * N[i];
  }
  
  sdrl = sqrt(P0/((1.0 - P0)*(1.0 - P0)));
  
  return sdrl;
}

double sdts_sprt(double g, double h, double gamma, double delta, double d, double C[], double R[], double N[])
{ double sdts;
  double pnorm(double x);
  void Rmat(double g, double h, double gamma, double delta, double R[]), Qvec(double g, double h, double gamma, double delta, double C[]), Nvec(int cases, double C[], double R[], double N[]), Cvec(double g, double h, double gamma, double delta, double C[]);
  
  Rmat(g, h, gamma, delta, R);
  
  Qvec(g, h, gamma, delta, C);
  
  Nvec(2, C, R, N);
  
  Cvec(g, h, gamma, delta, C);
  
  double P0 = pnorm(g + gamma - delta);
  for(int i = 0; i < M; ++i)
  { P0 += C[i] * N[i];
  }
  
  if(delta == 0.0){
    sdts = d * sqrt(P0/((1.0 - P0)*(1.0 - P0)));
  }else{
    sdts = d * sqrt(1.0/12.0 + P0/((1.0 - P0)*(1.0 - P0)));
  }
  
  return sdts;
}

double asn_sprt(double g, double h, double gamma, double delta, double C[], double R[], double N[])
{ 
  double pnorm(double x);
  void Rmat(double g, double h, double gamma, double delta, double R[]), Qvec(double g, double h, double gamma, double delta, double C[]), Nvec(int cases, double C[], double R[], double N[]), Cvec(double g, double h, double gamma, double delta, double C[]);
  
  Rmat(g, h, gamma, delta, R);
  
  Nvec(1, C, R, N);
  
  Cvec(g, h, gamma, delta, C);
  
  double asn = 1.0;
  for(int i = 0; i < M; ++i)
  { asn += C[i] * N[i];
  }
  
  return asn;
}

double anos_sprt(double g, double h, double gamma, double delta, double C[], double R[], double N[])
{
  return arl_sprt(g, h, gamma, delta, C, R, N) * asn_sprt(g, h, gamma, delta, C, R, N);
}

double aeql_sprt(double g, double h, double gamma, double d, double C[], double R[], double N[], double deltamin, double deltamax)
{
  double sum = 0.0;
  for(int i = 0; i < 9; ++i){
    double deltaval = 0.5 * ((deltamax - deltamin)*abscissa[i] + (deltamax + deltamin));
    sum += weight[i] * deltaval * deltaval * ats_sprt(g, h, gamma, deltaval, d, C, R, N);
  }
  sum *= 0.5;
  return sum;
}

double eanos_sprt(double g, double h, double gamma, double C[], double R[], double N[], double deltamin, double deltamax)
{
  double sum = 0.0;
  for(int i = 0; i < 9; ++i){
    double deltaval = 0.5 * ((deltamax - deltamin)*abscissa[i] + (deltamax + deltamin));
    sum += weight[i] * anos_sprt(g, h, gamma, deltaval, C, R, N);
  }
  sum *= 0.5;
  return sum;
}


/*
##############################
###   Vectors & Matrices   ###
##############################
 */

void Rmat(double g, double h, double gamma, double delta, double R[])
{ int rn, kk;
  double temp, last;
  double pnorm(double x);
  
  double r = (h - g) / M;
  double U = -0.5*r + gamma - delta;
  double start = pnorm(U);
  
  last = start;
  for(int i = 0; i < M; ++i)
  { U += r;
    temp = pnorm(U);
    R[i] = temp - last;
    last = temp;
  }
  
  last = start;
  U = -0.5*r + gamma - delta;
  for(int i = 1; i < M; ++i)
  { U -= r;
    temp = pnorm(U);
    R[i * M] = last - temp;
    last = temp;
  }
  
  for(int i = 1; i < M; ++i)
  { rn = i * M;
    for(int j = 1; j < M; ++j)
    { kk = rn + j;
      R[kk] = R[kk-M-1];
    }
  }
}

void Qvec(double g, double h, double gamma, double delta, double C[])
{ double U;
  double pnorm(double x);
  
  double r = (h - g) / M;
  
  for(int i = 0; i < M; ++i)
    { if(i == 0){
      U = gamma - 0.5*r  - delta;
    }else{
      U -= r;
    }
    
    C[i] = pnorm(U);
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

void Cvec(double g, double h, double gamma, double delta, double C[])
{
  double temp, last;
  double pnorm(double x);
  
  double r = (h - g) / M;
  double U = g + gamma - delta;
  last = pnorm(U);
  
  for(int i = 0; i < M; ++i)
  { U += r;
    temp = pnorm(U);
    C[i] = temp - last;
    last = temp;
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
  { for (i = n-1; i >=0; i--)  /* evaluate final solution vector */
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


// [[Rcpp::export]]
double ARL(double g, double h, double gamma, double delta) {
  double C[M], R[M*M], N[M];
  return arl_sprt(g, h, gamma, delta, C, R, N);
}

// [[Rcpp::export]]
double ATS(double g, double h, double gamma, double delta, double d) {
  double C[M], R[M*M], N[M];
  return ats_sprt(g, h, gamma, delta, d, C, R, N);
}

// [[Rcpp::export]]
double SDRL(double g, double h, double gamma, double delta) {
  double C[M], R[M*M], N[M];
  return sdrl_sprt(g, h, gamma, delta, C, R, N);
}

// [[Rcpp::export]]
double SDTS(double g, double h, double gamma, double delta, double d) {
  double C[M], R[M*M], N[M];
  return sdts_sprt(g, h, gamma, delta, d, C, R, N);
}

// [[Rcpp::export]]
double ASN(double g, double h, double gamma, double delta) {
  double C[M], R[M*M], N[M];
  return asn_sprt(g, h, gamma, delta, C, R, N);
}


// [[Rcpp::export]]
double ANOS(double g, double h, double gamma, double delta) {
  double C[M], R[M*M], N[M];
  return anos_sprt(g, h, gamma, delta, C, R, N);
}

// [[Rcpp::export]]
double AEQL(double g, double h, double gamma, double d, double deltamin, double deltamax){
  double C[M], R[M*M], N[M];
  return aeql_sprt(g, h, gamma, d, C, R, N, deltamin, deltamax);
}

// [[Rcpp::export]]
double EANOS(double g, double h, double gamma, double deltamin, double deltamax){
  double C[M], R[M*M], N[M];
  return eanos_sprt(g, h, gamma, C, R, N, deltamin, deltamax);
}

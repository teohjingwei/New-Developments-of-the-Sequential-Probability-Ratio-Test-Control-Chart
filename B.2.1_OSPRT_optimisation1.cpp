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
################################
###   Optimisation Builder   ###
################################
 */

double F1(int method, double g, double h, double k, double gamma, double d, double tau, double C[], double R[], double N[])
{
  return( ats_sprt(method, g, h, k, gamma, 0.0, 1.0, d, C, R, N) - tau );
}

double F1prime(int method, double g, double h, double dh, double k, double gamma, double d, double C[], double R[], double N[])
{
  return( ( ats_sprt(method, g, h + dh, k, gamma, 0.0, 1.0, d, C, R, N) - ats_sprt(method, g, h, k, gamma, 0.0, 1.0, d, C, R, N) )/dh);
}

double F2(int method, double g, double h, double k, double gamma, double n, double C[], double R[], double N[])
{
  return( asn_sprt(method, g, h, k, gamma, 0.0, 1.0, C, R, N) - n );
}

double F2prime(int method, double g, double dg, double h, double k, double gamma, double C[], double R[], double N[])
{
  return( ( asn_sprt(method, g + dg, h, k, gamma, 0.0, 1.0, C, R, N) - asn_sprt(method, g, h, k, gamma, 0.0, 1.0, C, R, N) )/dg);
}

double search_h(int method, double h, double hstep, double g, double k, double gamma, double n, double RATE, double tau, double C[], double R[], double N[])
{
  double temph, temp;
  double dh = pow(10,-6);
  double d = n/RATE;
  double last = ats_sprt(method, g, h, k, gamma, 0.0, 1.0, d, C, R, N);
  int counter = 0;
  
  while(counter < 50){
    if((last - tau) < -0.5){
      h += hstep;
    }else if((last - tau) > 0.5){
      h -= hstep;
    }else{
      break;
    }
    
    temp = ats_sprt(method, g, h, k, gamma, 0.0, 1.0, d, C, R, N);
    if(((last - tau)*(temp - tau)) < 0){
      if((last - tau) < -0.5){
        h += - hstep + hstep*(tau - last)/(temp - last);
      }else{
        h += hstep - hstep*(tau - last)/(temp - last);
      }
      break;
    }else{
      last = temp;
    }
    hstep *= 1.2;
    ++counter;
  }
  
  temph = h;
  counter = 0;
  while(counter < 20){
    h -= F1(method, g, h, k, gamma, d, tau, C, R, N)/F1prime(method, g, h, dh, k, gamma, d, C, R, N);
    if(fabs(h - temph) < 0.001){
      break;
    }else{
      temph = h;
    }
    ++counter;
  }
  return h;
}

double search_g(int method, double h, double g, double gstep, double k, double gamma, double n, double C[], double R[], double N[])
{
  double tempg, temp;
  double dg = pow(10,-6);
  double last = asn_sprt(method, g, h, k, gamma, 0.0, 1.0, C, R, N);
  int counter = 0;
  
  while(counter < 50){
    if((last - n)/n < -0.005){
      g -= gstep;
    }else if((last - n)/n > 0.005){
      g += gstep;
    }else{
      break;
    }
    
    temp = asn_sprt(method, g, h, k, gamma, 0.0, 1.0, C, R, N);
    if(((last - n)*(temp - n)) < 0){
      if((last - n)/n < -0.005){
        g -= gstep *(n - temp)/(temp - last);
      }else{
        g += gstep *(n - temp)/(temp - last);
      }
      break;
    }else{
      last = temp;
    }
    gstep *= 1.2;
    ++counter;
  }
  
  tempg = g;
  counter = 0;
  while(counter < 20){
    g -= F2(method, g, h, k, gamma, n, C, R, N)/F2prime(method, g, dg, h, k, gamma, C, R, N);
    if(fabs(g - tempg) < 0.001){
      break;
    }else{
      tempg = g;
    }
    ++counter;
  }
  return g;
}

void search_gh(int method, double &g, double &h, double k, double gamma, double n, double RATE, double tau, double C[], double R[], double N[])
{
  double inith = 15.0, initg = 1;
  int counter = 0;
  while(counter < 20){
    h = search_h(method, inith, 0.5, initg, k, gamma, n, RATE, tau, C, R, N);
    g = search_g(method, h, initg, 0.1, k, gamma, n, C, R, N);
    if(fabs(g - initg) < 0.001 && fabs(h - inith) < 0.001){
      break;
    }else{
      inith = h;
      initg = g;
    }
    ++counter;
  }
}

double ATSdgh_generator(const std::string& type, int method, double &g, double &h, double k, double gamma, double n, double RATE, double tau, double deltamumin, double deltamumax, double deltasigmamin, double deltasigmamax, double C[], double R[], double N[])
{
  double d = n/RATE;
  
  search_gh(method, g, h, k, gamma, n, RATE, tau, C, R, N);
  
  if(type == "ATS"){
    return ats_sprt(method, g, h, k, gamma, deltamumin, deltasigmamin, d, C, R, N);
  }else if(type == "AEQL"){
    return aeql_sprt(method, g, h, k, gamma, d, C, R, N, deltamumin, deltamumax, deltasigmamin, deltasigmamax);
  }
}

double search_gamma(const std::string& type, int method, double gammamin, double gammamax, double gammatol, double k, double n, double RATE, double tau, double deltamumin, double deltamumax, double deltasigmamin, double deltasigmamax, double C[], double R[], double N[]) {
  double invphi = (sqrt(5.0) - 1) / 2.0;
  double invphi2 = (3 - sqrt(5.0)) / 2.0;
  double a = gammamin;
  double b = gammamax;
  double r = (b-a);
  int nn;
  double c, d, fc, fd, h, g;
  
  nn = floor(ceil(log(gammatol / r) / log(invphi)));
  
  c = a + invphi2 * r;
  d = a + invphi * r;
  fc = ATSdgh_generator(type, method, g, h, k, c, n, RATE, tau, deltamumin, deltamumax, deltasigmamin, deltasigmamax, C, R, N);
  fd = ATSdgh_generator(type, method, g, h, k, d, n, RATE, tau, deltamumin, deltamumax, deltasigmamin, deltasigmamax, C, R, N);
  
  for(int j = 0; j < nn-1; ++j){
    if(fc < fd){
      b = d;
      d = c;
      fd = fc;
      r *= invphi;
      c = a + invphi2 * r;
      fc = ATSdgh_generator(type, method, g, h, k, c, n, RATE, tau, deltamumin, deltamumax, deltasigmamin, deltasigmamax, C, R, N);
    }else{
      a = c;
      c = d;
      fc = fd;
      r *= invphi;
      d = a + invphi * r;
      fd = ATSdgh_generator(type, method, g, h, k, d, n, RATE, tau, deltamumin, deltamumax, deltasigmamin, deltasigmamax, C, R, N);
    }
  }
  
  if(fc < fd){
    return (a+d)/2;
  }else{
    return (c+b)/2;
  }
}

double search_k(const std::string& type, int method, double kmin, double kmax, double ktol, double gamma, double n, double RATE, double tau, double deltamumin, double deltamumax, double deltasigmamin, double deltasigmamax, double C[], double R[], double N[]) {
  double invphi = (sqrt(5.0) - 1) / 2.0;
  double invphi2 = (3 - sqrt(5.0)) / 2.0;
  double a = kmin;
  double b = kmax;
  double r = (b-a);
  int nn;
  double c, d, fc, fd, h, g;
  
  nn = floor(ceil(log(ktol / r) / log(invphi)));
  
  c = a + invphi2 * r;
  d = a + invphi * r;
  fc = ATSdgh_generator(type, method, g, h, c, gamma, n, RATE, tau, deltamumin, deltamumax, deltasigmamin, deltasigmamax, C, R, N);
  fd = ATSdgh_generator(type, method, g, h, d, gamma, n, RATE, tau, deltamumin, deltamumax, deltasigmamin, deltasigmamax, C, R, N);
  
  for(int j = 0; j < nn-1; ++j){
    if(fc < fd){
      b = d;
      d = c;
      fd = fc;
      r *= invphi;
      c = a + invphi2 * r;
      fc = ATSdgh_generator(type, method, g, h, c, gamma, n, RATE, tau, deltamumin, deltamumax, deltasigmamin, deltasigmamax, C, R, N);
    }else{
      a = c;
      c = d;
      fc = fd;
      r *= invphi;
      d = a + invphi * r;
      fd = ATSdgh_generator(type, method, g, h, d, gamma, n, RATE, tau, deltamumin, deltamumax, deltasigmamin, deltasigmamax, C, R, N);
    }
  }
  
  if(fc < fd){
    return (a+d)/2;
  }else{
    return (c+b)/2;
  }
}

double search_kgamma(const std::string& type, int method, double &optk, double &optgamma, double &gg, double &hh, double kmin, double kmax, double gammamin, double gammamax, double ktol, double gammatol, double asn, double RATE, double tau, double deltamumin, double deltamumax, double deltasigmamin, double deltasigmamax, double C[], double R[], double N[])
{
  double initk = (kmin + kmax)/2;
  double initgamma = (gammamin + gammamax)/2;
  double temp, last = pow(10,7), kdiff = pow(10,7), gammadiff = pow(10,7), kmin2 = kmin, kmax2 = kmax, gammamin2 = gammamin, gammamax2 = gammamax;
  int counter = 0;
  
  
  while(counter < 20){
    if(gammadiff > gammatol){
      optgamma = search_gamma(type, method, gammamin2, gammamax2, gammatol, initk, asn, RATE, tau, deltamumin, deltamumax, deltasigmamin, deltasigmamax, C, R, N);
    }
    if(kdiff > ktol){
      optk = search_k(type, method, kmin2, kmax2, ktol, optgamma, asn, RATE, tau, deltamumin, deltamumax, deltasigmamin, deltasigmamax, C, R, N);
    }
    temp = ATSdgh_generator(type, method, gg, hh, optk, optgamma, asn, RATE, tau, deltamumin, deltamumax, deltasigmamin, deltasigmamax, C, R, N);
    if(fabs(temp - last)/last < 0.0001){
      break;
    }else{
      if(counter > 0){
        kdiff = initk - optk;
        gammadiff = initgamma - optgamma;
        kmin2 = (kdiff > 0) ? optk - kdiff : optk;
        kmax2 = (kdiff > 0) ? optk : optk - kdiff;
        gammamin2 = (gammadiff > 0) ? optgamma - gammadiff : optgamma;
        gammamax2 = (gammadiff > 0) ? optgamma : optgamma - gammadiff;
      }else{
        kmin2 = max(kmin, optk - 0.1);
        kmax2 = min(kmax, optk + 0.1);
        gammamin2 = max(gammamin, optgamma - 0.5);
        gammamax2 = min(gammamax, optgamma + 0.5);
      }
      initk = optk;
      initgamma = optgamma;
      last = temp;
      ++counter;
    }
  }
  
  return ATSdgh_generator(type, method, gg, hh, optk, optgamma, asn, RATE, tau, deltamumin, deltamumax, deltasigmamin, deltasigmamax, C, R, N);
}

double search_asn(const std::string& type, int method, double &optasn, double &gg, double &hh, double asnmin, double asnmax, double tol, double k, double gamma, double RATE, double tau, double deltamumin, double deltamumax, double deltasigmamin, double deltasigmamax, double C[], double R[], double N[]) {
  double invphi = (sqrt(5.0) - 1) / 2.0;
  double invphi2 = (3 - sqrt(5.0)) / 2.0;
  double a = asnmin;
  double b = asnmax;
  double r = (b-a);
  int nn;
  double c, d, fc, fd, h, g;
  
  nn = floor(ceil(log(tol / r) / log(invphi)));
  
  c = a + invphi2 * r;
  d = a + invphi * r;
  fc = ATSdgh_generator(type, method, g, h, k, gamma, c, RATE, tau, deltamumin, deltamumax, deltasigmamin, deltasigmamax, C, R, N);
  fd = ATSdgh_generator(type, method, g, h, k, gamma, d, RATE, tau, deltamumin, deltamumax, deltasigmamin, deltasigmamax, C, R, N);
  
  for(int j = 0; j < nn-1; ++j){
    if(fc < fd){
      b = d;
      d = c;
      fd = fc;
      r *= invphi;
      c = a + invphi2 * r;
      fc = ATSdgh_generator(type, method, g, h, k, gamma, c, RATE, tau, deltamumin, deltamumax, deltasigmamin, deltasigmamax, C, R, N);
    }else{
      a = c;
      c = d;
      fc = fd;
      r *= invphi;
      d = a + invphi * r;
      fd = ATSdgh_generator(type, method, g, h, k, gamma, d, RATE, tau, deltamumin, deltamumax, deltasigmamin, deltasigmamax, C, R, N);
    }
  }
  
  if(fc < fd){
    optasn = (a+d)/2;
  }else{
    optasn = (c+b)/2;
  }
  
  return ATSdgh_generator(type, method, gg, hh, k, gamma, optasn, RATE, tau, deltamumin, deltamumax, deltasigmamin, deltasigmamax, C, R, N);
}

double search_asnkgamma(const std::string& type, int method, double &optk, double &optgamma, double &optasn, double &gg, double &hh, double kmin, double kmax, double gammamin, double gammamax, double asnmin, double asnmax, double ktol, double gammatol, double ntol, double RATE, double tau, double deltamumin, double deltamumax, double deltasigmamin, double deltasigmamax, double C[], double R[], double N[]) {
  double invphi = (sqrt(5.0) - 1) / 2.0;
  double invphi2 = (3 - sqrt(5.0)) / 2.0;
  double a = asnmin;
  double b = asnmax;
  double r = (b-a);
  int nn;
  double c, d, fc, fd;
  
  nn = floor(ceil(log(ntol / r) / log(invphi)));
  
  c = a + invphi2 * r;
  d = a + invphi * r;
  fc = search_kgamma(type, method, optk, optgamma, gg, hh, kmin, kmax, gammamin, gammamax, ktol, gammatol, c, RATE, tau, deltamumin, deltamumax, deltasigmamin, deltasigmamax, C, R, N);
  fd = search_kgamma(type, method, optk, optgamma, gg, hh, kmin, kmax, gammamin, gammamax, ktol, gammatol, d, RATE, tau, deltamumin, deltamumax, deltasigmamin, deltasigmamax, C, R, N);
  
  for(int j = 0; j < nn-1; ++j){
    if(fc < fd){
      b = d;
      d = c;
      fd = fc;
      r *= invphi;
      c = a + invphi2 * r;
      fc = search_kgamma(type, method, optk, optgamma, gg, hh, kmin, kmax, gammamin, gammamax, ktol, gammatol, c, RATE, tau, deltamumin, deltamumax, deltasigmamin, deltasigmamax, C, R, N);
    }else{
      a = c;
      c = d;
      fc = fd;
      r *= invphi;
      d = a + invphi * r;
      fd = search_kgamma(type, method, optk, optgamma, gg, hh, kmin, kmax, gammamin, gammamax, ktol, gammatol, d, RATE, tau, deltamumin, deltamumax, deltasigmamin, deltasigmamax, C, R, N);
    }
  }
  
  if(fc < fd){
    optasn = (a+d)/2;
  }else{
    optasn = (c+b)/2;
  }
  
  return search_kgamma(type, method, optk, optgamma, gg, hh, kmin, kmax, gammamin, gammamax, ktol, gammatol, optasn, RATE, tau, deltamumin, deltamumax, deltasigmamin, deltasigmamax, C, R, N);
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
NumericVector optim_search_ats(List type, List parameter, List tol, int method, double RATE, double tau, double dmin) {
  double C[M], R[M*M], N[M];
  NumericVector vec(7);
  double optk, optgamma, optasn, g, h, AEQLs, deltamumin, deltamumax, deltasigmamin, deltasigmamax, kmin, kmax, gammamin, gammamax, asnmin, asnmax, ktol, gammatol, ntol;
  const std::string objfunc = type["type"];
  
  if(objfunc == "AEQL"){
    deltamumin = type["deltamin"];
    deltamumax = type["deltamax"];
    deltasigmamin = type["etamin"];
    deltasigmamax = type["etamax"];
    kmin = parameter["kmin"];
    kmax = parameter["kmax"];
    gammamin = parameter["gammamin"];
    gammamax = parameter["gammamax"];
    asnmin = parameter["asnmin"];
    asnmax = parameter["asnmax"];
    ktol = tol["ktol"];
    gammatol = tol["gammatol"];
    ntol = tol["asntol"];
  }else{
    deltamumin = type["delta"];
    deltasigmamin = type["eta"];
    asnmin = parameter["asnmin"];
    asnmax = parameter["asnmax"];
    ntol = tol["asntol"];
  }
  asnmin = max(dmin*RATE,asnmin);
  
  if(objfunc == "AEQL"){
    AEQLs = search_asnkgamma(objfunc, method, optk, optgamma, optasn, g, h, kmin, kmax, gammamin, gammamax, asnmin, asnmax, ktol, gammatol, ntol, RATE, tau, deltamumin, deltamumax, deltasigmamin, deltasigmamax, C, R, N);
  }else{
    optk = deltamumin/(deltasigmamin*deltasigmamin - 1.0);
    optgamma = optk*optk*deltasigmamin*deltasigmamin + 2*deltasigmamin*deltasigmamin*log(deltasigmamin)/(deltasigmamin*deltasigmamin - 1.0);
    
    AEQLs = search_asn(objfunc, method, optasn, g, h, asnmin, asnmax, ntol, optk, optgamma, RATE, tau, deltamumin, deltamumax, deltasigmamin, deltasigmamax, C, R, N);
  }
  
  vec(0) = AEQLs;
  vec(1) = optasn;
  vec(2) = optk;
  vec(3) = optgamma;
  vec(4) = optasn/RATE;
  vec(5) = g;
  vec(6) = h;
  return vec;
}




// [[Rcpp::export]]
double ATS(int method, double g, double h, double k, double gamma, double deltamu, double deltasigma, double d) {
  double C[M], R[M*M], N[M];
  return ats_sprt(method, g, h, k, gamma, deltamu, deltasigma, d, C, R, N);
}

// [[Rcpp::export]]
double ASN(int method, double g, double h, double k, double gamma, double deltamu, double deltasigma) {
  double C[M], R[M*M], N[M];
  return asn_sprt(method, g, h, k, gamma, deltamu, deltasigma, C, R, N);
}


// [[Rcpp::export]]
double AEQL(int method, double g, double h, double k, double gamma, double d, double deltamumin, double deltamumax, double deltasigmamin, double deltasigmamax){
  double C[M], R[M*M], N[M];
  return aeql_sprt(method, g, h, k, gamma, d, C, R, N, deltamumin, deltamumax, deltasigmamin, deltasigmamax);
}

// [[Rcpp::export]]
NumericVector ATSdgH_generator(const std::string& type, int method, double k, double gamma, double n, double RATE, double tau, double deltamumin, double deltamumax, double deltasigmamin, double deltasigmamax)
{
  double C[M], R[M*M], N[M];
  NumericVector vec(3);
  double g, h;
  
  vec(0) = ATSdgh_generator(type, method, g, h, k, gamma, n, RATE, tau, deltamumin, deltamumax, deltasigmamin, deltasigmamax, C, R, N);
  vec(1) = g;
  vec(2) = h;
  
  return vec;
}

/***
 
 // [[Rcpp::export]]
 double searchH(int method, double h, double hstep, double g, double k, double gamma, double n, double RATE, double tau){
 double C[M], R[M*M], N[M];
 
 return search_h(method, h, hstep, g, k, gamma, n, RATE, tau, C, R, N);
 }
 
 // [[Rcpp::export]]
 double searchG(int method, double h, double g, double gstep, double k, double gamma, double n){
 double C[M], R[M*M], N[M];
 
 return search_g(method, h, g, gstep, k, gamma, n, C, R, N);
 }
 
 // [[Rcpp::export]]
 NumericVector searchGH(int method, double k, double gamma, double n, double RATE, double tau){
 double C[M], R[M*M], N[M];
 double g, h;
 NumericVector vec(2);
 
 search_gh(method, g, h, k, gamma, n, RATE, tau, C, R, N);
 
 vec(0) = g;
 vec(1) = h;
 
 return vec;
 }
 
 
 // [[Rcpp::export]]
 NumericVector search_kgamma(const std::string& type, int method, double kmin, double kmax, double gammamin, double gammamax, double asn, double RATE, double tau, double deltamumin, double deltamumax, double deltasigmamin, double deltasigmamax){
 double C[M], R[M*M], N[M];
 double optk, optgamma, g, h;
 NumericVector vec(4);
 
 search_kgamma(type, method, optk, optgamma, g, h, kmin, kmax, gammamin, gammamax, asn, RATE, tau, deltamumin, deltamumax, deltasigmamin, deltasigmamax, C, R, N);
 
 vec(0) = optk;
 vec(1) = optgamma;
 vec(2) = g;
 vec(3) = h;
 
 return vec;
 }
 
 // [[Rcpp::export]]
 double searchK(const std::string& type, int method, double kmin, double kmax, double tol, double gamma, double n, double RATE, double tau, double deltamumin, double deltamumax, double deltasigmamin, double deltasigmamax){
 double C[M], R[M*M], N[M];
 
 return search_k(type, method, kmin, kmax, tol, gamma, n, RATE, tau, deltamumin, deltamumax, deltasigmamin, deltasigmamax, C, R, N);
 }
 
 // [[Rcpp::export]]
 double searchgamma(const std::string& type, int method, double gammamin, double gammamax, double tol, double k, double n, double RATE, double tau, double deltamumin, double deltamumax, double deltasigmamin, double deltasigmamax){
 double C[M], R[M*M], N[M];
 
 return search_gamma(type, method, gammamin, gammamax, tol, k, n, RATE, tau, deltamumin, deltamumax, deltasigmamin, deltasigmamax, C, R, N);
 }
 
 // [[Rcpp::export]]
 NumericVector ATSdgH_generator(const std::string& type, int method, double k, double gamma, double n, double RATE, double tau, double deltamumin, double deltamumax, double deltasigmamin, double deltasigmamax)
 {
 double C[M], R[M*M], N[M];
 NumericVector vec(3);
 double g, h;
 
 vec(0) = ATSdgh_generator(type, method, g, h, k, gamma, n, RATE, tau, deltamumin, deltamumax, deltasigmamin, deltasigmamax, C, R, N);
 vec(1) = g;
 vec(2) = h;
 
 return vec;
 }
 ***/
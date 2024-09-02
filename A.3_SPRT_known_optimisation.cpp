#include <Rcpp.h>
#include <iostream>
#include <cmath>
using namespace std;
using namespace Rcpp;

#define M 100

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

double eanos_sprt(double g, double h, double gamma, double B[], double R[], double N[], double deltamin, double deltamax)
{
  double sum = 0.0;
  for(int i = 0; i < 9; ++i){
    double deltaval = 0.5 * ((deltamax - deltamin)*abscissa[i] + (deltamax + deltamin));
    sum += weight[i] * anos_sprt(g, h, gamma, deltaval, B, R, N);
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
################################
###   Optimisation Builder   ###
################################
 */

double F1(double g, double h, double gamma, double d, double tau, double C[], double R[], double N[])
{
  return( ats_sprt(g, h, gamma, 0.0, d, C, R, N) - tau );
}

double F1prime(double g, double h, double dh, double gamma, double d, double C[], double R[], double N[])
{
  return( ( ats_sprt(g, h + dh, gamma, 0.0, d, C, R, N) - ats_sprt(g, h, gamma, 0.0, d, C, R, N) )/dh);
}

double F2(double g, double h, double gamma, double asn, double C[], double R[], double N[])
{
  return( asn_sprt(g, h, gamma, 0.0, C, R, N) - asn );
}

double F2prime(double g, double dg, double h, double gamma, double C[], double R[], double N[])
{
  return( ( asn_sprt(g + dg, h, gamma, 0.0, C, R, N) - asn_sprt(g, h, gamma, 0.0, C, R, N) )/dg);
}


double search_h(double h, double hstep, double g, double gamma, double d, double tau, double C[], double R[], double N[])
{
  double temph, temp;
  double dh = pow(10,-6);
  double F1(double g, double h, double gamma, double d, double tau, double C[], double R[], double N[]);
  double F1prime(double g, double h, double dh, double gamma, double d, double C[], double R[], double N[]);
  
  double last = ats_sprt(g, h, gamma, 0.0, d, C, R, N);
  int counter = 0;
  
  while(counter < 50){
    if((last - tau) < -0.5){
      h += hstep;
    }else if((last - tau) > 0.5){
      h -= hstep;
    }else{
      break;
    }
    
    temp = ats_sprt(g, h, gamma, 0.0, d, C, R, N);
    if(((last - tau)*(temp - tau)) < 0.0){
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
    h -= F1(g, h, gamma, d, tau, C, R, N)/F1prime(g, h, dh, gamma, d, C, R, N);
    if(fabs(h - temph) < 0.001){
      break;
    }else{
      temph = h;
    }
    ++counter;
  }
  return h;
}

double search_g(double g, double gstep, double h, double gamma, double n, double C[], double R[], double N[])
{
  double tempg, temp;
  double dg = pow(10,-6);
  double F2(double g, double h, double gamma, double asn, double C[], double R[], double N[]);
  double F2prime(double g, double dg, double h, double gamma, double C[], double R[], double N[]);
  double last = asn_sprt(g, h, gamma, 0.0, C, R, N);
  int counter = 0;
  
  while(counter < 50){
    if((last - n)/n < -0.005){
      g -= gstep;
    }else if((last - n)/n > 0.005){
      g += gstep;
    }else{
      break;
    }
    
    temp = asn_sprt(g, h, gamma, 0.0, C, R, N);
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
    g -= F2(g, h, gamma, n, C, R, N)/F2prime(g, dg, h, gamma, C, R, N);
    if(fabs(g - tempg) < 0.001){
      break;
    }else{
      tempg = g;
    }
    ++counter;
  }
  return g;
}

void search_gh(double &g, double &h, double gamma, double d, double tau, double n, double C[], double R[], double N[])
{
  double inith = 10.0, initg = 0.0;
  while(true){
    h = search_h(inith, 0.5, initg, gamma, d, tau, C, R, N);
    g = search_g(initg, 0.1, h, gamma, n, C, R, N);
    if(fabs(g - initg) < 0.001 && fabs(h - inith) < 0.001){
      break;
    }else{
      inith = h;
      initg = g;
    }
  }
}

double AEQLdgh_generator(const std::string& type, double &h, double &g, double gamma, double tau, double n, double RATE, double deltamin, double deltamax, double C[], double R[], double N[])
{
  double d = n/RATE;
  
  search_gh(g, h, gamma, d, tau, n, C, R, N);
  
  if(type == "AEQL"){
    return aeql_sprt(g, h, gamma, d, C, R, N, deltamin, deltamax);
  }
}

double EANOSgh_generator(const std::string& type, double &h, double &g, double gamma, double tau, double n, double deltamin, double deltamax, double C[], double R[], double N[])
{
  search_gh(g, h, gamma, 1.0, tau/n, n, C, R, N);
  
  if(type == "EANOS"){
    return eanos_sprt(g, h, gamma, C, R, N, deltamin, deltamax);
  }
}

double gss_1d(const std::string& type, double gammamin, double gammamax, double &optimgamma, double gammatol, double tau, double n, double delmin, double delmax, double B[], double R[], double N[], double &HH, double &gg) {
  double invphi = (sqrt(5.0) - 1) / 2.0;
  double invphi2 = (3 - sqrt(5.0)) / 2.0;
  double a = gammamin;
  double b = gammamax;
  double h = (b-a);
  int nn;
  double c, d, fc, fd, H, g;
  
  nn = floor(ceil(log(gammatol / h) / log(invphi)));
  
  c = a + invphi2 * h;
  d = a + invphi * h;
  
  fc = EANOSgh_generator(type, H, g, c, tau, n, delmin, delmax, B, R, N);
  fd = EANOSgh_generator(type, H, g, d, tau, n, delmin, delmax, B, R, N);
  
  
  for(int j = 0; j < nn-1; ++j){
    if(fc < fd){
      b = d;
      d = c;
      fd = fc;
      h *= invphi;
      c = a + invphi2 * h;
      fc = EANOSgh_generator(type, H, g, c, tau, n, delmin, delmax, B, R, N);
    }else{
      a = c;
      c = d;
      fc = fd;
      h *= invphi;
      d = a + invphi * h;
      fd = EANOSgh_generator(type, H, g, d, tau, n, delmin, delmax, B, R, N);
    }
  }
  
  if(fc < fd){
    optimgamma = (a+d)/2;
    return EANOSgh_generator(type, HH, gg, optimgamma, tau, n, delmin, delmax, B, R, N);
  }else{
    optimgamma = (c+b)/2;
    return EANOSgh_generator(type, HH, gg, optimgamma, tau, n, delmin, delmax, B, R, N);
  }
}

double gss_2d(const std::string& type, double gammamin, double gammamax, double asnmin, double asnmax, double &optimasn, double &optimgamma, double gammatol, double asntol, double RATE, double tau, double delmin, double delmax, double B[], double R[], double N[], double &HH, double &gg) {
  double invphi = (sqrt(5.0) - 1) / 2.0;
  double invphi2 = (3 - sqrt(5.0)) / 2.0;
  double gammawidth = (gammamax - gammamin);
  double asnwidth = (asnmax - asnmin);
  
  double gamma1, gamma2, asn1, asn2, fg1a1, fg1a2, fg2a1, fg2a2, g, H;
  
  int mm = ceil(log(gammatol / gammawidth) / log(invphi));
  int nn = ceil(log(asntol / asnwidth) / log(invphi));
  cout << "gamma iterations: " << mm << ", asn iterations: " << nn << endl;
  
  gamma1 = gammamin + invphi2 * gammawidth;
  gamma2 = gammamin + invphi * gammawidth;
  asn1 = asnmin + invphi2 * asnwidth;
  asn2 = asnmin + invphi * asnwidth;
  fg1a1 = AEQLdgh_generator(type, H, g, gamma1, tau, asn1, RATE, delmin, delmax, B, R, N);
  fg1a2 = AEQLdgh_generator(type, H, g, gamma1, tau, asn2, RATE, delmin, delmax, B, R, N);
  fg2a1 = AEQLdgh_generator(type, H, g, gamma2, tau, asn1, RATE, delmin, delmax, B, R, N);
  fg2a2 = AEQLdgh_generator(type, H, g, gamma2, tau, asn2, RATE, delmin, delmax, B, R, N);
  
  if(mm == nn){
    for(int j = 0; j < nn-1; ++j){
      if((fg1a1 < fg1a2) && (fg1a1 < fg2a1) && (fg1a1 < fg2a2)){
        gammamax = gamma2;
        asnmax = asn2;
        gamma2 = gamma1;
        asn2 = asn1;
        fg2a2 = fg1a1;
        gammawidth *= invphi;
        asnwidth *= invphi;
        gamma1 = gammamin + invphi2 * gammawidth;
        asn1 = asnmin + invphi2 * asnwidth;
        fg1a1 = AEQLdgh_generator(type, H, g, gamma1, tau, asn1, RATE, delmin, delmax, B, R, N);
        fg1a2 = AEQLdgh_generator(type, H, g, gamma1, tau, asn2, RATE, delmin, delmax, B, R, N);
        fg2a1 = AEQLdgh_generator(type, H, g, gamma2, tau, asn1, RATE, delmin, delmax, B, R, N);
      }else if((fg1a2 < fg1a1) && (fg1a2 < fg2a1) && (fg1a2 < fg2a2)){
        gammamax = gamma2;
        asnmin = asn1;
        gamma2 = gamma1;
        asn1 = asn2;
        fg2a1 = fg1a2;
        gammawidth *= invphi;
        asnwidth *= invphi;
        gamma1 = gammamin + invphi2 * gammawidth;
        asn2 = asnmin + invphi * asnwidth;
        fg1a1 = AEQLdgh_generator(type, H, g, gamma1, tau, asn1, RATE, delmin, delmax, B, R, N);
        fg1a2 = AEQLdgh_generator(type, H, g, gamma1, tau, asn2, RATE, delmin, delmax, B, R, N);
        fg2a2 = AEQLdgh_generator(type, H, g, gamma2, tau, asn2, RATE, delmin, delmax, B, R, N);
      }else if((fg2a1 < fg1a2) && (fg2a1 < fg1a1) && (fg2a1 < fg2a2)){
        gammamin = gamma1;
        asnmax = asn2;
        gamma1 = gamma2;
        asn2 = asn1;
        fg1a2 = fg2a1;
        gammawidth *= invphi;
        asnwidth *= invphi;
        gamma2 = gammamin + invphi * gammawidth;
        asn1 = asnmin + invphi2 * asnwidth;
        fg1a1 = AEQLdgh_generator(type, H, g, gamma1, tau, asn1, RATE, delmin, delmax, B, R, N);
        fg2a1 = AEQLdgh_generator(type, H, g, gamma2, tau, asn1, RATE, delmin, delmax, B, R, N);
        fg2a2 = AEQLdgh_generator(type, H, g, gamma2, tau, asn2, RATE, delmin, delmax, B, R, N);
      }else{
        gammamin = gamma1;
        asnmin = asn1;
        gamma1 = gamma2;
        asn1 = asn2;
        fg1a1 = fg2a2;
        gammawidth *= invphi;
        asnwidth *= invphi;
        gamma2 = gammamin + invphi * gammawidth;
        asn2 = asnmin + invphi * asnwidth;
        fg1a2 = AEQLdgh_generator(type, H, g, gamma1, tau, asn2, RATE, delmin, delmax, B, R, N);
        fg2a1 = AEQLdgh_generator(type, H, g, gamma2, tau, asn1, RATE, delmin, delmax, B, R, N);
        fg2a2 = AEQLdgh_generator(type, H, g, gamma2, tau, asn2, RATE, delmin, delmax, B, R, N);
      }
    }
  }else if(mm < nn){
    for(int j = 0; j < mm-1; ++j){
      if((fg1a1 < fg1a2) && (fg1a1 < fg2a1) && (fg1a1 < fg2a2)){
        gammamax = gamma2;
        asnmax = asn2;
        gamma2 = gamma1;
        asn2 = asn1;
        fg2a2 = fg1a1;
        gammawidth *= invphi;
        asnwidth *= invphi;
        gamma1 = gammamin + invphi2 * gammawidth;
        asn1 = asnmin + invphi2 * asnwidth;
        fg1a1 = AEQLdgh_generator(type, H, g, gamma1, tau, asn1, RATE, delmin, delmax, B, R, N);
        fg1a2 = AEQLdgh_generator(type, H, g, gamma1, tau, asn2, RATE, delmin, delmax, B, R, N);
        fg2a1 = AEQLdgh_generator(type, H, g, gamma2, tau, asn1, RATE, delmin, delmax, B, R, N);
      }else if((fg1a2 < fg1a1) && (fg1a2 < fg2a1) && (fg1a2 < fg2a2)){
        gammamax = gamma2;
        asnmin = asn1;
        gamma2 = gamma1;
        asn1 = asn2;
        fg2a1 = fg1a2;
        gammawidth *= invphi;
        asnwidth *= invphi;
        gamma1 = gammamin + invphi2 * gammawidth;
        asn2 = asnmin + invphi * asnwidth;
        fg1a1 = AEQLdgh_generator(type, H, g, gamma1, tau, asn1, RATE, delmin, delmax, B, R, N);
        fg1a2 = AEQLdgh_generator(type, H, g, gamma1, tau, asn2, RATE, delmin, delmax, B, R, N);
        fg2a2 = AEQLdgh_generator(type, H, g, gamma2, tau, asn2, RATE, delmin, delmax, B, R, N);
      }else if((fg2a1 < fg1a2) && (fg2a1 < fg1a1) && (fg2a1 < fg2a2)){
        gammamin = gamma1;
        asnmax = asn2;
        gamma1 = gamma2;
        asn2 = asn1;
        fg1a2 = fg2a1;
        gammawidth *= invphi;
        asnwidth *= invphi;
        gamma2 = gammamin + invphi * gammawidth;
        asn1 = asnmin + invphi2 * asnwidth;
        fg1a1 = AEQLdgh_generator(type, H, g, gamma1, tau, asn1, RATE, delmin, delmax, B, R, N);
        fg2a1 = AEQLdgh_generator(type, H, g, gamma2, tau, asn1, RATE, delmin, delmax, B, R, N);
        fg2a2 = AEQLdgh_generator(type, H, g, gamma2, tau, asn2, RATE, delmin, delmax, B, R, N);
      }else{
        gammamin = gamma1;
        asnmin = asn1;
        gamma1 = gamma2;
        asn1 = asn2;
        fg1a1 = fg2a2;
        gammawidth *= invphi;
        asnwidth *= invphi;
        gamma2 = gammamin + invphi * gammawidth;
        asn2 = asnmin + invphi * asnwidth;
        fg1a2 = AEQLdgh_generator(type, H, g, gamma1, tau, asn2, RATE, delmin, delmax, B, R, N);
        fg2a1 = AEQLdgh_generator(type, H, g, gamma2, tau, asn1, RATE, delmin, delmax, B, R, N);
        fg2a2 = AEQLdgh_generator(type, H, g, gamma2, tau, asn2, RATE, delmin, delmax, B, R, N);
      }
    }
    for(int j = 0; j < nn-mm; ++j){
      if(((fg1a1 < fg1a2) && (fg1a1 < fg2a1) && (fg1a1 < fg2a2)) || ((fg2a1 < fg1a2) && (fg2a1 < fg1a1) && (fg2a1 < fg2a2))){
        asnmax = asn2;
        asn2 = asn1;
        fg2a2 = fg2a1;
        fg1a2 = fg1a1;
        asnwidth *= invphi;
        asn1 = asnmin + invphi2 * asnwidth;
        fg1a1 = AEQLdgh_generator(type, H, g, gamma1, tau, asn1, RATE, delmin, delmax, B, R, N);
        fg2a1 = AEQLdgh_generator(type, H, g, gamma2, tau, asn1, RATE, delmin, delmax, B, R, N);
      }else{
        asnmin = asn1;
        asn1 = asn2;
        fg1a1 = fg1a2;
        fg2a1 = fg2a2;
        asnwidth *= invphi;
        asn2 = asnmin + invphi * asnwidth;
        fg1a2 = AEQLdgh_generator(type, H, g, gamma1, tau, asn2, RATE, delmin, delmax, B, R, N);
        fg2a2 = AEQLdgh_generator(type, H, g, gamma2, tau, asn2, RATE, delmin, delmax, B, R, N);
      }
    }
  }else{
    for(int j = 0; j < nn-1; ++j){
      if((fg1a1 < fg1a2) && (fg1a1 < fg2a1) && (fg1a1 < fg2a2)){
        gammamax = gamma2;
        asnmax = asn2;
        gamma2 = gamma1;
        asn2 = asn1;
        fg2a2 = fg1a1;
        gammawidth *= invphi;
        asnwidth *= invphi;
        gamma1 = gammamin + invphi2 * gammawidth;
        asn1 = asnmin + invphi2 * asnwidth;
        fg1a1 = AEQLdgh_generator(type, H, g, gamma1, tau, asn1, RATE, delmin, delmax, B, R, N);
        fg1a2 = AEQLdgh_generator(type, H, g, gamma1, tau, asn2, RATE, delmin, delmax, B, R, N);
        fg2a1 = AEQLdgh_generator(type, H, g, gamma2, tau, asn1, RATE, delmin, delmax, B, R, N);
      }else if((fg1a2 < fg1a1) && (fg1a2 < fg2a1) && (fg1a2 < fg2a2)){
        gammamax = gamma2;
        asnmin = asn1;
        gamma2 = gamma1;
        asn1 = asn2;
        fg2a1 = fg1a2;
        gammawidth *= invphi;
        asnwidth *= invphi;
        gamma1 = gammamin + invphi2 * gammawidth;
        asn2 = asnmin + invphi * asnwidth;
        fg1a1 = AEQLdgh_generator(type, H, g, gamma1, tau, asn1, RATE, delmin, delmax, B, R, N);
        fg1a2 = AEQLdgh_generator(type, H, g, gamma1, tau, asn2, RATE, delmin, delmax, B, R, N);
        fg2a2 = AEQLdgh_generator(type, H, g, gamma2, tau, asn2, RATE, delmin, delmax, B, R, N);
      }else if((fg2a1 < fg1a2) && (fg2a1 < fg1a1) && (fg2a1 < fg2a2)){
        gammamin = gamma1;
        asnmax = asn2;
        gamma1 = gamma2;
        asn2 = asn1;
        fg1a2 = fg2a1;
        gammawidth *= invphi;
        asnwidth *= invphi;
        gamma2 = gammamin + invphi * gammawidth;
        asn1 = asnmin + invphi2 * asnwidth;
        fg1a1 = AEQLdgh_generator(type, H, g, gamma1, tau, asn1, RATE, delmin, delmax, B, R, N);
        fg2a1 = AEQLdgh_generator(type, H, g, gamma2, tau, asn1, RATE, delmin, delmax, B, R, N);
        fg2a2 = AEQLdgh_generator(type, H, g, gamma2, tau, asn2, RATE, delmin, delmax, B, R, N);
      }else{
        gammamin = gamma1;
        asnmin = asn1;
        gamma1 = gamma2;
        asn1 = asn2;
        fg1a1 = fg2a2;
        gammawidth *= invphi;
        asnwidth *= invphi;
        gamma2 = gammamin + invphi * gammawidth;
        asn2 = asnmin + invphi * asnwidth;
        fg1a2 = AEQLdgh_generator(type, H, g, gamma1, tau, asn2, RATE, delmin, delmax, B, R, N);
        fg2a1 = AEQLdgh_generator(type, H, g, gamma2, tau, asn1, RATE, delmin, delmax, B, R, N);
        fg2a2 = AEQLdgh_generator(type, H, g, gamma2, tau, asn2, RATE, delmin, delmax, B, R, N);
      }
    }
    for(int j = 0; j < mm-nn; ++j){
      if(((fg1a1 < fg1a2) && (fg1a1 < fg2a1) && (fg1a1 < fg2a2)) || ((fg1a2 < fg1a1) && (fg1a2 < fg2a1) && (fg1a2 < fg2a2))){
        gammamax = gamma2;
        gamma2 = gamma1;
        fg2a1 = fg1a1;
        fg2a2 = fg1a2;
        gammawidth *= invphi;
        gamma1 = gammamin + invphi2 * gammawidth;
        fg1a1 = AEQLdgh_generator(type, H, g, gamma1, tau, asn1, RATE, delmin, delmax, B, R, N);
        fg1a2 = AEQLdgh_generator(type, H, g, gamma1, tau, asn2, RATE, delmin, delmax, B, R, N);
      }else{
        gammamin = gamma1;
        gamma1 = gamma2;
        fg1a1 = fg2a1;
        fg1a2 = fg2a2;
        gammawidth *= invphi;
        gamma2 = gammamin + invphi * gammawidth;
        fg2a1 = AEQLdgh_generator(type, H, g, gamma2, tau, asn1, RATE, delmin, delmax, B, R, N);
        fg2a2 = AEQLdgh_generator(type, H, g, gamma2, tau, asn2, RATE, delmin, delmax, B, R, N);
      }
    }
  }
  
  
  if((fg1a1 < fg1a2) && (fg1a1 < fg2a1) && (fg1a1 < fg2a2)){
    optimasn = (asnmin + asn2)/2;
    optimgamma = (gammamin + gamma2)/2;
  }else if((fg1a2 < fg1a1) && (fg1a2 < fg2a1) && (fg1a2 < fg2a2)){
    optimasn = (asnmax + asn1)/2;
    optimgamma = (gammamin + gamma2)/2;
  }else if((fg2a1 < fg1a2) && (fg2a1 < fg1a1) && (fg2a1 < fg2a2)){
    optimasn = (asnmin + asn2)/2;
    optimgamma = (gammamax + gamma1)/2;
  }else{
    optimasn = (asnmax + asn1)/2;
    optimgamma = (gammamax + gamma1)/2;
  }
  return AEQLdgh_generator(type, HH, gg, optimgamma, tau, optimasn, RATE, delmin, delmax, B, R, N);
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


// [[Rcpp::export]]
double ARL(double g, double h, double gamma, double delta) {
  double C[M], R[M*M], N[M];
  return arl_sprt(g, h, gamma, delta, C, R, N);
}

// [[Rcpp::export]]
double ASN(double g, double h, double gamma, double delta) {
  double C[M], R[M*M], N[M];
  return asn_sprt(g, h, gamma, delta, C, R, N);
}

// [[Rcpp::export]]
NumericVector optim_search_AEQL(double gammamin, double gammamax, double asnmin, double asnmax, double gammatol, double asntol, double RATE, double tau, double deltamin, double deltamax, double dmin) {
  double C[M], R[M*M], N[M];
  NumericVector vec(6);
  double optimn, optimgamma, AEQLmin, h, g;
  
  asnmin = max(dmin*RATE,asnmin);
  
  AEQLmin = gss_2d("AEQL", gammamin, gammamax, asnmin, asnmax, optimn, optimgamma, gammatol, asntol, RATE, tau, deltamin, deltamax, C, R, N, h, g);
  vec(0) = AEQLmin;
  vec(1) = optimn;
  vec(2) = optimgamma;
  vec(3) = optimn/RATE;
  vec(4) = g;
  vec(5) = h;
  return vec;
}

// [[Rcpp::export]]
NumericVector optim_search_EANOS(double gammamin, double gammamax, double gammatol, double asn, double tau, double deltamin, double deltamax) {
  double C[M], R[M*M], N[M];
  NumericVector vec(4);
  double optimgamma, EANOSmin, h, g;
  
  EANOSmin = gss_1d("EANOS", gammamin, gammamax, optimgamma, gammatol, tau, asn, deltamin, deltamax, C, R, N, h, g);
  vec(0) = EANOSmin;
  vec(1) = optimgamma;
  vec(2) = g;
  vec(3) = h;
  return vec;
}

// [[Rcpp::export]]
NumericVector SEARCHGH(double gamma, double d, double tau, double asn){
  double C[M], R[M*M], N[M], g, h;
  NumericVector vec(2);
  
  search_gh(g, h, gamma, d, tau, asn, C, R, N);
  
  vec(0) = g;
  vec(1) = h;
  
  return vec;
}


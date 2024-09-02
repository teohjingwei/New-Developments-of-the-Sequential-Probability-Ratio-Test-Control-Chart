#include <Rcpp.h>
#include <iostream>
#include <cmath>
#include <random>
using namespace std;
using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]
#ifdef _OPENMP
#include <omp.h>
#define thisThread omp_get_thread_num()
#else
#define thisThread 0
#endif

// [[Rcpp::depends(sitmo)]]
#include <sitmo.h>

// [[Rcpp::depends(BH)]]
#include <boost/random/normal_distribution.hpp>
#include <boost/math/distributions/normal.hpp>

#define threadRNG (parallel::rngs.r[thisThread])
namespace parallel {

struct RNGS {
  uint32_t nthreads, seed;
  std::vector< sitmo::prng_engine * > r;
  
  RNGS() : nthreads(1), seed(1) { 
    r.push_back(new sitmo::prng_engine(seed));
  }
  
  ~RNGS() {
    for (uint32_t i = 0; i<r.size(); ++i) delete r[i];
  }
  
  
  void setPThreads(uint32_t nt) {
#ifdef _OPENMP
    nthreads = nt;
    for (uint32_t i=r.size(); i<nthreads; ++i) {
      uint32_t thisSeed = (seed+i) % sitmo::prng_engine::max();
      r.push_back(new sitmo::prng_engine(thisSeed));
    }
#else
    Rcpp::warning("No openmp support");
#endif
  }
  
  void setPSeed(double s) {
    if ((s<=0.0) || (s>=1.0)) 
      Rcpp::stop("seed must be between 0 and 1");
    seed = s * sitmo::prng_engine::max();
    for (uint32_t i=0; i<r.size(); ++i) {
      uint32_t thisSeed = (seed+i) % sitmo::prng_engine::max();
      r[i]->seed(thisSeed);
    }
  }
  
};

RNGS rngs;

}

#define M 100

double weight9[9] = {0.0812743884, 0.1806481607, 0.2606106964, 0.3123470770, 0.3302393550, 0.3123470770, 0.2606106964, 0.1806481607, 0.0812743884};
double abscissa9[9] = {-0.9681602395, -0.8360311073, -0.6133714327, -0.3242534234, 0., 0.3242534234,  0.6133714327,  0.8360311073,  0.9681602395};

double weight15[15] = {0.2025782419, 0.1984314853, 0.1984314853, 0.1861610000, 0.1861610000, 0.1662692058, 0.1662692058, 0.1395706779, 0.1395706779, 0.1071592205, 0.1071592205, 0.0703660475, 0.0703660475, 0.0307532420, 0.0307532420};
double abscissa15[15] = {0.0, -0.2011940940, 0.2011940940, -0.3941513471, 0.3941513471, -0.5709721726, 0.5709721726, -0.7244177314, 0.7244177314, -0.8482065834, 0.8482065834, -0.9372733924, 0.9372733924, -0.9879925180, 0.9879925180};

/*
#############################################
###   Conditional Run length properties   ###
#############################################
 */

double arl_sprt(double g, double h, double gamma, double delta, double m, double V, double W, double C[], double R[], double N[])
{ 
  double arl;
  void Rmat(double g, double h, double gamma, double delta, double m, double V, double W, double R[]), Qvec(double g, double h, double gamma, double delta, double m, double V, double W, double C[]), Nvec(int cases, double C[], double R[], double N[]), Cvec(double g, double h, double gamma, double delta, double m, double V, double W, double C[]);
  double pnorm(double x);
  
  Rmat(g, h, gamma, delta, m, V, W, R);
  
  Qvec(g, h, gamma, delta, m, V, W, C);
  
  Nvec(2, C, R, N);
  
  Cvec(g, h, gamma, delta, m, V, W, C);
  
  double P0 = pnorm(V*(g + gamma) + W/sqrt(m) - delta);
  for(int i = 0; i < M; ++i){
    P0 += C[i] * N[i];
  }
  
  arl = 1.0 / (1.0 - P0);
  
  return arl;
}

double asn_sprt(double g, double h, double gamma, double delta, double m, double V, double W, double C[], double R[], double N[])
{ 
  void Rmat(double g, double h, double gamma, double delta, double m, double V, double W, double R[]), Qvec(double g, double h, double gamma, double delta, double m, double V, double W, double C[]), Nvec(int cases, double C[], double R[], double N[]), Cvec(double g, double h, double gamma, double delta, double m, double V, double W, double C[]);
  
  Rmat(g, h, gamma, delta, m, V, W, R);
  
  Qvec(g, h, gamma, delta, m, V, W, C);
  
  Nvec(1, C, R, N);
  
  Cvec(g, h, gamma, delta, m, V, W, C);
  
  double asn = 1.0;
  for(int i = 0; i < M; ++i){
    asn += C[i] * N[i];
  }
  
  return asn;
}

double anos_sprt(double g, double h, double gamma, double delta, double m, double V, double W, double C[], double R[], double N[])
{
  return arl_sprt(g, h, gamma, delta, m, V, W, C, R, N) * asn_sprt(g, h, gamma, delta, m, V, W, C, R, N);
}

double eanos_sprt(double g, double h, double gamma, double m, double deltamin, double deltamax, double V, double W, double C[], double R[], double N[])
{
  double sum = 0.0;
  for(int i = 0; i < 9; ++i){
    double deltaval = 0.5 * ((deltamax - deltamin)*abscissa9[i] + (deltamax + deltamin));
    sum += weight9[i] * anos_sprt(g, h, gamma, deltaval, m, V, W, C, R, N);
  }
  sum *= 0.5;
  return sum;
}


/*
###############################################
###   Unconditional Run length properties   ###
###############################################
 */

double aanos_sprt(double g, double h, double gamma, double delta, double m, double C[], double R[], double N[])
{
  double dgamma(double value, double a, double b), dnorm(double x, double mu, double sigma);
  
  double c = (m-1.0)/2.0;
  double b = 2.0/(m-1.0);
  double sV = sqrt((c-exp(2.0*(lgamma(c+0.5)-lgamma(c))))*b);
  double sW = 1.0;
  
  double sum = 0.0;
  double Vdiff = 1.0 + 5.0*sV - max(0.0,1.0 - 5.0*sV);
  double Vsum = 1.0 + 5.0*sV + max(0.0,1.0 - 5.0*sV);
  double Wdiff = 10.0*sW;
  for(int j = 0; j < 15; ++j){
    double Vval = 0.5 * ((Vdiff)*abscissa15[j] + Vsum);
    for(int i = 0; i < 15; ++i){
      double Wval = 0.5 * (Wdiff)*abscissa15[i] ;
      sum += weight15[j] * weight15[i] * anos_sprt(g, h, gamma, delta, m, Vval, Wval, C, R, N) * 2.0 * Vval * dgamma(Vval*Vval,c,c) * dnorm(Wval,0.0,sW);
    }
  }
  sum *= (0.25*Vdiff*Wdiff);
  return sum;
}

double aasn_sprt(double g, double h, double gamma, double delta, double m, double C[], double R[], double N[])
{
  double dgamma(double value, double a, double b), dnorm(double x, double mu, double sigma);
  
  double c = (m-1.0)/2.0;
  double b = 2.0/(m-1.0);
  double sV = sqrt((c-exp(2.0*(lgamma(c+0.5)-lgamma(c))))*b);
  double sW = 1.0;
  
  double sum = 0.0;
  double Vdiff = 1.0 + 5.0*sV - max(0.0,1.0 - 5.0*sV);
  double Vsum = 1.0 + 5.0*sV + max(0.0,1.0 - 5.0*sV);
  double Wdiff = 10.0*sW;
  for(int j = 0; j < 15; ++j){
    double Vval = 0.5 * ((Vdiff)*abscissa15[j] + Vsum);
    for(int i = 0; i < 15; ++i){
      double Wval = 0.5 * (Wdiff)*abscissa15[i] ;
      sum += weight15[j] * weight15[i] * asn_sprt(g, h, gamma, delta, m, Vval, Wval, C, R, N) * 2.0 * Vval * dgamma(Vval*Vval,c,c) * dnorm(Wval,0.0,sW);
    }
  }
  sum *= (0.25*Vdiff*Wdiff);
  return sum;
}

double eaanos_sprt(double g, double h, double gamma, double m, double deltamin, double deltamax, double C[], double R[], double N[])
{
  double dgamma(double value, double a, double b), dnorm(double x, double mu, double sigma);
  
  double c = (m-1.0)/2.0;
  double b = 2.0/(m-1.0);
  double sV = sqrt((c-exp(2.0*(lgamma(c+0.5)-lgamma(c))))*b);
  double sW = 1.0;
  
  double sum = 0.0;
  double Vdiff = 1.0 + 5.0*sV - max(0.0,1.0 - 5.0*sV);
  double Vsum = 1.0 + 5.0*sV + max(0.0,1.0 - 5.0*sV);
  double Wdiff = 10.0*sW;
  for(int j = 0; j < 15; ++j){
    double Vval = 0.5 * ((Vdiff)*abscissa15[j] + Vsum);
    for(int i = 0; i < 15; ++i){
      double Wval = 0.5 * (Wdiff)*abscissa15[i] ;
      sum += weight15[j] * weight15[i] * eanos_sprt(g, h, gamma, m, deltamin, deltamax, Vval, Wval, C, R, N) * 2.0 * Vval * dgamma(Vval*Vval,c,c) * dnorm(Wval,0.0,sW);
    }
  }
  sum *= (0.25*Vdiff*Wdiff);
  return sum;
}


/*
##############################
###   Vectors & Matrices   ###
##############################
 */

void Rmat(double g, double h, double gamma, double delta, double m, double V, double W, double R[])
{ 
  int rn, kk;
  double temp, last;
  double pnorm(double x);
  
  double r = (h - g) / M;
  double U = V*(-0.5*r + gamma) + W/sqrt(m) - delta;
  double start = pnorm(U);
  
  last = start;
  for(int i = 0; i < M; ++i){
    U += (V*r);
    temp = pnorm(U);
    R[i] = temp - last;
    last = temp;
  }
  
  last = start;
  U = V*(-0.5*r + gamma) + W/sqrt(m) - delta;
  for(int i = 1; i < M; ++i){
    U -= (V*r);
    temp = pnorm(U);
    R[i * M] = last - temp;
    last = temp;
  }
  
  for(int i = 1; i < M; ++i){
    rn = i * M;
    for(int j = 1; j < M; ++j){
      kk = rn + j;
      R[kk] = R[kk-M-1];
    }
  }
}

void Qvec(double g, double h, double gamma, double delta, double m, double V, double W, double C[])
{ 
  double U;
  double pnorm(double x);
  
  double r = (h - g) / M;
  
  for(int i = 0; i < M; ++i){
    if(i == 0){
      U = V*(gamma - 0.5 * r) + W/sqrt(m) - delta;
    }else{
      U -= (V*r);
    }
    C[i] = pnorm(U);
  }
}

void Nvec(int cases, double C[], double R[], double N[])
{ 
  double g_equation(double *a, double *b, double *x, int n, int kk);
  
  // determine matrix (I - R)  
  for(int i = 0; i < M; ++i){
    int ij = i * M - 1;
    
    for(int j = 0; j < M; ++j){
      ++ij;
      if(i == j){
        R[ij] = 1.0 - R[ij];
      }else{
        R[ij] = -R[ij];
      }
    }
    
    if(cases == 1)
      C[i] = 1.0;
  }
  
  g_equation(R, C, N, M, -1);
}

void Cvec(double g, double h, double gamma, double delta, double m, double V, double W, double C[])
{ 
  double U, temp, last;
  double pnorm(double x);
  
  double r = (h - g) / M;
  U = V*(g + gamma) + W/sqrt(m) - delta;
  last = pnorm(U);
  
  for(int i = 0; i < M; ++i)
  { U += (V*r);
    temp = pnorm(U);
    C[i] = temp - last;
    last = temp;
  }
}

/*
####################################
###   EPC Optimisation Builder   ###
####################################
 */

double F1(double g, double h, double gamma, double m, double tau, double C[], double R[], double N[])
{
  return( anos_sprt(g, h, gamma, 0.0, m, 1.0, 0.0, C, R, N) - tau );
}

double F1prime(double g, double h, double dh, double gamma, double m, double C[], double R[], double N[])
{
  return( ( anos_sprt(g, h + dh, gamma, 0.0, m, 1.0, 0.0, C, R, N) - anos_sprt(g, h, gamma, 0.0, m, 1.0, 0.0, C, R, N) )/dh);
}

double F2(double g, double h, double gamma, double m, double n, double C[], double R[], double N[])
{
  return( aasn_sprt(g, h, gamma, 0.0, m, C, R, N) - n );
}

double F2prime(double g, double dg, double h, double gamma, double m, double C[], double R[], double N[])
{
  return( ( aasn_sprt(g + dg, h, gamma, 0.0, m, C, R, N) - aasn_sprt(g, h, gamma, 0.0, m, C, R, N) )/dg);
}

double F3(double g, double h, double gamma, double m, double tau, double C[], double R[], double N[])
{
  return( aanos_sprt(g, h, gamma, 0.0, m, C, R, N) - tau );
}

double F3prime(double g, double h, double dh, double gamma, double m, double C[], double R[], double N[])
{
  return( ( aanos_sprt(g, h + dh, gamma, 0.0, m, C, R, N) - aanos_sprt(g, h, gamma, 0.0, m, C, R, N) )/dh);
}

double search_g(double g, double h, double gstep, double gamma, double m, double n, double C[], double R[], double N[])
{
  double tempg, temp;
  double dg = pow(10,-6);
  double F2(double g, double h, double gamma, double m, double n, double C[], double R[], double N[]);
  double F2prime(double g, double dg, double h, double gamma, double m, double C[], double R[], double N[]);
  double last = aasn_sprt(g, h, gamma, 0.0, m, C, R, N);
  int counter = 0;
  
  while(counter < 50){
    if((last - n)/n < -0.005){
      g -= gstep;
    }else if((last - n)/n > 0.005){
      g += gstep;
    }else{
      break;
    }
    
    temp = aasn_sprt(g, h, gamma, 0.0, m, C, R, N);
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
    g -= F2(g, h, gamma, m, n, C, R, N)/F2prime(g, dg, h, gamma, m, C, R, N);
    if(fabs(g - tempg) < 0.001){
      break;
    }else{
      tempg = g;
    }
    ++counter;
  }
  return g;
}

double search_h_known(double g, double h, double hstep, double gamma, double m, double tau, double C[], double R[], double N[])
{
  double temph, temp;
  double dh = pow(10,-6);
  double F1(double g, double h, double gamma, double m, double tau, double C[], double R[], double N[]);
  double F1prime(double g, double h, double dh, double gamma, double m, double C[], double R[], double N[]);
  double last = anos_sprt(g, h, gamma, 0.0, m, 1.0, 0.0, C, R, N);
  int counter = 0;
  
  while(counter < 50){
    if((last - tau) < -0.5){
      h += hstep;
    }else if((last - tau) > 0.5){
      h -= hstep;
    }else{
      break;
    }
    
    temp = anos_sprt(g, h, gamma, 0.0, m, 1.0, 0.0, C, R, N);
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
    h -= F1(g, h, gamma, m, tau, C, R, N)/F1prime(g, h, dh, gamma, m, C, R, N);
    if(fabs(h - temph) < 0.001){
      break;
    }else{
      temph = h;
    }
    ++counter;
  }
  return h;
}

double search_h1(double g, double h, double hstep, double gamma, double m, double tau, double C[], double R[], double N[])
{
  double temph, temp;
  double dh = pow(10,-6);
  double F3(double g, double h, double gamma, double m, double tau, double C[], double R[], double N[]);
  double F3prime(double g, double h, double dh, double gamma, double m, double C[], double R[], double N[]);
  double last = aanos_sprt(g, h, gamma, 0.0, m, C, R, N);
  int counter = 0;
  
  while(counter < 50){
    if((last - tau) < -0.5){
      h += hstep;
    }else if((last - tau) > 0.5){
      h -= hstep;
    }else{
      break;
    }
    
    temp = aanos_sprt(g, h, gamma, 0.0, m, C, R, N);
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
    h -= F3(g, h, gamma, m, tau, C, R, N)/F3prime(g, h, dh, gamma, m, C, R, N);
    if(fabs(h - temph) < 0.001){
      break;
    }else{
      temph = h;
    }
    ++counter;
  }
  return h;
}

double search_h2(double g, double h, double hstep, double gamma, double m, double tau, double p, double epsilon, int N_init, int N_final)
{
  double C[M], R[M*M], N[M];
  double h0 = search_h_known(g, h, hstep, gamma, m, tau, C, R, N);
  double rnorm(double deltamu, double deltasigma);
  double rchisq(double df);
  double adjustment = 0.0;
  double y;
  
  double h_guess = 1.1*h0;
  double gain = h0/100.0;
  
  h = h_guess;
  
  for(int i = (-N_init + 1); i < (N_final + 1); ++i){
    double v = sqrt(rchisq(m - 1.0)/(m - 1.0));
    double w = rnorm(0.0,1.0);
    
    y = (anos_sprt(g, h, gamma, 0.0, m, v, w, C, R, N) < (1.0-epsilon)*tau) ? (p - 1.0) : p;
    h = max(h0, h - gain*y);
    if(i > 0){
      adjustment += h;
    }
  }
  
  h = adjustment/static_cast<double>(N_final);
  
  return h;
}

void search_gh(const std::string& approach, double &g, double &h, double gamma, double asn, double m, double tau, double p, double epsilon, double C[], double R[], double N[])
{
  double inith = 10.0, initg = 0.0;
  int counter = 0;
  
  if(approach == "GICP"){
    while(counter < 5){
      h = search_h2(initg, inith, 0.5, gamma, m, tau, p, epsilon, 2000, max(counter*20000 + 20000, 60000));
      g = search_g(initg, h, 0.2, gamma, m, asn, C, R, N);
      if(fabs(g - initg) < 0.001 && fabs(h - inith)/inith < 0.005){
        break;
      }else{
        inith = h;
        initg = g;
      }
      ++counter;
    }
  }else if(approach == "AANOS-matching"){
    while(counter < 20){
      h = search_h1(initg, inith, 0.5, gamma, m, tau, C, R, N);
      g = search_g(initg, h, 0.2, gamma, m, asn, C, R, N);
      if(fabs(g - initg) < 0.001 && fabs(h - inith) < 0.001){
        break;
      }else{
        inith = h;
        initg = g;
      }
      ++counter;
    }
  }
}

double AANOSgh_generator(const std::string& type, const std::string& approach, double &g, double &h, double gamma, double asn, double m, double tau, double p, double epsilon, double deltamin, double deltamax, double C[], double R[], double N[])
{
  search_gh(approach, g, h, gamma, asn, m, tau, p, epsilon, C, R, N);
  
  // it is possible to build optimal SPRT chart with different objective functions in the future
  if(type == "AANOS"){
    return aanos_sprt(g, h, gamma, deltamin, m, C, R, N);
  }else if(type == "EAANOS"){
    return eaanos_sprt(g, h, gamma, m, deltamin, deltamax, C, R, N);
  }
}

double search_r_golden(const std::string& type, const std::string& approach, double gammamin, double gammamax, double &optimgamma, double tol, double asn, double m, double tau, double p, double epsilon, double deltamin, double deltamax, double B[], double R[], double N[], double &GG, double &HH) {
  double invphi = (sqrt(5.0) - 1) / 2.0;
  double invphi2 = (3 - sqrt(5.0)) / 2.0;
  double a = gammamin;
  double b = gammamax;
  double width = (b-a);
  int nn;
  double c, d, fc, fd, h, g;
  
  nn = floor(ceil(log(tol / width) / log(invphi)));
  
  c = a + invphi2 * width;
  d = a + invphi * width;
  fc = AANOSgh_generator(type, approach, g, h, c, asn, m, tau, p, epsilon, deltamin, deltamax, B, R, N);
  fd = AANOSgh_generator(type, approach, g, h, d, asn, m, tau, p, epsilon, deltamin, deltamax, B, R, N);
  
  for(int j = 0; j < nn-1; ++j){
    if(fc < fd){
      b = d;
      d = c;
      fd = fc;
      width *= invphi;
      c = a + invphi2 * width;
      fc = AANOSgh_generator(type, approach, g, h, c, asn, m, tau, p, epsilon, deltamin, deltamax, B, R, N);
    }else{
      a = c;
      c = d;
      fc = fd;
      width *= invphi;
      d = a + invphi * width;
      fd = AANOSgh_generator(type, approach, g, h, d, asn, m, tau, p, epsilon, deltamin, deltamax, B, R, N);
    }
  }
  
  if(fc < fd){
    optimgamma = (a+d)/2;
  }else{
    optimgamma = (c+b)/2;
  }
  return AANOSgh_generator(type, approach, GG, HH, optimgamma, asn, m, tau, p, epsilon, deltamin, deltamax, B, R, N);
}


/*
###########################################
###   Standard distribution functions   ###
###########################################
 */

double dgamma(double value, double a, double b)
{
  return pow(value, (a-1.0))*exp(a*log(b) - b*value - lgamma(a));
}

double dnorm(double x, double mu, double sigma)
{
  static const double inv_sqrt_2pi = 1.0/sqrt(2*M_PI);
  double y = (x - mu) / sigma;
  
  return inv_sqrt_2pi / sigma * exp(-0.5 * y * y);
}

double pnorm(double x)
{
  return 0.5 * erfc(-x * M_SQRT1_2);
}

double rnorm(double deltamu, double deltasigma)
{
  double x;
  boost::random::normal_distribution<double> normal;
  
  sitmo::prng_engine *rng = threadRNG;
  x = normal(*rng);
  x = x*deltasigma + deltamu;
  
  return x;
}

double runif(double min, double max)
{
  std::uniform_real_distribution<> unif(min,max);
  
  sitmo::prng_engine *rng = threadRNG;
  
  return unif(*rng);
}

double rchisq(double df)
{
  chi_squared_distribution<double> chisq(df);
  
  sitmo::prng_engine *rng = threadRNG;
  return chisq(*rng);
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


// [[Rcpp::export]]
NumericVector optim_search_anos(List type, List approach, double gammamin, double gammamax, double gammatol, double asn, double m, double tau) {
  double C[M], R[M*M], N[M];
  NumericVector vec(4);
  double optimr, g, h, EAANOSs, p, epsilon, deltamin, deltamax;
  const std::string method = approach["approach"];
  const std::string objfunc = type["type"];
  if(method == "GICP"){
    p = approach["p"];
    epsilon = approach["epsilon"];
  }else{
    p = 0.0;
    epsilon = 0.0;
  }
  if(objfunc == "EAANOS"){
    deltamin = type["deltamin"];
    deltamax = type["deltamax"];
  }else{
    deltamin = type["delta"];
  }

  EAANOSs = search_r_golden(objfunc, method, gammamin, gammamax, optimr, gammatol, asn, m, tau, p, epsilon, deltamin, deltamax, C, R, N, g, h);
  vec(0) = EAANOSs;
  vec(1) = optimr;
  vec(2) = g;
  vec(3) = h;
  return vec;
} 

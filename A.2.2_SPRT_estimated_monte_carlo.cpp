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

double ats_sprt(double g, double h, double gamma, double delta, double d, double m, double V, double W, double C[], double R[], double N[])
{ 
  double ats;
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
  
  ats = d / (1.0 - P0);
  
  if(delta > 0.0){
    ats -= d * 0.5; 
  }
  
  return ats;
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
###############################
###   Simulation Programs   ###
###############################
 */

int RLgenerator(double g, double h, double gamma, double delta, double m)
{
  int RL = 0;
  double rnorm(double deltamu, double deltasigma), rchisq(double df), runif(double min, double max);
  
  double xbar = rnorm(0.0,1.0/sqrt(m));
  double s = sqrt(rchisq(m - 1.0)/(m - 1.0));
  
  while(true){
    ++RL;
    double U = 0.0;
    while(true){
      double x = rnorm(delta, 1.0);
      U += (x - xbar)/s - gamma;
      if((U < g) || (U > h)){ break; }
    }
    if(U > h) { break; }
  }
  
  return RL;
}

int SNgenerator(double g, double h, double gamma, double delta, double m)
{
  int SN = 0;
  double rnorm(double deltamu, double deltasigma), rchisq(double df), runif(double min, double max);
  
  double xbar = rnorm(0.0,1.0/sqrt(m));
  double s = sqrt(rchisq(m - 1.0)/(m - 1.0));
  
  double U = 0.0;
  while(true){
    ++SN;
    double x = rnorm(delta, 1.0);
    U += (x - xbar)/s - gamma;
    if((U < g) || (U > h)){ break; }
  }
  
  return SN;
}

double TSgenerator(double g, double h, double gamma, double delta, double d, double m)
{
  double TS = 0.0;
  double rnorm(double deltamu, double deltasigma), rchisq(double df), runif(double min, double max);
  
  double xbar = rnorm(0.0,1.0/sqrt(m));
  double s = sqrt(rchisq(m - 1.0)/(m - 1.0));
  
  while(true){
    TS += 1.0;
    double U = 0.0;
    while(true){
      double x = rnorm(delta, 1.0);
      U += (x - xbar)/s - gamma;
      if((U < g) || (U > h)){ break; }
    }
    if(U > h) { break; }
  }
  
  TS -= (delta == 0.0) ? 0.0 : runif(0.0,1.0);
  
  return d * TS;
}

int NOSgenerator(double g, double h, double gamma, double delta, double m)
{
  int NOS = 0;
  double rnorm(double deltamu, double deltasigma), rchisq(double df), runif(double min, double max);
  
  double xbar = rnorm(0.0,1.0/sqrt(m));
  double s = sqrt(rchisq(m - 1.0)/(m - 1.0));
  
  while(true){
    double U = 0.0;
    while(true){
      ++NOS;
      double x = rnorm(delta, 1.0);
      U += (x - xbar)/s - gamma;
      if((U < g) || (U > h)){ break; }
    }
    if(U > h) { break; }
  }
  
  return NOS;
}

inline void CATS_EP(double g, double h, double gamma, double d, double m, double tau, double epsilon, int nsim, int temp[])
{
  double rnorm(double deltamu, double deltasigma), rchisq(double df), runif(double min, double max);
  
#ifdef _OPENMP
  size_t nt = parallel::rngs.nthreads;
#pragma omp parallel for num_threads(nt)
#endif
  for(int j = 0; j < nsim; ++j){
    double C[M], R[M*M], N[M];
    double w = rnorm(0.0,1.0);
    double v = sqrt(rchisq(m - 1.0)/(m - 1.0));
    
    double cats = ats_sprt(g, h, gamma, 0.0, d, m, v, w, C, R, N);
    if(cats < tau*(1.0-epsilon)){
      temp[j] = 0;
    }else{
      temp[j] = 1;
    }
  }
}

inline void CANOS_EP(double g, double h, double gamma, double m, double tau, double epsilon, int nsim, int temp[])
{
  double rnorm(double deltamu, double deltasigma), rchisq(double df), runif(double min, double max);
  
#ifdef _OPENMP
  size_t nt = parallel::rngs.nthreads;
#pragma omp parallel for num_threads(nt)
#endif
  for(int j = 0; j < nsim; ++j){
    double C[M], R[M*M], N[M];
    double w = rnorm(0.0,1.0);
    double v = sqrt(rchisq(m - 1.0)/(m - 1.0));
    
    double canos = anos_sprt(g, h, gamma, 0.0, m, v, w, C, R, N);
    if(canos < tau*(1.0-epsilon)){
      temp[j] = 0;
    }else{
      temp[j] = 1;
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

/*
###########################################
###   Standard distribution functions   ###
###########################################
 */

// generate random sample from a normal distribution
double rnorm(double deltamu, double deltasigma)
{
  double x;
  boost::random::normal_distribution<double> normal;
  
  sitmo::prng_engine *rng = threadRNG;
  x = normal(*rng);
  x = x*deltasigma + deltamu;
  
  return x;
}

// generate random sample from a uniform distribution
double runif(double min, double max)
{
  std::uniform_real_distribution<> unif(min,max);
  
  sitmo::prng_engine *rng = threadRNG;
  
  return unif(*rng);
}

// generate random sample from a chi-squared distribution
double rchisq(double df)
{
  chi_squared_distribution<double> chisq(df);
  
  sitmo::prng_engine *rng = threadRNG;
  return chisq(*rng);
}

double pnorm(double x)
{
  return 0.5 * erfc(-x * M_SQRT1_2);
}



// [[Rcpp::export]]
bool hasOMP() {
  bool ans;
#ifdef _OPENMP
  ans = true;
#else
  ans = false;
#endif
  return ans;
}

// [[Rcpp::export]]
void setOMPThreads(uint32_t nthreads) {
  parallel::rngs.setPThreads(nthreads);
}

// [[Rcpp::export]]
void setSITMOSeeds(double seed) {
  parallel::rngs.setPSeed(seed);
}

// [[Rcpp::export]]
uint32_t getOMPThreads() {
  return parallel::rngs.nthreads;
}

// [[Rcpp::export]]
double AARL_Simulation(double g, double h, double gamma, double delta, double m, int nsim){
  double AARL = 0.0;
  
#ifdef _OPENMP
  size_t nt = parallel::rngs.nthreads;
#pragma omp parallel for num_threads(nt)
#endif
  for(int i = 0; i < nsim; ++i){
    int RL = RLgenerator(g, h, gamma, delta, m);
    AARL += static_cast<double>(RL);
  }
  AARL /= static_cast<double>(nsim);
  
  return AARL;
}

// [[Rcpp::export]]
double ASDRL_Simulation(double g, double h, double gamma, double delta, double m, int nsim){
  double RLvec[nsim];
  double AARL = 0.0, ASDRL = 0.0;
  
#ifdef _OPENMP
  size_t nt = parallel::rngs.nthreads;
#pragma omp parallel for num_threads(nt)
#endif
  for(int i = 0; i < nsim; ++i){
    int RL = RLgenerator(g, h, gamma, m, delta);
    RLvec[i] = static_cast<double>(RL);
    AARL += RLvec[i];
  }
  AARL /= static_cast<double>(nsim);
  
  for(int i = 0; i < nsim; ++i){
    ASDRL += (RLvec[i] - AARL)*(RLvec[i] - AARL);
  }
  ASDRL /= static_cast<double>(nsim - 1);
  ASDRL = sqrt(ASDRL);
  
  return ASDRL;
}

// [[Rcpp::export]]
double AATS_Simulation(double g, double h, double gamma, double delta, double d, double m, int nsim){
  double AATS = 0.0;
  
#ifdef _OPENMP
  size_t nt = parallel::rngs.nthreads;
#pragma omp parallel for num_threads(nt)
#endif
  for(int i = 0; i < nsim; ++i){
    double TS = TSgenerator(g, h, gamma, delta, d, m);
    AATS += TS;
  }
  AATS /= static_cast<double>(nsim);
  
  return AATS;
}

// [[Rcpp::export]]
double ASDTS_Simulation(double g, double h, double gamma, double delta, double d, double m, int nsim){
  double TSvec[nsim];
  double AATS = 0.0, ASDTS = 0.0;
  
#ifdef _OPENMP
  size_t nt = parallel::rngs.nthreads;
#pragma omp parallel for num_threads(nt)
#endif
  for(int i = 0; i < nsim; ++i){
    double TS = TSgenerator(g, h, gamma, delta, d, m);
    TSvec[i] = TS;
    AATS += TS;
  }
  AATS /= static_cast<double>(nsim);
  
  for(int i = 0; i < nsim; ++i){
    ASDTS += (TSvec[i] - AATS)*(TSvec[i] - AATS);
  }
  ASDTS /= static_cast<double>(nsim - 1);
  ASDTS = sqrt(ASDTS);
  
  return ASDTS;
}

// [[Rcpp::export]]
double AASN_Simulation(double g, double h, double gamma, double delta, double m, int nsim){
  double AASN = 0.0;
  
#ifdef _OPENMP
  size_t nt = parallel::rngs.nthreads;
#pragma omp parallel for num_threads(nt)
#endif
  for(int i = 0; i < nsim; ++i){
    int SN = SNgenerator(g, h, gamma, delta, m);
    AASN += static_cast<double>(SN);
  }
  AASN /= static_cast<double>(nsim);
  
  return AASN;
}

// [[Rcpp::export]]
double AANOS_Simulation(double g, double h, double gamma, double delta, double m, int nsim){
  double AANOS = 0.0;
  
  for(int i = 0; i < nsim; ++i){
    int NOS = NOSgenerator(g, h, gamma, delta, m);
    AANOS += static_cast<double>(NOS);
  }
  AANOS /= static_cast<double>(nsim);
  
  return AANOS;
}

// [[Rcpp::export]]
double CATSEP(double g, double h, double gamma, double d, double m, double tau, double epsilon, int nsim){
  
  int temp[nsim];
  int count = 0;
  
  CATS_EP(g,h,gamma,d,m,tau,epsilon,nsim,temp);
  
  for(int j = 0; j < nsim; ++j){
    count += temp[j];
  }
  return static_cast<double>(count)/static_cast<double>(nsim);
}

// [[Rcpp::export]]
double CANOSEP(double g, double h, double gamma, double m, double tau, double epsilon, int nsim){
  
  int temp[nsim];
  int count = 0;
  
  CANOS_EP(g,h,gamma,m,tau,epsilon,nsim,temp);
  
  for(int j = 0; j < nsim; ++j){
    count += temp[j];
  }
  return static_cast<double>(count)/static_cast<double>(nsim);
}

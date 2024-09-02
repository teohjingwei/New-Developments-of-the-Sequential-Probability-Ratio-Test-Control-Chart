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

int RLgenerator(double g, double h, double k, double gamma, double deltamu, double deltasigma){
  int RL = 0;
  double rnorm(double deltamu, double deltasigma);
  
  while(true){
    ++RL;
    double U = 0.0;
    while(true){
      double x = rnorm(deltamu, deltasigma);
      U += (x + k)*(x + k) - gamma;
      if((U < g) || (U > h)){ break; }
    }
    if(U > h) { break; }
  }
  return RL;
}

int SNgenerator(double g, double h, double k, double gamma, double deltamu, double deltasigma)
{
  int SN = 0;
  double rnorm(double deltamu, double deltasigma);
  
  double U = 0.0;
  while(true){
    ++SN;
    double x = rnorm(deltamu, deltasigma);
    U += (x + k)*(x + k) - gamma;
    if((U < g) || (U > h)){ break; }
  }
  
  return SN;
}

int NOSgenerator(double g, double h, double k, double gamma, double deltamu, double deltasigma)
{
  int NOS = 0;
  double rnorm(double deltamu, double deltasigma);
  
  while(true){
    double U = 0.0;
    while(true){
      ++NOS;
      double x = rnorm(deltamu, deltasigma);
      U += (x + k)*(x + k) - gamma;
      if((U < g) || (U > h)){ break; }
    }
    if(U > h) { break; }
  }
  
  return NOS;
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
double ARL_Simulation(double g, double h, double k, double gamma, double deltamu, double deltasigma, int nsim){
  double ARL = 0.0;
  
#ifdef _OPENMP
  size_t nt = parallel::rngs.nthreads;
#pragma omp parallel for num_threads(nt)
#endif
  for(int i = 0; i < nsim; ++i){
    ARL += static_cast<double>(RLgenerator(g, h, k, gamma, deltamu, deltasigma));
  }
  ARL /= static_cast<double>(nsim);
  
  return ARL;
}

// [[Rcpp::export]]
double SDRL_Simulation(double g, double h, double k, double gamma, double deltamu, double deltasigma, int nsim){
  double RLvec[nsim];
  double ARL = 0.0, SDRL = 0.0;
  
#ifdef _OPENMP
  size_t nt = parallel::rngs.nthreads;
#pragma omp parallel for num_threads(nt)
#endif
  for(int i = 0; i < nsim; ++i){
    int RL = RLgenerator(g, h, k, gamma, deltamu, deltasigma);
    RLvec[i] = static_cast<double>(RL);
    ARL += RLvec[i];
  }
  ARL /= static_cast<double>(nsim);
  
  for(int i = 0; i < nsim; ++i){
    SDRL += (RLvec[i] - ARL)*(RLvec[i] - ARL);
  }
  SDRL /= static_cast<double>(nsim - 1);
  SDRL = sqrt(SDRL);
  
  return SDRL;
}

// [[Rcpp::export]]
double ATS_Simulation(double g, double h, double k, double gamma, double deltamu, double deltasigma, double d, int nsim){
  double ATS = 0.0;
  
  if((deltamu == 0.0) && (deltasigma == 1.0)){
#ifdef _OPENMP
    size_t nt = parallel::rngs.nthreads;
#pragma omp parallel for num_threads(nt)
#endif
    for(int i = 0; i < nsim; ++i){
      double TS = static_cast<double>(RLgenerator(g, h, k, gamma, deltamu, deltasigma));
      ATS += TS;
    }
    ATS /= static_cast<double>(nsim);
  }else{
#ifdef _OPENMP
    size_t nt = parallel::rngs.nthreads;
#pragma omp parallel for num_threads(nt)
#endif
    for(int i = 0; i < nsim; ++i){
      double TS = static_cast<double>(RLgenerator(g, h, k, gamma, deltamu, deltasigma)) - runif(0.0,1.0);
      ATS += TS;
    }
    ATS /= static_cast<double>(nsim);
  }
  
  return d*ATS;
}

// [[Rcpp::export]]
double SDTS_Simulation(double g, double h, double k, double gamma, double deltamu, double deltasigma, double d, int nsim){
  double TSvec[nsim];
  double ATS = 0.0, SDTS = 0.0;
  
  if((deltamu == 0.0) && (deltasigma == 1.0)){
#ifdef _OPENMP
    size_t nt = parallel::rngs.nthreads;
#pragma omp parallel for num_threads(nt)
#endif
    for(int i = 0; i < nsim; ++i){
      double TS = static_cast<double>(RLgenerator(g, h, k, gamma, deltamu, deltasigma));
      TSvec[i] = TS;
      ATS += TS;
    }
    ATS /= static_cast<double>(nsim);
  }else{
#ifdef _OPENMP
    size_t nt = parallel::rngs.nthreads;
#pragma omp parallel for num_threads(nt)
#endif
    for(int i = 0; i < nsim; ++i){
      double TS = static_cast<double>(RLgenerator(g, h, k, gamma, deltamu, deltasigma)) - runif(0.0,1.0);
      TSvec[i] = TS;
      ATS += TS;
    }
    ATS /= static_cast<double>(nsim);
  }
  
  for(int i = 0; i < nsim; ++i){
    SDTS += (TSvec[i] - ATS)*(TSvec[i] - ATS);
  }
  SDTS /= static_cast<double>(nsim - 1);
  SDTS = sqrt(SDTS);
  
  return d*SDTS;
}

// [[Rcpp::export]]
double ASN_Simulation(double g, double h, double k, double gamma, double deltamu, double deltasigma, int nsim){
  double ASN = 0.0;

#ifdef _OPENMP
  size_t nt = parallel::rngs.nthreads;
#pragma omp parallel for num_threads(nt)
#endif
  for(int i = 0; i < nsim; ++i){
    int SN = SNgenerator(g, h, k, gamma, deltamu, deltasigma);
    ASN += static_cast<double>(SN);
  }
  ASN /= static_cast<double>(nsim);
  
  return ASN;
}

// [[Rcpp::export]]
double ANOS_Simulation(double g, double h, double k, double gamma, double deltamu, double deltasigma, int nsim){
  double ANOS = 0.0;
  
#ifdef _OPENMP
  size_t nt = parallel::rngs.nthreads;
#pragma omp parallel for num_threads(nt)
#endif
  for(int i = 0; i < nsim; ++i){
    int NOS = NOSgenerator(g, h, k, gamma, deltamu, deltasigma);
    ANOS += static_cast<double>(NOS);
  }
  ANOS /= static_cast<double>(nsim);
  
  return ANOS;
}

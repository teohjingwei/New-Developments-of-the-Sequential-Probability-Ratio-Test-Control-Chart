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


// calculate the in-control standard deviation
double sigma(const std::string distr, double param1, double param2){
  if(distr == "weibull"){
    return param2 * sqrt(tgamma(1.0 + 2.0/param1) - pow(tgamma(1.0 + 1.0/param1),2));
  }else if(distr == "lognormal"){
    return sqrt(exp(2*param2 + param1*param1)*(exp(param1*param1) - 1.0));
  }else if(distr == "gamma"){
    return sqrt(param1) * param2;
  }
}

int RLgenerator(const std::string distr, double param1, double param2, double g, double h, double gamma, double delta, double xbar, double s, double sigma1){
  int RL = 0;
  double U, x;
  double rdistr(const std::string distr, double param1, double param2);
  
  while(true){
    ++RL;
    U = 0.0;
    while(true){
      x = rdistr(distr, param1, param2) + delta*sigma1;
      U += (x - xbar)/s - gamma;
      if((U < g) || (U > h)){ break; }
    }
    if(U > h) { break; }
  }
  return RL;
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

double rdistr(const std::string distr, double param1, double param2)
{
  double x;
  
  if(distr == "gamma"){
    x = R::rgamma(param1, param2);
  }else if(distr == "lognormal"){
    x = exp(rnorm(param2, param1));
  }else if(distr == "weibull"){
    x = runif(0.0, 1.0);
    x = param2 * pow(-log(x),1.0/param1);
  }
  
  return x;
}

void xbarsgenerator(const std::string distr, double param1, double param2, double &xbar, double &s, int m){
  double PhaseIsample, PhaseIvec[m];
  
  xbar = 0.0;
  s = 0.0;
  for(int j = 0; j < m; ++j){
    PhaseIsample = rdistr(distr, param1, param2);
    PhaseIvec[j] = PhaseIsample;
    xbar += (PhaseIsample - xbar)/(static_cast<double>(j+1));
  }
  for(int j = 0; j < m; ++j){
    s += (PhaseIvec[j] - xbar)*(PhaseIvec[j] - xbar);
  }
  s = sqrt(s/static_cast<double>(m - 1));
}

double CATS0generator(const std::string distr, double param1, double param2, double g, double h, double gamma, double delta, double d, double xbar, double s, double sigma1){
  double TS, CATS;
  for(int i = 0; i < 1000; ++i){
    TS = static_cast<double>(RLgenerator(distr, param1, param2, g, h, gamma, delta, xbar, s, sigma1));
    CATS += (TS - CATS)/(static_cast<double>(i+1));
  }
  return d*CATS;
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
double AATS_Simulation(List distribution, double g, double h, double gamma, double delta, double d, int m, int nsim){
  double param1, param2, sigma1;
  double ATS = 0.0;
  double PhaseIvec[m], PhaseIsample;
  const std::string distr = distribution["distribution"];
  
  if(distr == "gamma"){
    param1 = distribution["shape"];
    param2 = distribution["scale"];
  }else if(distr == "lognormal"){
    param1 = distribution["sigma"];
    param2 = distribution["mu"];
  }else if(distr == "weibull"){
    param1 = distribution["shape"];
    param2 = distribution["scale"];
  }

  sigma1 = sigma(distr, param1, param2);
  
  if(delta == 0.0){
#ifdef _OPENMP
    size_t nt = parallel::rngs.nthreads;
#pragma omp parallel for num_threads(nt)
#endif
    for(int i = 0; i < nsim; ++i){
      double xbar, s;
      xbarsgenerator(distr, param1, param2, xbar, s, m);
      double TS = static_cast<double>(RLgenerator(distr, param1, param2, g, h, gamma, delta, xbar, s, sigma1));
      ATS += TS;
    }
  }else{
#ifdef _OPENMP
    size_t nt = parallel::rngs.nthreads;
#pragma omp parallel for num_threads(nt)
#endif
    for(int i = 0; i < nsim; ++i){
      double xbar, s;
      xbarsgenerator(distr, param1, param2, xbar, s, m);
      double TS = static_cast<double>(RLgenerator(distr, param1, param2, g, h, gamma, delta, xbar, s, sigma1)) - runif(0.0,1.0);
      ATS += TS;
    }
  }
  ATS /= static_cast<double>(nsim);
  
  return d*ATS;
}

// [[Rcpp::export]]
double ASDTS_Simulation(List distribution, double g, double h, double gamma, double delta, double d, int m, int nsim){
  double TSvec[nsim];
  double ATS = 0.0, SDTS = 0.0;
  double param1, param2, mu1, sigma1;
  const std::string distr = distribution["distribution"];
  
  if(distr == "gamma"){
    param1 = distribution["shape"];
    param2 = distribution["scale"];
  }else if(distr == "lognormal"){
    param1 = distribution["sigma"];
    param2 = distribution["mu"];
  }else if(distr == "weibull"){
    param1 = distribution["shape"];
    param2 = distribution["scale"];
  }
  
  sigma1 = sigma(distr, param1, param2);
  
  
  if(delta == 0.0){
#ifdef _OPENMP
    size_t nt = parallel::rngs.nthreads;
#pragma omp parallel for num_threads(nt)
#endif
    for(int i = 0; i < nsim; ++i){
      double xbar, s;
      xbarsgenerator(distr, param1, param2, xbar, s, m);
      double TS = static_cast<double>(RLgenerator(distr, param1, param2, g, h, gamma, delta, xbar, s, sigma1));
      TSvec[i] = TS;
      ATS += TS;
    }
  }else{
#ifdef _OPENMP
    size_t nt = parallel::rngs.nthreads;
#pragma omp parallel for num_threads(nt)
#endif
    for(int i = 0; i < nsim; ++i){
      double xbar, s;
      xbarsgenerator(distr, param1, param2, xbar, s, m);
      double TS = static_cast<double>(RLgenerator(distr, param1, param2, g, h, gamma, delta, xbar, s, sigma1)) - runif(0.0,1.0);
      TSvec[i] = TS;
      ATS += TS;
    }
  }
  ATS /= static_cast<double>(nsim);
  
  for(int i = 0; i < nsim; ++i){
    SDTS += (TSvec[i] - ATS)*(TSvec[i] - ATS);
  }
  SDTS /= static_cast<double>(nsim - 1);
  SDTS = sqrt(SDTS);
  
  return d*SDTS;
}

// [[Rcpp::export]]
double CATSEP_Simulation(List distribution, double g, double h, double gamma, double d, int m, double tau, int nsim){
  double param1, param2, sigma1;
  int count = 0;
  const std::string distr = distribution["distribution"];
  
  if(distr == "gamma"){
    param1 = distribution["shape"];
    param2 = distribution["scale"];
  }else if(distr == "lognormal"){
    param1 = distribution["sigma"];
    param2 = distribution["mu"];
  }else if(distr == "weibull"){
    param1 = distribution["shape"];
    param2 = distribution["scale"];
  }
  
  sigma1 = sigma(distr, param1, param2);
  
#ifdef _OPENMP
  size_t nt = parallel::rngs.nthreads;
#pragma omp parallel for num_threads(nt)
#endif
  for(int i = 0; i < nsim; ++i){
    double xbar, s;
    xbarsgenerator(distr, param1, param2, xbar, s, m);
    if(CATS0generator(distr, param1, param2, g, h, gamma, 0.0, d, xbar, s, sigma1) > tau) count += 1;
  }
  return static_cast<double>(count)/static_cast<double>(nsim);
}

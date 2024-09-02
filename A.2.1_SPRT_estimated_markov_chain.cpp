#include <Rcpp.h>
#include <iostream>
#include <cmath>
#include <random>
using namespace std;
using namespace Rcpp;

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

double sdrl_sprt(double g, double h, double gamma, double delta, double m, double V, double W, double C[], double R[], double N[])
{ 
  double sdrl;
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
  
  sdrl = sqrt(P0 / ((1.0 - P0)*(1.0 - P0)));
  
  return sdrl;
}

double sdts_sprt(double g, double h, double gamma, double delta, double d, double m, double V, double W, double C[], double R[], double N[])
{ 
  double sdts;
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
  
  if(delta == 0.0){
    sdts = d * sqrt(P0 / ((1.0 - P0)*(1.0 - P0)));
  }else{
    sdts = d * sqrt(1.0/12.0 + P0 / ((1.0 - P0)*(1.0 - P0)));
  } 
  
  return sdts;
}

double oc_sprt(double g, double h, double gamma, double delta, double m, double V, double W, double C[], double R[], double N[])
{
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
  
  return P0;
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

double aeql_sprt(double g, double h, double gamma, double d, double m, double deltamin, double deltamax, double V, double W, double C[], double R[], double N[])
{
  double sum = 0.0;
  for(int i = 0; i < 9; ++i){
    double deltaval = 0.5 * ((deltamax - deltamin)*abscissa9[i] + (deltamax + deltamin));
    sum += weight9[i] * deltaval * deltaval * ats_sprt(g, h, gamma, deltaval, d, m, V, W, C, R, N);
  }
  sum *= 0.5;
  return sum;
}


/*
###############################################
###   Unconditional Run length properties   ###
###############################################
 */

double aarl_sprt(double g, double h, double gamma, double delta, double m, double C[], double R[], double N[])
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
      sum += weight15[j] * weight15[i] * arl_sprt(g, h, gamma, delta, m, Vval, Wval, C, R, N) * 2.0 * Vval * dgamma(Vval*Vval,c,c) * dnorm(Wval,0.0,sW);
    }
  }
  sum *= (0.25*Vdiff*Wdiff);
  return sum;
}

double aats_sprt(double g, double h, double gamma, double delta, double d, double m, double C[], double R[], double N[])
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
      sum += weight15[j] * weight15[i] * ats_sprt(g, h, gamma, delta, d, m, Vval, Wval, C, R, N) * 2.0 * Vval * dgamma(Vval*Vval,c,c) * dnorm(Wval,0.0,sW);
    }
  }
  sum *= (0.25*Vdiff*Wdiff);
  return sum;
}

double earl2_sprt(double g, double h, double gamma, double delta, double m, double C[], double R[], double N[])
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
      double temp = arl_sprt(g, h, gamma, delta, m, Vval, Wval, C, R, N);
      sum += (weight15[j] * weight15[i] * temp * 2.0 * Vval * dgamma(Vval*Vval,c,c) * dnorm(Wval,0.0,sW)) * temp;
    }
  }
  sum *= (0.25*Vdiff*Wdiff);
  return sum;
}

double eats2_sprt(double g, double h, double gamma, double delta, double d, double m, double C[], double R[], double N[])
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
      double temp = ats_sprt(g, h, gamma, delta, d, m, Vval, Wval, C, R, N);
      sum += (weight15[j] * weight15[i] * temp * 2.0 * Vval * dgamma(Vval*Vval,c,c) * dnorm(Wval,0.0,sW)) * temp;
    }
  }
  sum *= (0.25*Vdiff*Wdiff);
  return sum;
}

double sdarl_sprt(double g, double h, double gamma, double delta, double m, double C[], double R[], double N[])
{
  return sqrt(earl2_sprt(g, h, gamma, delta, m, C, R, N) - pow(aarl_sprt(g, h, gamma, delta, m, C, R, N),2));
}

double sdats_sprt(double g, double h, double gamma, double delta, double d, double m, double C[], double R[], double N[])
{
  return sqrt(eats2_sprt(g, h, gamma, delta, d, m, C, R, N) - pow(aats_sprt(g, h, gamma, delta, d, m, C, R, N),2));
}

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

double eanos2_sprt(double g, double h, double gamma, double delta, double m, double C[], double R[], double N[])
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
      double temp = anos_sprt(g, h, gamma, delta, m, Vval, Wval, C, R, N);
      sum += (weight15[j] * weight15[i] * temp * 2.0 * Vval * dgamma(Vval*Vval,c,c) * dnorm(Wval,0.0,sW)) * temp;
    }
  }
  sum *= (0.25*Vdiff*Wdiff);
  return sum;
}

double sdanos_sprt(double g, double h, double gamma, double delta, double m, double C[], double R[], double N[])
{
  return sqrt(eanos2_sprt(g, h, gamma, delta, m, C, R, N) - pow(aanos_sprt(g, h, gamma, delta, m, C, R, N),2));
}


double eurl_sprt(double g, double h, double gamma, double delta, double m, double C[], double R[], double N[])
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
      sum += weight15[j] * weight15[i] * 1.0/(1.0 - oc_sprt(g, h, gamma, delta, m, Vval, Wval, C, R, N)) * 2.0 * Vval * dgamma(Vval*Vval,c,c) * dnorm(Wval,0.0,sW);
    }
  }
  
  sum *= (0.25*Vdiff*Wdiff);
  return sum;
}

double eurl22_sprt(double g, double h, double gamma, double delta, double m, double C[], double R[], double N[])
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
      double p = oc_sprt(g, h, gamma, delta, m, Vval, Wval, C, R, N);
      sum += weight15[j] * weight15[i] * (1.0 + p)/((1.0 - p)*(1.0 - p)) * 2.0 * Vval * dgamma(Vval*Vval,c,c) * dnorm(Wval,0.0,sW);
    }
  }
  
  sum *= (0.25*Vdiff*Wdiff);
  return sum;
}

double eurl2_sprt(double g, double h, double gamma, double delta, double m, double C[], double R[], double N[])
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
  
  if(delta == 0.0){
    for(int j = 0; j < 15; ++j){
      double Vval = 0.5 * ((Vdiff)*abscissa15[j] + Vsum);
      for(int i = 0; i < 15; ++i){
        double Wval = 0.5 * (Wdiff)*abscissa15[i] ;
        double p = oc_sprt(g, h, gamma, delta, m, Vval, Wval, C, R, N);
        sum += weight15[j] * weight15[i] * (1.0 + p)/((1.0 - p)*(1.0 - p)) * 2.0 * Vval * dgamma(Vval*Vval,c,c) * dnorm(Wval,0.0,sW);
      }
    }
  }else{
    for(int j = 0; j < 15; ++j){
      double Vval = 0.5 * ((Vdiff)*abscissa15[j] + Vsum);
      for(int i = 0; i < 15; ++i){
        double Wval = 0.5 * (Wdiff)*abscissa15[i] ;
        double p = oc_sprt(g, h, gamma, delta, m, Vval, Wval, C, R, N);
        sum += weight15[j] * weight15[i] * ((2.0*p)/((1.0 - p)*(1.0 - p)) + 1.0/3.0) * 2.0 * Vval * dgamma(Vval*Vval,c,c) * dnorm(Wval,0.0,sW);
      }
    }
  }
  
  sum *= (0.25*Vdiff*Wdiff);
  return sum;
}

double asdrl_sprt(double g, double h, double gamma, double delta, double m, double C[], double R[], double N[])
{
  return sqrt(eurl22_sprt(g, h, gamma, delta, m, C, R, N) - pow(eurl_sprt(g, h, gamma, delta, m, C, R, N),2));
}

double asdts_sprt(double g, double h, double gamma, double delta, double d, double m, double C[], double R[], double N[])
{
  return d * sqrt(eurl2_sprt(g, h, gamma, delta, m, C, R, N) - pow(eurl_sprt(g, h, gamma, delta, m, C, R, N),2));
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

double aaeql_sprt(double g, double h, double gamma, double d, double m, double deltamin, double deltamax, double C[], double R[], double N[])
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
      sum += weight15[j] * weight15[i] * aeql_sprt(g, h, gamma, d, m, deltamin, deltamax, Vval, Wval, C, R, N) * 2.0 * Vval * dgamma(Vval*Vval,c,c) * dnorm(Wval,0.0,sW);
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
double AARL(double g, double h, double gamma, double delta, double m) {
  double C[M], R[M*M], N[M];
  return aarl_sprt(g, h, gamma, delta, m, C, R, N);
}

// [[Rcpp::export]]
double AATS(double g, double h, double gamma, double delta, double d, double m) {
  double C[M], R[M*M], N[M];
  return aats_sprt(g, h, gamma, delta, d, m, C, R, N);
}

// [[Rcpp::export]]
double AANOS(double g, double h, double gamma, double delta, double m) {
  double C[M], R[M*M], N[M];
  return aanos_sprt(g, h, gamma, delta, m, C, R, N);
}

// [[Rcpp::export]]
double SDARL(double g, double h, double gamma, double delta, double m) {
  double C[M], R[M*M], N[M];
  return sdarl_sprt(g, h, gamma, delta, m, C, R, N);
}

// [[Rcpp::export]]
double SDATS(double g, double h, double gamma, double delta, double d, double m) {
  double C[M], R[M*M], N[M];
  return sdats_sprt(g, h, gamma, delta, d, m, C, R, N);
}

// [[Rcpp::export]]
double SDANOS(double g, double h, double gamma, double delta, double m) {
  double C[M], R[M*M], N[M];
  return sdanos_sprt(g, h, gamma, delta, m, C, R, N);
}

// [[Rcpp::export]]
double ASDRL(double g, double h, double gamma, double delta, double m) {
  double C[M], R[M*M], N[M];
  return asdrl_sprt(g, h, gamma, delta, m, C, R, N);
}

// [[Rcpp::export]]
double ASDTS(double g, double h, double gamma, double delta, double d, double m) {
  double C[M], R[M*M], N[M];
  return asdts_sprt(g, h, gamma, delta, d, m, C, R, N);
}

// [[Rcpp::export]]
double AASN(double g, double h, double gamma, double delta, double m) {
  double C[M], R[M*M], N[M];
  return aasn_sprt(g, h, gamma, delta, m, C, R, N);
}

// [[Rcpp::export]]
double AAEQL(double g, double h, double gamma, double d, double m, double deltamin, double deltamax) {
  double C[M], R[M*M], N[M];
  return aaeql_sprt(g, h, gamma, d, m, deltamin, deltamax, C, R, N);
}

// [[Rcpp::export]]
double EAANOS(double g, double h, double gamma, double m, double deltamin, double deltamax) {
  double C[M], R[M*M], N[M];
  return eaanos_sprt(g, h, gamma, m, deltamin, deltamax, C, R, N);
}


// [[Rcpp::export]]
double CARL(double g, double h, double gamma, double delta, double m, double V, double W) {
  double C[M], R[M*M], N[M];
  return arl_sprt(g, h, gamma, delta, m, V, W, C, R, N);
}

// [[Rcpp::export]]
double CATS(double g, double h, double gamma, double delta, double d, double m, double V, double W) {
  double C[M], R[M*M], N[M];
  return ats_sprt(g, h, gamma, delta, d, m, V, W, C, R, N);
}

// [[Rcpp::export]]
double CSDRL(double g, double h, double gamma, double delta, double m, double V, double W) {
  double C[M], R[M*M], N[M];
  return sdrl_sprt(g, h, gamma, delta, m, V, W, C, R, N);
}

// [[Rcpp::export]]
double CSDTS(double g, double h, double gamma, double delta, double d, double m, double V, double W) {
  double C[M], R[M*M], N[M];
  return sdts_sprt(g, h, gamma, delta, d, m, V, W, C, R, N);
}

// [[Rcpp::export]]
double CASN(double g, double h, double gamma, double delta, double m, double V, double W) {
  double C[M], R[M*M], N[M];
  return asn_sprt(g, h, gamma, delta, m, V, W, C, R, N);
}

// [[Rcpp::export]]
double CANOS(double g, double h, double gamma, double delta, double m, double V, double W) {
  double C[M], R[M*M], N[M];
  return anos_sprt(g, h, gamma, delta, m, V, W, C, R, N);
}


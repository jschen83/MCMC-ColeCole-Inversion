/**********************************************************
Author:

  Jinsong Chen
  Lawrence Berkeley National Lab
  Earth Sciences Division
  Berkeley, CA 94720
  jchen@lbl.gov

Notices: 

  (1) The software was developed at Lawrence Berkeley 
      Laboratory by Jinsong Chen of Earth Sciences Division.
      The work is supported by the U. S. Department of 
      Energy. Any commercial use of the software requires
      PRIOR AGREEMENT with Lawrence Berkeley Laboratory. For
      further information, please contact Jinsong Chen at
      jchen@lbl.gov.

  (2) Any use of the software in source and binary forms
     (with or without modification) is permitted with
      agreement of citing the following paper:

  "Chen, J., A. Kemna, and S. Hubbard (2008), A comparison
  between Gauss-Newton and Markov Chain Monte Carlo based
  methods for inverting spectral induced polarization data,
  for Cole-Cole parameters, Geophysics, Vol. 73, No. 6."
***********************************************************/
#include "basefuns.h"
#include "randfuns.h"

/********************************************************
  Long period (>2e+18) uniform random generator based 
  on the algorithm of Pierre L'ecuyer (1988) 
*********************************************************/
double runif(long *seed)
{
  const long m1=2147483563, m2=2147483399,
             a1=40014, a2=40692,
             q1=53668, q2=52774,
             r1=12211, r2=3791;
  static long seed2=123456789;
  long znum,k;
  double rv;
 
  k=(*seed)/q1;
  (*seed)=a1*((*seed)-k*q1)-k*r1;
  if((*seed)<0) (*seed) += m1;

  k=seed2/q2;
  seed2=a2*(seed2-k*q2)-k*r2;
  if(seed2<0) seed2 += m2;

  znum=(*seed)-seed2;
  if(znum<1) znum += (m1-1);

  rv=znum*(1.0/m1);
  return(rv);
}

/********************************************************
  Transformation method for the exponential distribution 
*********************************************************/
double rexp(long *seed)
{
  double tmp,rv;

  while(1) {
    tmp=runif(seed);
    if(tmp<1.0) break;
  };

  rv=(-1.0)*log(1.0-tmp);
  return(rv);
}

/*********************************************************
  Random generator of the standard normal distribution
  using Box-Muller transformation.
**********************************************************/
double rnorm(long *seed)
{
  static int istart=0;
  static double keep_value;
  double u1,u2,r2,coef,nv1,nv2,tmp;

  if(istart==0) {
    do {
       u1=2.0*runif(seed)-1.0;
       u2=2.0*runif(seed)-1.0;
       r2=u1*u1+u2*u2;
    } while (r2>=1.0 || r2==0.0);

    coef=sqrt(-2.0*log(r2)/r2);
    nv1=u1*coef;
    nv2=u2*coef;
  };

  if(istart==0) {
    istart=1;
    tmp=nv1;
    keep_value=nv2;
  } else {
    istart=0;
    tmp=keep_value;
  }

  return(tmp);
}

/*********************************************************
  Generate a random number with gamma distribution
  with shape parameter of a and scale parameter of b.
  
  Reference:
    Marsaglia, G. and W. W. Tsang (2001), A simple method
    for generating gamma variables, ACM Transaction on
    Mathematical Software, Vol. 26, No. 3, P363-372.
**********************************************************/
double rgamma(double a,double b,long *seed)
// a is the shape parameter
// b is the scale parameter
{
  double a1,tmp,d,c,x,v,u,rv;

  //Generate random variable with shape parameter of a
  //and scale parameter of 1.
  if(a<1.0) a1=1.0+a;
  else a1=a;

  d=a1-1.0/3.0;
  c=1.0/sqrt(9.0*d);

  for(;;) {
    do {
       x=rnorm(seed);
       v=1.0+c*x;
    } while (v<=0.0);

    v=v*v*v;
    u=runif(seed);
    if(u<(1.0-0.0331*x*x*x*x)) break;
    if(log(u)<0.5*x*x+d*(1.0-v+log(v))) break;
  };

  rv=d*v;

  if(a<1.0) {
    u=runif(seed);
    tmp=log(rv)+log(u)/a;
    rv=exp(tmp);
  };

  //Get random variables of Gamma(a,b) by transformation
  rv=rv*b;

  return(rv);
}

double dqfs_normal_walk_finite(double x,double *y,double sd,
                               double a,double b,long *seed)
{
  double tmp,log_ratio;

  log_ratio=0.0;
  while(1) {
    tmp=x+sd*rnorm(seed);
    if(tmp>=a && tmp<=b) break;
  };
  (*y)=tmp;
  return(log_ratio);
}

double dqfs_uniform_walk_finite(double x,double *s,
                 double sd,double a,double b,long *seed)
{
  double u,lbd,ubd,a1,b1;
  double log_ratio=0.0;

  lbd=x-sd;
  ubd=x+sd;

  if(lbd<a) {
    u=(double)runif(seed);
    if(u<(ubd-a)/(2*sd)) {
      a1=a;
      b1=ubd;
      (*s)=a1+(double)runif(seed)*(b1-a1);
    }
    else {
      a1=b-(a-lbd);
      b1=b;
      (*s)=a1+(double)runif(seed)*(b1-a1);
    };
  }
 else if(ubd>b) {
    u=runif(seed);
    if(u<(b-lbd)/(2*sd)) {
      a1=lbd;
      b1=b;
      (*s)=a1+(double)runif(seed)*(b1-a1);
    }
    else {
      a1=a;
      b1=a+(ubd-b);
      (*s)=a1+(double)runif(seed)*(b1-a1);
    };
  }
  else (*s)=lbd+(double)runif(seed)*(ubd-lbd);

  return(log_ratio);
}

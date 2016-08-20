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
#include "sisip.h"

extern int Debug,initRand;

ColePost::ColePost(){}
ColePost::~ColePost() 
{
 free(lrho),free(urho);
 free(lm),free(um);
 free(ltau),free(utau);
 free(lc),free(uc);
}

void ColePost::InputData(FILE *fp,long *num)
{
  long i;
  nmdl=num[0];
  nfrq=num[1];
  nvar=3*nmdl+3;
  cole.ReadData(fp,num);
}

void ColePost::InitModel(double *t1,double *a1,double *b1,
	    	          double *t2,double *a2,double *b2,
		          double *t3,double *a3,double *b3,
		          double *t4,double *a4,double *b4,
	 	          double *x,long *idum)
{
  long i,k;
  double tmp,logpdf;

  lrho=dmalloc(1);
  urho=dmalloc(1);
  lm=dmalloc(nmdl);
  um=dmalloc(nmdl);
  ltau=dmalloc(nmdl);
  utau=dmalloc(nmdl);
  lc=dmalloc(nmdl);
  uc=dmalloc(nmdl);

  lrho[0]=(double)a1[0]; 
  urho[0]=(double)b1[0];
  for(i=0;i<nmdl;i++) {
     lm[i]=log((double)a2[i]);
     um[i]=log((double)b2[i]);
     ltau[i]=log(pow(10.0,a3[i]));
     utau[i]=log(pow(10.0,b3[i]));
     lc[i]=(double)a4[i];
     uc[i]=(double)b4[i];
  };

  if(initRand) {
    x[0]=lrho[0]+(urho[0]-lrho[0])*runif(idum);

    for(i=0;i<nmdl;i++) {
      k=1;
      x[i+k]=lm[i]+(um[i]-lm[i])*runif(idum);
    };
    for(i=0;i<nmdl;i++) {
      k=nmdl+1;
      x[i+k]=ltau[i]+(utau[i]-ltau[i])*runif(idum);
    };
    for(i=0;i<nmdl;i++) {
      k=2*nmdl+1;
      x[i+k]=lc[i]+(uc[i]-lc[i])*runif(idum);
    };

    x[3*nmdl+1]=100.0;
    x[3*nmdl+2]=100.0;
  }
  else {
    x[0]=t1[0];
    for(i=0;i<nmdl;i++) {
      x[i+1]=log(t2[i]);
      x[i+nmdl+1]=log(pow(10.0,t3[i]));
      x[i+2*nmdl+1]=t4[i];
    };

    x[3*nmdl+1]=100.0;
    x[3*nmdl+2]=100.0;
  };

  logpdf=LogFunc(x,&setup);
  printf("setup=%f,logpdf=%f\n",setup,logpdf);
}

void ColePost::UpdateModel(void)
{
  mse=setup;
}

double ColePost::LogFunc(double *x,double *msft)
{
  long i,num;
  double coleLKHD,sumf;
  double *rho,*m,*tau,*c,myerr[2];
  
  rho=dmalloc(1);
  m=dmalloc(nmdl);
  tau=dmalloc(nmdl);
  c=dmalloc(nmdl);
  GetParaBack(x,rho,m,tau,c,myerr);

  coleLKHD=cole.LogLkhd(rho,m,tau,c,myerr,msft);
  setup=(*msft);
  sumf=coleLKHD;

  free(rho);
  free(m);
  free(tau);
  free(c);
  return(sumf);
}

void ColePost::Output(double *x,long niter,long M,
                      float *outall,double *misfit)
{
  long i,k;
  double *rho,*m,*tau,*c,myerr[2];
  float *mytau;

  rho=dmalloc(1);
  m=dmalloc(nmdl);
  tau=dmalloc(nmdl);
  mytau=fmalloc(nmdl);
  c=dmalloc(nmdl);
  GetParaBack(x,rho,m,tau,c,myerr);
  for(i=0;i<nmdl;i++) mytau[i]=(float)log10(tau[i]);

  outall[niter]=(float)rho[0];
  for(i=0;i<nmdl;i++) 
     outall[(i+1)*M+niter]=(float)m[i];
  for(i=0;i<nmdl;i++) 
     outall[(i+nmdl+1)*M+niter]=(float)mytau[i];
  for(i=0;i<nmdl;i++) 
     outall[(i+2*nmdl+1)*M+niter]=(float)c[i];

  for(i=0;i<2;i++) 
     outall[(i+3*nmdl+1)*M+niter]=(float)myerr[i];

  misfit[niter]=mse;

  if(Debug==1) { 
    printf("rho:");
    printf("%8.4f ",rho[0]);
    printf("\n");

    printf("m:");
    for(i=0;i<nmdl;i++) printf("%10.6le ",m[i]);
    printf("\n");

    printf("tau:");
    for(i=0;i<nmdl;i++) printf("%10.6le ",mytau[i]);
    printf("\n");

    printf("c:");
    for(i=0;i<nmdl;i++) printf("%10.6le ",c[i]);
    printf("\n");

    printf("Relative_error in real part=%f\n",sqrt(1.0/myerr[0]));
    printf("Relative_error in imag part=%f\n",sqrt(1.0/myerr[1]));

    printf("Misfit=%f\n",mse);
  };

  free(mytau);
  free(rho),free(m),free(tau),free(c);
}

void ColePost::GetParaBack(double *x,double *rho,
                 double *m,double *tau,double *c,
                 double *myerr)
{
  long i;

  rho[0]=x[0];

  for(i=0;i<nmdl;i++) {
//     m[i]=x[i+1];
     m[i]=exp(x[i+1]); //Using log-scale for chargeability
     tau[i]=exp(x[i+nmdl+1]);
     c[i]=x[i+2*nmdl+1];
  };

 myerr[0]=x[3*nmdl+1];
 myerr[1]=x[3*nmdl+2];
}

void ColePost::GetBounds(double *lower,double *upper,
                          double *d,long *num)
{
  long i,k;

  lower[0]=lrho[0];
  upper[0]=urho[0];
  d[0]=(upper[0]-lower[0])/num[0];

  for(i=0;i<nmdl;i++) {
    k=1;
    lower[i+k]=lm[i];
    upper[i+k]=um[i];
    d[i+k]=(upper[i+k]-lower[i+k])/num[1];

    k=nmdl+1;
    lower[k+i]=ltau[i];
    upper[k+i]=utau[i];
    d[k+i]=(upper[k+i]-lower[k+i])/num[2];

    k=2*nmdl+1;
    lower[k+i]=lc[i];
    upper[k+i]=uc[i];
    d[k+i]=(upper[k+i]-lower[k+i])/num[3];
  };
}

void ColePost::GetCoef(double *rtn,double *x)
{
  double *rho,*m,*tau,*c,myerr[2];

  rho=dmalloc(1);
  m=dmalloc(nmdl);
  tau=dmalloc(nmdl);
  c=dmalloc(nmdl);

  GetParaBack(x,rho,m,tau,c,myerr);
  cole.GetCoef(rtn,myerr,rho,m,tau,c);

  free(rho);
  free(m);
  free(tau);
  free(c);
}

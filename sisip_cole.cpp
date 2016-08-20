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

extern int emxyYES;

COLE::COLE() {}
COLE::~COLE() 
{ 
  if(readID==1) {
    free(freq);
    free(rldt);
    free(igdt);}
  };

void COLE::ReadData(FILE *fp,long *num)
{
  long i,myDebug=0;
  double amp,phase,tmp;

  readID=1;
  nmdl=num[0];
  nfrq=num[1];
  freq=dmalloc(nfrq);
  rldt=dmalloc(nfrq);
  igdt=dmalloc(nfrq);

  /*Real and imaginary components of data*/

  for(i=0;i<nfrq;i++) {
    if(emxyYES)
      fscanf(fp,"%lf %lf %lf",freq+i,rldt+i,igdt+i);
    else {
      fscanf(fp,"%lf %lf %lf",freq+i,&amp,&phase);
      rldt[i]=amp*cos(-0.001*phase);
      igdt[i]=amp*sin(-0.001*phase);
    };

    if(myDebug)  
      printf("%f %f %f\n",freq[i],rldt[i],igdt[i]);
  };
}

void COLE::Forward(FILE *fp,long *num,float *freq,
           double *rho,double *m,double *tau,double *c,
           float ns_level,long *idum)
{
  int AddNoise=1,myDebug=1;
  const double pi=3.1415926;
  double w,R,I,a,b,Rx,Ry,tmp1,tmp2;
  double sumrl,sumig;
  long i,j;

  readID=0;
  nmdl=num[0];
  nfrq=num[1];

  if(myDebug) {
    printf("rho=%f\n",rho[0]);
    for(i=0;i<nmdl;i++) printf("m=%f\n",m[i]);
    for(i=0;i<nmdl;i++) printf("tau=%le\n",tau[i]);
    for(i=0;i<nmdl;i++) printf("c=%f\n",c[i]);
  };

  for(i=0;i<nfrq;i++) {
    w=2.0*pi*(double)freq[i];

    sumrl=sumig=0.0;
    for(j=0;j<nmdl;j++) {
      R=pow(w*tau[j],c[j])*cos(0.5*pi*c[j])+1.0;
      I=pow(w*tau[j],c[j])*sin(0.5*pi*c[j]);
      a=R/(R*R+I*I);
      b=I/(R*R+I*I);
      sumrl += m[j]*(1.0-a);
      sumig += m[j]*b;
    };

    Rx=rho[0]*(1.0-sumrl);
    Ry=(-1.0)*rho[0]*sumig;

    if(AddNoise) {
      Rx=(1.0+ns_level*rnorm(idum))*Rx;
      Ry=(1.0+ns_level*rnorm(idum))*Ry;
    };

    fprintf(fp,"%f  %f  %f\n",freq[i],Rx,Ry);
  };
}

double COLE::LogLkhd(double *rho,double *m,
                     double *tau,double *c,
                     double *myerr,double *msft)
{
  int myDebug=0;
  const double pi=3.1415926;
  double w,R,I,a,b,Rx,Ry,sumrl,sumig,RSS;
  double tmp1,tmp2,sum11,sum12,sum21,sum22,
         lkhd1,lkhd2,BIC1,BIC2,loglkhd,AIC1,AIC2;
  long i,j,num,nfd;

  if(myDebug) {
    printf("rho=%f\n",rho[0]);
    for(i=0;i<nmdl;i++) printf("m=%f\n",m[i]);
    for(i=0;i<nmdl;i++) printf("tau=%le\n",tau[i]);
    for(i=0;i<nmdl;i++) printf("c=%f\n",c[i]);
  };

  sum11=sum12=sum21=sum22=0.0;
  for(i=0;i<nfrq;i++) {
    w=2.0*pi*(double)freq[i];
    sumrl=sumig=0.0;

    for(j=0;j<nmdl;j++) {
      R=pow(w*tau[j],c[j])*cos(0.5*pi*c[j])+1.0;
      I=pow(w*tau[j],c[j])*sin(0.5*pi*c[j]);

      a=R/(R*R+I*I);
      b=I/(R*R+I*I);
      sumrl += m[j]*(1.0-a);
      sumig += m[j]*b;
    };

    Rx=rho[0]*(1.0-sumrl);
    Ry=(-1.0)*rho[0]*sumig;

    tmp1=(Rx-rldt[i])/rldt[i];
    sum11 += tmp1*tmp1;
    sum21 += tmp1*tmp1*myerr[0];

    tmp2=(Ry-igdt[i])/igdt[i];
    sum12 += tmp2*tmp2;
    sum22 += tmp2*tmp2*myerr[1];
  };

  (*msft)= sqrt((sum11+sum12)/(2*nfrq-1));
  loglkhd=0.5*nfrq*(log(myerr[0])+log(myerr[1]))
          -0.5*(sum21+sum22);

  if(myDebug) {
    printf("misfit=%f, loglkhd=%f\n",(*msft),loglkhd);
    exit(1);
  };

  return(loglkhd);
}

void COLE::GetCoef(double *rtn,double *myerr,
                   double *rho,double *m,
                   double *tau,double *c)
{
  int myDebug=0;
  const double pi=3.1415926;
  double w,R,I,a,b,sumrl,sumig,
         Aj,Bj,Cj,Dj,tmpmu,tmptau,
         sum1,sum2,sumr,sumr2,sums,sums2;
  long i,j;

  sumr=sums=0.0;
  sumr2=sums2=0.0;
  sum1=sum2=0.0;

  for(i=0;i<nfrq;i++) {

    w=2.0*pi*freq[i];
    sumrl=sumig=0.0;

    for(j=0;j<nmdl;j++) {
      R=pow(w*tau[j],c[j])*cos(0.5*pi*c[j])+1.0;
      I=pow(w*tau[j],c[j])*sin(0.5*pi*c[j]);
      a=R/(R*R+I*I);
      b=I/(R*R+I*I);

      sumrl += m[j]*(1.0-a);
      sumig += m[j]*b;
    };

    Aj=rldt[i];
    Bj=igdt[i];
    Cj=(1.0-sumrl);
    Dj=(-1.0)*sumig;

    sumr += Cj/Aj;
    sumr2 += (Cj/Aj)*(Cj/Aj);
    sums += Dj/Bj;
    sums2 += (Dj/Bj)*(Dj/Bj);

    sum1 += (1.0-Cj*rho[0]/Aj)*(1.0-Cj*rho[0]/Aj);
    sum2 += (1.0-Dj*rho[0]/Bj)*(1.0-Dj*rho[0]/Bj);
  };

  if(myDebug) {
    printf("sumr=%f,sumr2=%f\n",sumr,sumr2);
    printf("sums=%f,sums2=%f\n",sums,sums2);
  };

  tmptau=sumr2*myerr[0]+sums2*myerr[1];
  tmpmu=(sumr*myerr[0]+sums*myerr[1])/tmptau;
   
  rtn[0]=tmpmu;
  rtn[1]=sqrt(1.0/tmptau);
  rtn[2]=sum1;
  rtn[3]=sum2;
  rtn[4]=(double)nfrq;
}

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

extern int Debug;

SAMP::SAMP() {}
SAMP::~SAMP() {free(x);}

void SAMP::InputData(FILE *fp,long *num,
		double *irho,double *lrho,double *urho,
		double *im,double *lm,double *um,
	        double *itau,double *ltau,double *utau,
	        double *ic,double *lc,double *uc,
                long *idum)
{
  post.InputData(fp,num);
  xsize=post.nvar;
  nmdl=post.nmdl;
  x=dmalloc(xsize);
  post.InitModel(irho,lrho,urho,im,lm,um,itau,ltau,utau,
                 ic,lc,uc,x,idum);
}

void SAMP::TunePara(long *tmpnum)
{
  long i;
  for(i=0;i<4;i++) num[i]=tmpnum[i];
}

void SAMP::UpdateRho(long *idum)
{
  double oldsum,newsum,msft;
  double r1,r2,sum_ratio,my,*x1,tmpval,
         mymu,mystd,rtn[5];
  double *w,*lower,*upper;

  post.GetCoef(rtn,x);
  mymu=rtn[0];
  mystd=rtn[1];

  w=dmalloc(xsize);
  lower=dmalloc(xsize);
  upper=dmalloc(xsize);
  post.GetBounds(lower,upper,w,num);

  //Check consistency of the initial rho data
  if(mymu<lower[0] | mymu>upper[0]) {
    printf("mymu=%f, mystd=%f\n",mymu,mystd);
    printf("The initial values are not consistent!\n");
    exit(1);
  };

  r2=dqfs_normal_walk_finite(mymu,&tmpval,mystd,
                 lower[0],upper[0],idum);
  x[0]=tmpval;
  post.UpdateModel();

  free(lower);
  free(upper);
  free(w);
}

void SAMP::UpdateTauAB(long *idum)
{
  int i,k,n;
  double rtn[5];
  double alpha0=1e-3,lambda0=1e-3; //priors
  double alpha,beta,lambda;

  post.GetCoef(rtn,x);
  n=(int)rtn[4];

  k=3*nmdl+1;
  alpha=alpha0+0.5*n;
  for(i=0;i<2;i++) {
    lambda=lambda0+0.5*rtn[i+2];
    beta=1.0/lambda;
    x[k+i]=rgamma(alpha,beta,idum);
  };

  post.UpdateModel();
}

void SAMP::MultipleMH(int *myset,int ns,long *idum)
{
  int i,j,mixID;
  double oldsum,newsum,misfit;
  double r1,r2,sum_ratio,my,*x1,sd,tmp;
  double *w,*lower,*upper,u,p[]={0.5,0.5};

  w=dmalloc(xsize);
  lower=dmalloc(xsize);
  upper=dmalloc(xsize);
  post.GetBounds(lower,upper,w,num);

  x1=dmalloc(xsize);
  for(i=0;i<xsize;i++) {
    x1[i]=x[i];
    w[i] *=0.1;
  };

  r2=0.0;
  for(j=0;j<ns;j++) {
    i=myset[j]-1;
    if(mixID==1)
      r2 +=dqfs_normal_walk_finite(x[i],&tmp,
                w[i],lower[i],upper[i],idum);
    else
      r2 +=dqfs_uniform_walk_finite(x[i],&tmp,
                w[i],lower[i],upper[i],idum);

    x1[i]=tmp;
  };

  oldsum=post.LogFunc(x,&misfit);
  newsum=post.LogFunc(x1,&misfit);

  r1=newsum-oldsum;
  sum_ratio=r1+r2;
  my=(sum_ratio>=0 ? 1.0:exp(sum_ratio));
  u=runif(idum);
  if(u<my) {
     for(j=0;j<ns;j++) {
       i=myset[j]-1;
       x[i]=x1[i];
     };
    post.UpdateModel();
  };

  free(w);
  free(lower);
  free(upper);
  free(x1);
}

void SAMP::MultipleSS0(int *myset,int ns,long *idum)
{
  int i,j,myDebug=0;
  double *x1,z,gx0,gx1,misfit;
  double *w,*R,*L,*lower,*upper;

  w=dmalloc(xsize);
  lower=dmalloc(xsize);
  upper=dmalloc(xsize);
  post.GetBounds(lower,upper,w,num);

  L=dmalloc(xsize);
  R=dmalloc(xsize);
  for(i=0;i<xsize;i++) {
    L[i]=lower[i];
    R[i]=upper[i];
    if(myDebug)
      printf("lower[%d]=%f,upper[%d]=%f\n",
              i,lower[i],i,upper[i]);
  };
  gx0=post.LogFunc(x,&misfit);
  z=gx0-rexp(idum);
  if(myDebug) printf("z=%f,gx0=%f\n",z,gx0);

  x1=dmalloc(xsize);
  for(i=0;i<xsize;i++) x1[i]=x[i];

  while(1) {
    for(j=0;j<ns;j++) {
       i=myset[j]-1;
       x1[i]=(double)(L[i]+(R[i]-L[i])*runif(idum));
       if(myDebug) printf("x[%d]=%f,x1[%d]=%f\n",
                   i,x[i],i,x1[i]);
    };

    gx1=post.LogFunc(x1,&misfit);
    if(myDebug) printf("z=%f,gx1=%f\n",z,gx1);

    if(z<gx1) {
      for(j=0;j<ns;j++) {
       i=myset[j]-1;
       x[i]=x1[i];
      };
      post.UpdateModel();
      break;
    }
    else
      for(j=0;j<ns;j++) {
       i=myset[j]-1;
       if(x1[i]<x[i]) L[i]=(float)x1[i];
       else R[i]=(float)x1[i];
      };
  };

  free(x1),free(L),free(R);
  free(w),free(lower),free(upper);
}

void SAMP::MultipleSS1(int *myset,int ns,long *idum)
{
  int i,j;
  double *x1,z,gx0,gx1,misfit;
  double *w,*R,*L,*lower,*upper;

  w=dmalloc(xsize);
  lower=dmalloc(xsize);
  upper=dmalloc(xsize);
  post.GetBounds(lower,upper,w,num);

  L=dmalloc(xsize);
  R=dmalloc(xsize);
  for(i=0;i<xsize;i++) {
     L[i]=x[i]-w[i]*runif(idum);
     if(L[i]<lower[i]) L[i]=lower[i];
     R[i]=L[i]+w[i];
     if(R[i]>upper[i]) R[i]=upper[i];
  };

  gx0=post.LogFunc(x,&misfit);
  z=gx0-rexp(idum);

  x1=dmalloc(xsize);
  for(i=0;i<xsize;i++) x1[i]=x[i];

  while(1) {
    for(j=0;j<ns;j++) {
      i=myset[j]-1;
      x1[i]=(double)(L[i]+(R[i]-L[i])*runif(idum));
    };
    gx1=post.LogFunc(x1,&misfit);

    if(z<=gx1) {
      for(j=0;j<ns;j++) {
        i=myset[j]-1;
        x[i]=x1[i];
      };
      post.UpdateModel();
      break;
    }
    else
      for(j=0;j<ns;j++) {
        i=myset[j]-1;
        if(x1[i]<x[i]) L[i]=(float)x1[i];
        else R[i]=(float)x1[i];
      };
  };

  free(x1),free(L),free(R);
  free(w),free(lower),free(upper);
}

void SAMP::UpdateChain(float *p,long *idum)
{
  int set1[]={1},ns1=1;
  int *set2,*set3,*set4,ns2,ns3,ns4;
  int *allset,allns;
  long i,sampID;
  double u;

  ns2=ns3=ns4=nmdl;
  set2=imalloc(ns2);
  set3=imalloc(ns3);
  set4=imalloc(ns4);
  for(i=0;i<ns2;i++) set2[i]=2+3*i;
  for(i=0;i<ns3;i++) set3[i]=3+3*i;
  for(i=0;i<ns4;i++) set4[i]=4+3*i;

  allns=3*nmdl+1;
  allset=imalloc(allns);
  for(i=0;i<allns;i++) allset[i]=i+1;

  u=runif(idum);
  if(u<=p[0]) sampID=1;
  else if(u<=p[0]+p[1]) sampID=2;
  else sampID=3;

  switch(sampID) {
  case 1:
    u=runif(idum);
    if(u<0.2) 
      MultipleMH(set2,ns2,idum);
    else if(u<0.4)
      MultipleMH(set3,ns3,idum);
    else if(u<0.6)
      MultipleMH(set4,ns4,idum);
    else if(u<0.8)
      MultipleMH(allset,allns,idum);
    else
      UpdateRho(idum);
    break;

  case 2:
    u=runif(idum);
    if(u<0.2) 
      MultipleSS0(set2,ns2,idum);
    else if(u<0.4)
      MultipleSS0(set3,ns3,idum);
    else if(u<0.6)
      MultipleSS0(set4,ns4,idum);
    else if(u<0.8)
      MultipleSS0(allset,allns,idum);
    else 
      UpdateRho(idum);
    break;

  case 3: 
    u=runif(idum);
    if(u<0.2)
      MultipleSS1(set2,ns2,idum);
    else if(u<0.4)
      MultipleSS1(set3,ns3,idum);
    else if(u<0.6)
      MultipleSS1(set4,ns4,idum);
    else if(u<0.8)
      MultipleSS1(allset,allns,idum);
    else
      UpdateRho(idum);
    break;

  default :
    printf("No method for updating\n");
    exit(1);
  };

  UpdateTauAB(idum);

  free(set2),free(set3),free(set4),free(allset);
}

void SAMP::SaveChain(long niter,long M,
                      float *outall,double *misfit)
{
  post.Output(x,niter,M,outall,misfit);
}

void SAMP::OutputChain(FILE *fp1,FILE *fp2,FILE *fp3,long M,
		        float *outall,double *msft)
{
  long i;
  
  fwrite(outall,sizeof(float),xsize*M,fp1);
  for(i=1;i<=xsize;i++)
    fprintf(fp2,"Var%d    %6d    %6d\n",i,(i-1)*M+1,i*M);
  for(i=1;i<=M;i++) 
    fprintf(fp3,"%le\n",msft[i]);
}

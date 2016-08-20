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

void ReadModelParam(FILE *fp,FILE *fpout,long *num,
	  double *irho,double *lrho,double *urho,  
	  double *im,double *lm,double *um, 
	  double *itau,double *ltau,double *utau, 
	  double *ic,double *lc,double *uc) 
{
  long i,nmdl,nfrq;
  char buf[101],title[101],name[101];

  fgets(title,100,fp); //Title of the data file
  if(fgets(buf,100,fp)!=NULL) sscanf(buf,"%ld", num+0);
  if(fgets(buf,100,fp)!=NULL) sscanf(buf,"%ld", num+1);

  nmdl=num[0];
  nfrq=num[1];

  fgets(name,100,fp); 
  for(i=0;i<1;i++) {
    if(fgets(buf,100,fp)!=NULL)
      sscanf(buf,"%lf %lf %lf",irho+i,lrho+i,urho+i);
  };

  fgets(name,100,fp); 
  for(i=0;i<nmdl;i++) {
    if(fgets(buf,100,fp)!=NULL)
      sscanf(buf,"%lf %lf %lf",im+i,lm+i,um+i);
  };

  fgets(name,100,fp); 
  for(i=0;i<nmdl;i++) {
    if(fgets(buf,100,fp)!=NULL)
      sscanf(buf,"%lf %lf %lf",itau+i,ltau+i,utau+i);
  };

  fgets(name,100,fp); 
  for(i=0;i<nmdl;i++) {
    if(fgets(buf,100,fp)!=NULL)
      sscanf(buf,"%lf %lf %lf",ic+i,lc+i,uc+i);
  };

  if(Debug) {
    fprintf(fpout,"%s",title);
    fprintf(fpout,"nmdl=%d\n",num[0]);
    fprintf(fpout,"nfrq=%d\n",num[1]);

    fprintf(fpout,"irho   lrho  urho\n");
    for(i=0;i<1;i++)
      fprintf(fpout,"%le %le %le\n", 
        irho[i],lrho[i],urho[i]);

    fprintf(fpout,"im   lm  um\n");
    for(i=0;i<nmdl;i++)
      fprintf(fpout,"%le %le %le\n", im[i],lm[i],um[i]);

    fprintf(fpout,"itau   ltau  utau\n");
    for(i=0;i<nmdl;i++)
      fprintf(fpout,"%le %le %le\n", itau[i],ltau[i],utau[i]);

    fprintf(fpout,"ic   lc  uc\n");
    for(i=0;i<nmdl;i++)
      fprintf(fpout,"%le %le %le\n", ic[i],lc[i],uc[i]);
  }
}

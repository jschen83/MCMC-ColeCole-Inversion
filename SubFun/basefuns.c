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

int *imalloc(int num)
{
  int *lptr;

  if((lptr=(int*)malloc(num*sizeof(int)))==NULL) {
    fprintf(stderr, "No enough memories for int pointers.\n");
    exit(1);
    };

  return(lptr);
}

float *fmalloc(int num)
{
  float *dptr;

  dptr=(float*)malloc(num*sizeof(float));
  if(dptr==NULL) {
    fprintf(stderr, "No enough memories for float pointers.\n");
    exit(1);
    };
  return(dptr);
}

double *dmalloc(int num)
{
  double *dptr;

  if((dptr=(double*)malloc(num*sizeof(double)))==NULL) {
    fprintf(stderr, "No enough memories for double pointers.\n");
    exit(1);
    };
  return(dptr);
}

FILE *myfopen(char *filename, char *mode)
{
  FILE *tmp_fp;

  tmp_fp=fopen(filename, mode);
  if(tmp_fp==NULL) {
    fprintf(stderr, "The file %s cannot be opened.\n",filename);
    exit(1);
    };

  return(tmp_fp);
}

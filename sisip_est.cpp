/****************************************************************
Author:

  Jinsong Chen
  Lawrence Berkeley National Lab
  Earth Sciences Division
  Berkeley, CA 94720
  jchen@lbl.gov

Developed: January, 2004
Last updated: February, 2010

Usage:

  (1) This is the main function for calling the forward model
      to calculate Cole-Cole responses for a given set of Cole
      Cole parameters. The Cole-Cole parameters can be provided
      in any file (or in any form), for example *.stat here.

  (2) When estID=1, you can modify the code and read your own
      parameter files.

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

WARRANTY DISCLAIMER: 

  THE SOFTWARE IS SUPPLIED "AS IS" WITHOUT WARRANTY OF ANY
  KIND. BERKELEY LAB, ITS LICENSORS, THE UNITED STATES, THE
  UNITED STATES DEPARTMENT OF ENERGY, AND THEIR EMPLOYEES:
  (1) DISCLAIM ANY WARRANTIES, EXPRESS OR IMPLIED, INCLUDING
  BUT NOT LIMITED TO ANY IMPLIED WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE, TITLE OR NON-INFRINGEMENT,
  (2) DO NOT ASSUME ANY LEGAL LIABILITY OR RESPONSIBILITY FOR
  THE ACCURACY, COMPLETENESS, OR USEFULNESS OF THE SOFTWARE,
  (3) DO NOT REPRESENT THAT USE OF THE SOFTWARE WOULD NOT 
  INFRINGE PRIVATELY OWNED RIGHTS, (4) DO NOT WARRANT THAT
  THE SOFTWARE WILL FUNCTION UNINTERRUPTED, THAT IT IS ERROR
  -FREE OR THAT ANY ERRORS WILL BE CORRECTED.

LIMITATION OF LIABILITY:

  IN NO EVENT WILL BERKELEY LAB OR ITS LICENSORS BE LIABLE FOR 
  ANY INDIRECT, INCIDENTAL, CONSEQUENTIAL, SPECIAL OR PUNITIVE
  DAMAGES OF ANY KIND OR NATURE, INCLUDING BUT NOT LIMITED TO
  LOSS OF PROFITS OR LOSS OF DATA, FOR ANY REASON WHATSOEVER, 
  WHETHER SUCH LIABILITY IS ASSERTED ON THE BASIS OF CONTRACT,
  TORT (INCLUDING NEGLIGENCE OR STRICT LIABILITY), OR OTHERWISE,
  EVEN IF BERKELEY LAB HAS BEEN WARNED OF THE POSSIBILITY OF 
  SUCH LOSS OR DAMAGES. IN NO EVENT SHALL BERKELEY'S LIABILITY
  FOR DAMAGES ARISING FROM OR IN CONNECTION WITH THIS AGREEMENT
  EXCEED THE AMOUNT PAID BY YOU FOR THE SOFTWARE.
****************************************************************/
#include "sisip.h"

int Debug=1;
int dataTYPE=1;
int lkhdID=1;
int optionID=1;
int emxyYES=1;

int main(void)
{
  FILE *fp1,*fp2,*fp3;
  char ndata[60],nstat[60],nest[60];
  int estID=2;//1-Read .txt file and 2-Read .stat files 
  int myID=1; //1-init1, 2-init, 3-best
  int typeID=1;//1-median and 2-mode
  long i,nmdl,nfrq,nvar,num[2],idum=-17L;
  COLE cole;
  float *freq,noise_level=0.0;
  double tmp1,tmp2,*rho,*m,*tau,*c;
  char title[201],buf[201];

  nmdl=2;
  rho=dmalloc(1);
  m=dmalloc(nmdl);
  tau=dmalloc(nmdl);
  c=dmalloc(nmdl);

  strcpy(ndata,"");
  strcat(ndata,"./input/example.txt");
  strcpy(nstat,"");
  strcat(nstat,"./output/example_mod2.stat");
  strcpy(nest,"");
  strcat(nest,"./output/example_mod2.est");

  fp1=myfopen(ndata,"r");
  freq=fmalloc(MAXFREQ);
  i=0;
  for(;;) {
    fscanf(fp1,"%f %lf %lf",freq+i,&tmp1,&tmp2);
    if(feof(fp1)==1) break;
    else i++;
  };
  fclose(fp1);
  nfrq=i;
  num[0]=nmdl;
  num[1]=nfrq;

  if(estID==1) {//Replace it with your own file
    fp3=myfopen("./DATak/lab1_ak_summary.txt","r");
    fgets(buf,200,fp3);

    fgets(buf,200,fp3);
    if(myID==1)
      sscanf(buf,"%lf %lf %lf %lf",rho,&tmp1,&tmp1,&tmp1);
    else if(myID==2)
      sscanf(buf,"%lf %lf %lf %lf",&tmp1,rho,&tmp1,&tmp1);
    else if(myID==3)
      sscanf(buf,"%lf %lf %lf %lf",&tmp1,&tmp1,&tmp1,rho);
    printf("rho=%f\n",rho[0]);

    for(i=0;i<nmdl;i++) {
      fgets(buf,200,fp3);
      if(myID==1)
        sscanf(buf,"%lf %lf %lf %lf",m+i,&tmp1,&tmp1,&tmp1);
      else if(myID==2)
        sscanf(buf,"%lf %lf %lf %lf",&tmp1,m+i,&tmp1,&tmp1);
      else if(myID==3)
        sscanf(buf,"%lf %lf %lf %lf",&tmp1,&tmp1,&tmp1,m+i);
      printf("m=%f\n",m[i]);
    };
    for(i=0;i<nmdl;i++) {
      fgets(buf,200,fp3);
      if(myID==1)
      sscanf(buf,"%lf %lf %lf %lf",tau+i,&tmp1,&tmp1,&tmp1);
      else if(myID==2)
        sscanf(buf,"%lf %lf %lf %lf",&tmp1,tau+i,&tmp1,&tmp1);
      else if(myID==3)
        sscanf(buf,"%lf %lf %lf %lf",&tmp1,&tmp1,&tmp1,tau+i);
      tau[i]=pow(10,tau[i]);
      printf("tau=%le\n",tau[i]);
    };
    for(i=0;i<nmdl;i++) {
      fgets(buf,200,fp3);
      if(myID==1)
        sscanf(buf,"%lf %lf %lf %lf",c+i,&tmp1,&tmp1,&tmp1);
      else if(myID==2)
        sscanf(buf,"%lf %lf %lf %lf",&tmp1,c+i,&tmp1,&tmp1);
      else if(myID==3)
        sscanf(buf,"%lf %lf %lf %lf",&tmp1,&tmp1,&tmp1,c+i);
      printf("c=%le\n",c[i]);
    };
    fclose(fp3);

    if(myID==1)
      fp1=myfopen("./DATak/lab1_ak_med_best.est","w");
    else if(myID==2)
      fp1=myfopen("./estMysyn/mysyn10_init2.est","w");
    else if(myID==3)
      fp1=myfopen("./DATak/lab3_ak_best.est","w");
   }

  else {
    fp1=myfopen(nstat,"r");
    fgets(title,200,fp1);
    nvar=3*nmdl+1;
    fgets(buf,200,fp1);
    if(typeID==1) //median
      sscanf(buf,"%lf %lf %lf %lf %lf %lf %lf %lf",&tmp1,
             rho,&tmp1,&tmp1,&tmp1,&tmp1,&tmp1,&tmp1);
    else if(typeID==2) //mode
      sscanf(buf,"%lf %lf %lf %lf %lf %lf %lf %lf",
             &tmp1,&tmp1,&tmp1,&tmp1,&tmp1,rho,&tmp1,&tmp1);
    printf("rho=%f\n",rho[0]);

    for(i=0;i<nmdl;i++) {
       fgets(buf,200,fp1);
       if(typeID==1) //median
       sscanf(buf,"%lf %lf %lf %lf %lf %lf %lf %lf",&tmp1,m+i,
                  &tmp1,&tmp1,&tmp1,&tmp1,&tmp1,&tmp1);
       else if(typeID==2) //mode
       sscanf(buf,"%lf %lf %lf %lf %lf %lf %lf %lf",
                  &tmp1,&tmp1,&tmp1,&tmp1,&tmp1,m+i,
                  &tmp1,&tmp1);
       printf("m=%f\n",m[i]);
    };

    for(i=0;i<nmdl;i++) {
      fgets(buf,200,fp1);
      if(typeID==1)
      sscanf(buf,"%lf %lf %lf %lf %lf %lf %lf %lf",&tmp1,tau+i,
                &tmp1,&tmp1,&tmp1,&tmp1,&tmp1,&tmp1);
      else if(typeID==2)
      sscanf(buf,"%lf %lf %lf %lf %lf %lf %lf %lf",&tmp1,&tmp1,
                &tmp1,&tmp1,&tmp1,tau+i,&tmp1,&tmp1);
      tau[i]=pow(10.0,tau[i]);
      printf("tau=%le\n",tau[i]);
    };

    for(i=0;i<nmdl;i++) {
      fgets(buf,200,fp1);
      if(typeID==1)
      sscanf(buf,"%lf %lf %lf %lf %lf %lf %lf %lf",&tmp1,c+i,
                &tmp1,&tmp1,&tmp1,&tmp1,&tmp1,&tmp1);
      else if(typeID==2)
      sscanf(buf,"%lf %lf %lf %lf %lf %lf %lf %lf",&tmp1,&tmp1,
                &tmp1,&tmp1,&tmp1,c+i,&tmp1,&tmp1);
      printf("c=%f\n",c[i]);
    };
    fclose(fp1);

    fp1=myfopen(nest,"w");
  };

  cole.Forward(fp1,num,freq,rho,m,tau,c,noise_level,&idum);
  fclose(fp1);
  free(rho);
  free(m);
  free(tau);
  free(c);
  free(freq);
  return(0);
}

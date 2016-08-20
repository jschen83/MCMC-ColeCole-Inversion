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

  (1) This is the main function for setting up and for 
      drawing samples. The default parameters are good for
      most cases. You can run single or multiple chains.
      The output includes three files in the format of R 
      package "coda". You can use it to analyze the chains.
      *.bin is a binary file, and *.ind and *.msft are text
      files. Please only use the later half of the samples for
      analysis.

  (2) To achieve fast convergence, you can adjust the 
      following parameters:  
      (a) Fraction of using different sampling methods 
          (MMH, MSS0, and MSS1). Typically, the portion 
          of MSS0 should be small (say less than 10%).
      (b) Random seed. Change the random seed and rerun.
      (c) The sampling size "num[]". The default values
          are good for many cases.

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
#include "./SubFun/timer.h"

int Debug=1;
int initRand=0;
int emxyYES=0;

int main(void)
{
  FILE *fp1,*fp2,*fp3;
  double *irho,*lrho,*urho,*im,*lm,*um,
         *itau,*ltau,*utau,*ic,*lc,*uc;
  double *misfit;
  float *outall;
  long numID=500,saveID=10,numM=400000;
  long i,j,k,nvar,niter,npar[2];
  long M,Debug2;
  int sampID;
  SAMP mcmc;

  /*Parameters can be adjusted to achieve fast convergence*/
  float p[]={0.4,0.0,0.6};//1st-MMH, 2nd-MSS0, 3rd-MSS1
  long idum=-11L; //random seed
  long num[]={10,10,20,10}; //interval size


  irho=dmalloc(MAXLAYER);
  lrho=dmalloc(MAXLAYER);
  urho=dmalloc(MAXLAYER);
  im=dmalloc(MAXLAYER);
  lm=dmalloc(MAXLAYER);
  um=dmalloc(MAXLAYER);
  itau=dmalloc(MAXLAYER);
  ltau=dmalloc(MAXLAYER);
  utau=dmalloc(MAXLAYER);
  ic=dmalloc(MAXLAYER);
  lc=dmalloc(MAXLAYER);
  uc=dmalloc(MAXLAYER);

  fp1=myfopen("./input/mod1.invert","r");
  fp2=myfopen("param_check.out","w");
  ReadModelParam(fp1,fp2,npar,irho,lrho,urho,
                 im,lm,um,itau,ltau,utau,ic,lc,uc);
  fclose(fp1),fclose(fp2);

  mcmc.TunePara(num);
  fp1=myfopen("./input/example.txt","r");
  mcmc.InputData(fp1,npar,irho,lrho,urho,im,lm,um,
                 itau,ltau,utau,ic,lc,uc,&idum);
  fclose(fp1);
  M=numM/saveID;
  nvar=mcmc.xsize;
  outall=fmalloc(nvar*M+1);
  misfit=dmalloc(M+1);

  TIME0
  for(i=1;i<=numM;i++) {
    Debug=((i-(i/numID)*numID)==0 ? 1:0);
    if(Debug==1) printf("Iteration:%d\n",i);
    mcmc.UpdateChain(p,&idum);
    Debug2=((i-(i/saveID)*saveID)==0 ? 1:0);
    if(Debug2==1) {
      j=i/saveID; 
      mcmc.SaveChain(j,M,outall,misfit);
    };
  };
  TIME1("Time for MCMC sampling:");
  fp1=myfopen("./output/example_mod1.bin","w");
  fp2=myfopen("./output/example_mod1.ind","w");
  fp3=myfopen("./output/example_mod1.msft","w");

  mcmc.OutputChain(fp1,fp2,fp3,M,outall,misfit);
  fclose(fp1),fclose(fp2),fclose(fp3);
  free(outall),free(misfit);
}

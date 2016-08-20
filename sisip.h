/***********************************************************
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

#include "./SubFun/basefuns.h"
#include "./SubFun/randfuns.h"

#define MAXLAYER 10
#define MAXFREQ  2000

class COLE
{  
  private:
    double *freq,*rldt,*igdt;
    long nmdl,nfrq;
    int readID;

  public:
    COLE();
    ~COLE();
    void ReadData(FILE*,long*);
    void Forward(FILE*,long*,float*,
                 double*,double*,double*,double*,
                 float,long*);
    void GetCoef(double*,double*,double*,double*,
                 double*,double*);
    double LogLkhd(double*,double*,double*,
                   double*,double*,double*);
    void SelectModel(double*,double*,double*,
                   double*,double*,double*);
};

class ColePost
{
  private:
    COLE cole;
    double *lrho,*urho,*lm,*um,*ltau,*utau,*lc,*uc;
    double mse,setup;
    void GetParaBack(double*,double*,double*,
                     double*,double*,double*);

  public:  
    ColePost();
    ~ColePost();
    long nfrq,nmdl,nvar;
    void InputData(FILE*,long*);
    void InitModel(double*,double*,double*,
                   double*,double*,double*,
                   double*,double*,double*,
                   double*,double*,double*,
                   double*,long*);
    void UpdateModel(void);
    void GetCoef(double*,double*);
    void GetBounds(double*,double*,double*,long*);
    double LogFunc(double*,double*);
    void Output(double*,long,long,float*,double*);
};


class SAMP
{
  private:
    long nmdl,num[4];
    double *x;
    ColePost post;
    void UpdateRho(long*);
    void UpdateTauAB(long*);
    void MultipleMH(int*,int,long*);
    void MultipleSS0(int*,int,long*);
    void MultipleSS1(int*,int,long*);
		       
  public:
    SAMP();
    ~SAMP();
    long xsize;
    void TunePara(long*);
    void InputData(FILE*,long*,double*,double*,double*,
                   double*,double*,double*,double*,
                   double*,double*,double*,double*,
                   double*,long*); 
    void UpdateChain(float*,long*);
    void SaveChain(long,long,float*,double*);
    void OutputChain(FILE*,FILE*,FILE*,long,float*,double*);
};

void ReadModelParam(FILE*,FILE*,long*,double*,double*,
                    double*,double*,double*,double*,
                    double*,double*,double*,double*,
                    double*,double*);

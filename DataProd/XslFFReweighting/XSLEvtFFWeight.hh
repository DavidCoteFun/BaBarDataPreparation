//*****************************************************************
//   
//   Interface class inherited by XSLPseudoScalarFF and XSLVectorFF.
//   XSLPseudoScalarFF and XSLVectorFF themselves are abstract classes from
//   which the various XSLModelName classes inherit.
//   If you add new models, please make them inherit from XSLPseudoScalarFF or XSLVectorFF
//
//   These class should be used in this fashion:
//
//   XSLEvtFFWeight* S = new XSLPseudoScalarISGW2( myXSLKin );
//   double w1,w2;
//   w1 = S->ISGW2ToThisModel();
//   w2 = S->DGammaDq2();
//   delete S;
//
//   XSLEvtFFWeight* V = new XSLVectorISGW2(myXSLKin);
//   //stuff...
//   delete V;
//
//  The 4-Vectors given has input should be in the LAB frame.
//  In case of large resonnances (rho), one should use the measured mass instead of the mean PDG mass.
//  We're following conventions from: P. Burchat and J. Richman, Rev.Mod.Phys.67:893-976,1995
//                                                               e-Print Archive: hep-ph/9508250 
//
//  NOTE: Burchat/Richman convention is NOT the same as Gilman and Singleton, Phys.Rev.D41:142,1990  

//  Creation: David Cote, Universite de Montreal, 12/11/03
//  Complete documentation available in BAD #809

#ifndef XSLEVTFFWEIGHT
#define XSLEVTFFWEIGHT

class XSLKin;

class XSLEvtFFWeight {
public:
  XSLEvtFFWeight(){ 
    PI=3.141592654; 
    G_F=1.166371e-5; //GeV^-2 * (hbar*c)^3
    Vub=3.67e-3;     //PDG 2004: |Vub|=(3.67 +/- 0.47) x 10^-3
  }
  virtual ~XSLEvtFFWeight(){};

  virtual double FromSP4ToThisModel()=0;
  virtual double FromSP5ToThisModel()=0;
  virtual double FromSP6ToThisModel()=0;
  virtual double FromSP7ToThisModel()=0;
  virtual double FromISGW2ToThisModel()=0;
  virtual double FromPHSPToThisModel()=0;
  virtual double FromFLATQ2ToThisModel()=0;
  virtual double dGammadQ2()=0;
  virtual double dGammadQ2dCtl()=0;
  virtual double dGammadQ2dCtv()=0;
  virtual double dGammadQ2dCtldCtv()=0;
  virtual double dGammadQ2dCtldCtvdChi()=0;
  virtual double Gamma()=0;  //in ps^-1

  double PI;
  double G_F;
  double Vub;
};

#endif





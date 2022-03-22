//*****************************************************************
//   
//   Creation: David Cote, 08/24/05, Universite de Montreal
//   see comments in the XSLEvtFFWeight.hh file or BAD #809
//
//   This class can produce reweighting with the functions of:
//     -R. Hill (hep-ph/0505129)        or
//     -Becirevic-Kaidalov (Phys. Lett. B478 (2000))

#include "BaBar/BaBar.hh"

#include "XslFFReweighting/XSLGenericFplus.hh"
#include "CLHEP/Vector/LorentzVector.h"

using std::cout;
using std::endl;
using std::string;

XSLGenericFplus::XSLGenericFplus(double mB,double mXu,double q2,double theta_l,string mode,double alpha,double delta) : 
  XSLPseudoScalarFF(mB,mXu,q2,theta_l) 
{  
  _alpha=alpha;
  _delta=delta;
  SetNormalizations(mode);  
  Compute(); 
}

XSLGenericFplus::XSLGenericFplus(XSLKin* DecayKin,string mode,double alpha,double delta) : XSLPseudoScalarFF(DecayKin) 
{  
  _alpha=alpha;
  _delta=delta;
  SetNormalizations(mode);  
  Compute(); 
}

XSLGenericFplus::XSLGenericFplus(HepLorentzVector BLab,HepLorentzVector LepLab,HepLorentzVector XuLab,string mode,double alpha,double delta) : 
  XSLPseudoScalarFF( BLab, LepLab, XuLab) 
{  
  _alpha=alpha;
  _delta=delta;
  SetNormalizations(mode);  
  Compute(); 
}

double
XSLGenericFplus::GetFplus(double q2)
{
  //Note: The R. Hill formula becomes like BK when delta=0
  double z = q2/5.32/5.32;
  float fp=(1.0-_delta*z)/(1-z)/(1-(_alpha+_delta*(1-_alpha))*z);
  return fp;
}


void
XSLGenericFplus::SetNormalizations(string mode)
{
  //Note 1: These normalization constants were determined empirically from MC samples of 1000k events -> 0.3% stat err.
  //Note 2: Only FLATQ2 is currently supported. Please contact cote@slac.stanford.edu is another generator is needed.
  _FLATQ2Normalization=1.0; 
  _PHSPNormalization=1.0; 
  _ISGW2Normalization=1.0;  

  if(mode=="pilnu"||mode=="pi0lnu"||mode=="pi0lnu_o") { 
    if(_alpha==0.45 && _delta==0.0){ _FLATQ2Normalization=1.0/8.722711; }
    else if(_alpha==0.52 && _delta==0.0){ _FLATQ2Normalization=1.0/9.367713; }
    else if(_alpha==0.60 && _delta==0.0){ _FLATQ2Normalization=1.0/10.29179; }
    else if(_alpha==0.61 && _delta==0.0){ _FLATQ2Normalization=1.0/10.407812; }
    else if(_alpha==0.70 && _delta==0.0){ _FLATQ2Normalization=1.0/11.770405; }
    else if(_alpha==0.75 && _delta==0.0){ _FLATQ2Normalization=1.0/12.755433; }
    else{cout<<"Unknown normalization for alpha="<<_alpha<<" and delta="<<_delta<<"...  Results won't be reliable!"<<endl;}
  }
  else if(mode=="rhoClnu"||mode=="rho0lnu"||mode=="omegalnu"||mode=="etalnu"||mode=="eta2lnu"||mode=="eta3lnu"
	  ||mode=="etaplnu"||mode=="etaplnuE2PP"||mode=="etaplnuE3PP"||mode=="etaplnuRG") {
    cout<<"XSLGenericFplus::SetNormalizations  -- This model is not intended to be used with mode: "<<mode<<" !!!!"<<endl;
    cout<<"We still fix the normalizations to 1.0, but these results might be meaningless."<<endl;
  }
  else if(mode!="Ignore"){
    cout<<"XSLGenericFplus::SetNormalizations  -- Unknown mode: "<<mode<<" !!!    Normalizations set to 1..."<<endl; 
  }
  
  return;
}

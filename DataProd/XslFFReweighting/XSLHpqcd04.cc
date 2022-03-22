//*****************************************************************
//   
//   Creation: Jochen Dingfelder 11/15/04 SLAC
//   see comments in the XSLEvtFFWeight.hh file or BAD #809
//
//   From the paper: hep-lat/0408019

#include "BaBar/BaBar.hh"

#include "XslFFReweighting/XSLHpqcd04.hh"
#include "CLHEP/Vector/LorentzVector.h"

using std::cout;
using std::endl;
using std::string;

XSLHpqcd04::XSLHpqcd04( double mB, double mXu, double q2, double theta_l,string mode) : 
  XSLPseudoScalarFF(mB,mXu,q2,theta_l) 
{  SetNormalizations(mode); Compute(); }

XSLHpqcd04::XSLHpqcd04( XSLKin* DecayKin,string mode ) : XSLPseudoScalarFF(DecayKin) 
{  SetNormalizations(mode); Compute(); }

XSLHpqcd04::XSLHpqcd04( HepLorentzVector BLab, HepLorentzVector LepLab, HepLorentzVector XuLab,string mode) : XSLPseudoScalarFF( BLab, LepLab, XuLab) 
{  SetNormalizations(mode);  Compute(); }

double
XSLHpqcd04::GetFplus(double q2)
{

  const double mBs2=5.32*5.32;
  double zp=q2/mBs2;
  
  const double cB =  0.42;
  const double alpha = 0.41; 

  double fval = cB*(1-alpha)/(1-zp)/(1-alpha*zp);
  
  return fval;
}


void
XSLHpqcd04::SetNormalizations(string mode)
{
  //Note: These normalization constants were determined empirically from MC samples of 1000k events (for most of the modes)
  //The quoted error is simply the statistical one, but there can be a bigger error coming from rare high weight events for
  //some region of the phase space, mainly for the Vector-ISGW2 generator.

  if(mode=="Ignore") { _FLATQ2Normalization=1.0; _PHSPNormalization=1.0; _ISGW2Normalization=1.0; }
  else if(mode=="pilnu") { 
    _FLATQ2Normalization=1.935984992; //determined from 100k events: 0.3% error stat
    _PHSPNormalization=2.660944804; //determined from 100k events: 0.3% error stat
    _ISGW2Normalization=1.67626031; } //determined from 100k events: 0.3% error stat
  else if(mode=="pi0lnu"||mode=="pi0lnu_o") { 
    _FLATQ2Normalization=1.932387425; //determined from 100k events: 0.3% error stat
    _PHSPNormalization=2.658378065; //determined from 100k events: 0.3% error stat
    _ISGW2Normalization=1.6738413; } //determined from 100k events: 0.3% error stat
  else if(mode=="rhoClnu"||mode=="rho0lnu"||mode=="omegalnu"||mode=="etalnu"||mode=="eta2lnu"||mode=="eta3lnu"
	  ||mode=="etaplnu"||mode=="etaplnuE2PP"||mode=="etaplnuE3PP"||mode=="etaplnuRG") {
    cout<<"XSLHpqcd04::SetNormalizations  -- This model is not intended to be used with mode: "<<mode<<" !!!!"<<endl;
    cout<<"We still fix the normalizations to 1.0, but these results might be meaningless."<<endl;
    _FLATQ2Normalization=1.0; _PHSPNormalization=1.0; _ISGW2Normalization=1.0;  }
  else {
    cout<<"XSLHpqcd04::SetNormalizations  -- Unknown mode: "<<mode<<" !!!    Normalizations set to 1..."<<endl;
    _FLATQ2Normalization=1.0; _PHSPNormalization=1.0; _ISGW2Normalization=1.0;  }
  
  return;
}


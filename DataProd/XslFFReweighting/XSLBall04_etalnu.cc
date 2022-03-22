//*****************************************************************
//   
//   Creation: David Cote, Universite de Montreal, 12/11/03
//   Source: hep-ph/0406232 
//   see comments in the XSLEvtFFWeight.hh file or BAD #809
#include "BaBar/BaBar.hh"

#include "XslFFReweighting/XSLBall04_etalnu.hh"
#include "CLHEP/Vector/LorentzVector.h"
using std::cout;
using std::endl;
using std::string;

XSLBall04_etalnu::XSLBall04_etalnu( double mB, double mXu, double q2, double theta_l, string mode) : 
  XSLPseudoScalarFF(mB,mXu,q2,theta_l) 
{  SetNormalizations(mode);  Compute(); }

XSLBall04_etalnu::XSLBall04_etalnu( XSLKin* DecayKin, string mode) : XSLPseudoScalarFF(DecayKin) 
{  SetNormalizations(mode);  Compute(); }

XSLBall04_etalnu::XSLBall04_etalnu( HepLorentzVector BLab, HepLorentzVector LepLab, HepLorentzVector XuLab, string mode) : XSLPseudoScalarFF( BLab, LepLab, XuLab) 
{  SetNormalizations(mode);  Compute(); }

double
XSLBall04_etalnu::GetFplus(double q2)
{
  // hep-ph/0406232
  // from Table 3
  const double r1 = 0.122;
  const double r2 = 0.155;
  const double m1 = 5.32; // B* mass
  
  // from Eqn (31)
  double fplus = r1/(1.-q2/m1/m1) + r2/pow((1.-q2/m1/m1),2);
  
  return fplus;
}


void
XSLBall04_etalnu::SetNormalizations(string mode)
{
  //Note: These normalization constants were determined empirically from MC samples of 1000k events (for most of the modes)
  //The quoted error is simply the statistical one, but there can be a bigger error coming from rare high weight events for
  //some region of the phase space, mainly for the Vector-ISGW2 generator.

  if(mode=="Ignore") { _FLATQ2Normalization=1.0; _PHSPNormalization=1.0; _ISGW2Normalization=1.0; }
  else if(mode=="etalnu"||mode=="eta2lnu"||mode=="eta3lnu"||mode=="etalnu_o") { 
    _FLATQ2Normalization=1.244341164; //determined from 100k events: 0.3% error stat
    _PHSPNormalization=1.844734199;  //determined from 100k events: 0.3% error stat
    _ISGW2Normalization=1.806008959; } //determined from 100k events: 0.3% error stat
  else if(mode=="etaplnu"||mode=="etaplnuE2PP"||mode=="etaplnuE3PP"||mode=="etaplnuRG"||mode=="etaplnu_o") { 
    _FLATQ2Normalization=1.596165412; //determined from 100k events: 0.3% error stat
    _PHSPNormalization=2.428592963; //determined from 100k events: 0.3% error stat
    _ISGW2Normalization=3.39687119; } //determined from 100k events: 0.3% error stat
  else if(mode=="rhoClnu"||mode=="rho0lnu"||mode=="omegalnu"||mode=="pi0lnu"||mode=="pilnu") {
    cout<<"XSLBall04_etalnu::SetNormalizations  -- This model is not intended to be used with mode: "<<mode<<" !!!!"<<endl;
    cout<<"We still fix the normalizations to 1.0, but these results might be meaningless."<<endl;
    _FLATQ2Normalization=1.0; _PHSPNormalization=1.0; _ISGW2Normalization=1.0;  }
  else {
    cout<<"XSLBall04_etalnu::SetNormalizations  -- Unknown mode: "<<mode<<" !!!    Normalizations set to 1..."<<endl;
    _FLATQ2Normalization=1.0; _PHSPNormalization=1.0; _ISGW2Normalization=1.0;  }
  
  return;
}

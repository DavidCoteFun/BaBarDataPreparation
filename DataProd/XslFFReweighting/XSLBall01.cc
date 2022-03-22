//*****************************************************************
//   
//   Creation: David Cote, Universite de Montreal, 12/11/03
//   Source: hep-ph/0110115 v1
//   see comments in the XSLEvtFFWeight.hh file or BAD #809
#include "BaBar/BaBar.hh"

#include "XslFFReweighting/XSLBall01.hh"
#include "CLHEP/Vector/LorentzVector.h"
using std::cout;
using std::endl;
using std::string;

XSLBall01::XSLBall01( double mB, double mXu, double q2, double theta_l, string mode) : 
  XSLPseudoScalarFF(mB,mXu,q2,theta_l) 
{  SetNormalizations(mode);  Compute(); }

XSLBall01::XSLBall01( XSLKin* DecayKin, string mode) : XSLPseudoScalarFF(DecayKin) 
{  SetNormalizations(mode);  Compute(); }

XSLBall01::XSLBall01( HepLorentzVector BLab, HepLorentzVector LepLab, HepLorentzVector XuLab, string mode) : XSLPseudoScalarFF( BLab, LepLab, XuLab) 
{  SetNormalizations(mode);  Compute(); }

double
XSLBall01::GetFplus(double q2)
{
  double z=q2/(5.2794*5.2794);
  double LCSR=0;
  if(q2<15.7) LCSR=(0.261/(1-(2.03*z)+(1.29*z*z)));
  else{ 
    double mBs2=28.4;
    double zp=q2/mBs2;
    LCSR=(0.439/(1-zp));
  }

  return LCSR;
}


void
XSLBall01::SetNormalizations(string mode)
{
  //Note: These normalization constants were determined empirically from MC samples of 1000k events (for most of the modes)
  //The quoted error is simply the statistical one, but there can be a bigger error coming from rare high weight events for
  //some region of the phase space, mainly for the Vector-ISGW2 generator.

  if(mode=="Ignore") { _FLATQ2Normalization=1.0; _PHSPNormalization=1.0; _ISGW2Normalization=1.0; }
  else if(mode=="pilnu") { 
    _FLATQ2Normalization=1.232487888; //determined from 1000k events: 0.1% error stat
    _PHSPNormalization=1.690646194; //determined from 1000k events: 0.1% error stat
    _ISGW2Normalization=1.066761395; } //determined from 1000k events: 0.1% error stat
  else if(mode=="pi0lnu"||mode=="pi0lnu_o") { 
    _FLATQ2Normalization=1.232202567; //determined from 1000k events: 0.1% error stat
    _PHSPNormalization=1.693948053; //determined from 1000k events: 0.1% error stat
    _ISGW2Normalization=1.064911267; } //determined from 1000k events: 0.1% error stat
  else if(mode=="etalnu"||mode=="eta2lnu"||mode=="eta3lnu"||mode=="etalnu_o") { 
    _FLATQ2Normalization=1.30865483; //determined from 1000k events: 0.1% error stat
    _PHSPNormalization=1.936813517;  //determined from 1000k events: 0.1% error stat
    _ISGW2Normalization=1.897036675; } //determined from 1000k events: 0.1% error stat
  else if(mode=="etaplnu"||mode=="etaplnuE2PP"||mode=="etaplnuE3PP"||mode=="etaplnuRG"||mode=="etaplnu_o") { 
    _FLATQ2Normalization=1.550885774; //determined from 1000k events: 0.1% error stat
    _PHSPNormalization=2.360587919; //determined from 1000k events: 0.1% error stat
    _ISGW2Normalization=3.3028835; } //determined from 1000k events: 0.1% error stat
  else if(mode=="rhoClnu"||mode=="rho0lnu"||mode=="omegalnu") {
    cout<<"XSLBall01::SetNormalizations  -- This model is not intended to be used with mode: "<<mode<<" !!!!"<<endl;
    cout<<"You probably want to use XSLBall98 instead... We still fix the normalizations to 1.0, but these results might be meaningless."<<endl;
    _FLATQ2Normalization=1.0; _PHSPNormalization=1.0; _ISGW2Normalization=1.0;  }
  else {
    cout<<"XSLBall01::SetNormalizations  -- Unknown mode: "<<mode<<" !!!    Normalizations set to 1..."<<endl;
    _FLATQ2Normalization=1.0; _PHSPNormalization=1.0; _ISGW2Normalization=1.0;  }
  
  return;
}

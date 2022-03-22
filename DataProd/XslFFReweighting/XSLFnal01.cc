//*****************************************************************
//   
//   Creation: David Cote, Universite de Montreal, 12/11/03
//   see comments in the XSLEvtFFWeight.hh file or BAD #809
//
//   A Quenched Fermilab calculation  -  hep-ph/0101023v1
#include "BaBar/BaBar.hh"


#include "XslFFReweighting/XSLFnal01.hh"
#include "CLHEP/Vector/LorentzVector.h"
using std::cout;
using std::endl;
using std::string;

XSLFnal01::XSLFnal01( double mB, double mXu, double q2, double theta_l, string mode) : 
  XSLPseudoScalarFF(mB,mXu,q2,theta_l) 
{  SetNormalizations(mode);  Compute(); }

XSLFnal01::XSLFnal01( XSLKin* DecayKin, string mode) : XSLPseudoScalarFF(DecayKin) 
{  SetNormalizations(mode);  Compute(); }

XSLFnal01::XSLFnal01( HepLorentzVector BLab, HepLorentzVector LepLab, HepLorentzVector XuLab, string mode) : XSLPseudoScalarFF( BLab, LepLab, XuLab) 
{  SetNormalizations(mode);  Compute(); }

double
XSLFnal01::GetFplus(double q2)
{
  double w=0;
  
  if(q2<16.71) w=0;
  else if(q2<17.75) w=1.13;
  else if(q2<18.79) w=1.36;
  else if(q2<19.83) w=1.59;
  else if(q2<20.865) w=1.72;
  else if(q2<21.895) w=1.84;
  else if(q2<22.91) w=1.96;
  else if(q2<23.90) w=2.10;
  else w=0;
 
  return w;
}

void
XSLFnal01::SetNormalizations(string mode)
{
  if(mode=="Ignore") { _FLATQ2Normalization=1.0; _PHSPNormalization=1.0; _ISGW2Normalization=1.0; }
  else if(mode=="pilnu") { _FLATQ2Normalization=1.0; _PHSPNormalization=1.0; _ISGW2Normalization=1.0; }
  else if(mode=="pi0lnu") { _FLATQ2Normalization=1.0; _PHSPNormalization=1.0; _ISGW2Normalization=1.0; }
  else if(mode=="etalnu"||mode=="eta2lnu"||mode=="eta3lnu") { _FLATQ2Normalization=1.0; _PHSPNormalization=1.0; _ISGW2Normalization=1.0; }
  else if(mode=="etaplnu"||mode=="etaplnuE2PP"||mode=="etaplnuE3PP"||mode=="etaplnuRG") 
    { _FLATQ2Normalization=1.0; _PHSPNormalization=1.0; _ISGW2Normalization=1.0; }
  else if(mode=="rhoClnu") { _FLATQ2Normalization=1.0; _PHSPNormalization=1.0; _ISGW2Normalization=1.0; }
  else if(mode=="rho0lnu") { _FLATQ2Normalization=1.0; _PHSPNormalization=1.0; _ISGW2Normalization=1.0; }
  else if(mode=="omegalnu") { _FLATQ2Normalization=1.0; _PHSPNormalization=1.0; _ISGW2Normalization=1.0; }
  else {
    cout<<"XSLFnal01::SetNormalizations  -- Unknown mode: "<<mode<<" !!!    Normalizations set to 1..."<<endl;
    _FLATQ2Normalization=1.0; _PHSPNormalization=1.0; _ISGW2Normalization=1.0;  }
  
  return;
}

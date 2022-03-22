//*****************************************************************
//   
//   Creation: David Cote, Universite de Montreal, 12/11/03
//   see comments in the XSLEvtFFWeight.hh file or BAD #809
//
//   From the LCSR calculation of P. Ball -  hep-ph/9805422
#include "BaBar/BaBar.hh"

#include "XslFFReweighting/XSLBall98.hh"
#include "CLHEP/Vector/LorentzVector.h"
using std::cout;
using std::endl;
using std::string;

XSLBall98::XSLBall98( double mB, double mXu, double q2, double theta_l, double theta_V, double chi, string mode) : 
  XSLVectorFF(mB,mXu,q2,theta_l,theta_V,chi) 
{  SetNormalizations(mode);  Compute(); }

XSLBall98::XSLBall98( XSLKin* DecayKin, string mode) : XSLVectorFF(DecayKin) 
{  SetNormalizations(mode);  Compute(); }

XSLBall98::XSLBall98(HepLorentzVector BLab, HepLorentzVector LepLab, HepLorentzVector XuLab, HepLorentzVector XuDaughterLab, string mode) : 
  XSLVectorFF(BLab, LepLab, XuLab, XuDaughterLab) 
{  SetNormalizations(mode);  Compute(); }

void 
XSLBall98::GetAllFF(double *A1, double *A2, double *V)
{
  *A1=GetA1(_q2);
  *A2=GetA2(_q2);
  *V=GetV(_q2);
}


double XSLBall98::GetA1(double q2)
{
  double F0=0.261;
  double aF=0.29;
  double bF=-0.415;
  double mB=5.2794;
  double mB2=mB*mB;

  double A1= F0/(1 - aF*q2/mB2 + bF*q2*q2/(mB2*mB2));
  return A1;
}

double XSLBall98::GetA2(double q2)
{
  double F0=0.223;
  double aF=0.93;
  double bF=-0.092;
  double mB=5.2794;
  double mB2=mB*mB;

  double A2= F0/(1 - aF*q2/mB2 + bF*q2*q2/(mB2*mB2));
  return A2;
}

double XSLBall98::GetV(double q2)
{ 
  double F0=0.338;
  double aF=1.37;
  double bF=0.315;
  double mB=5.2794;
  double mB2=mB*mB;

  double V= F0/(1 - aF*q2/mB2 + bF*q2*q2/(mB2*mB2));
  return V;
}


void
XSLBall98::SetNormalizations(string mode)
{
  //Note: These normalization constants were determined empirically from MC samples of 1000k events.
  //The quoted error is simply the statistical one, but there can be a bigger error coming from rare high weight events for
  //some region of the phase space, mainly for the Vector-ISGW2 generator. A simple comparison of the rho+/rho0/omega 
  //normalizations suggest a ~10% error for the ISGW2 case.

  if(mode=="Ignore") { _FLATQ2Normalization=1.0; _PHSPNormalization=1.0; _ISGW2Normalization=1.0; }
  else if(mode=="rhoClnu") {
    _FLATQ2Normalization=0.09416582451; //determined from 1000k events: 0.1% error stat
    _PHSPNormalization=0.1415359452;    //determined from 1000k events: 0.1% error stat
    _ISGW2Normalization=0.8566228565; } //determined from 1000k events: 0.1% error stat 
                                //but the effect of rare high weight event seem to introduce much larger uncertainties from ISGW2...
  else if(mode=="rho0lnu") {
    _FLATQ2Normalization=0.09426052383;  //determined from 1000k events: 0.1% error stat
    _PHSPNormalization=0.1418544045;   //determined from 1000k events: 0.1% error stat
    _ISGW2Normalization=0.7907258923; }  //determined from 1000k events: 0.1% error stat
                                //but the effect of rare high weight event seem to introduce much larger uncertainties from ISGW2...
  else if(mode=="omegalnu" || mode=="omegalnu_o") {
    _FLATQ2Normalization=0.09569658212; //determined from 1000k events: 0.1% error stat
    _PHSPNormalization=0.1446713526; //determined from 1000k events: 0.1% error stat
    _ISGW2Normalization=0.7540312395; } //determined from 1000k events: 0.1% error stat
                                //but the effect of rare high weight event seem to introduce much larger uncertainties from ISGW2...
  else if(mode=="pilnu"||mode=="pi0lnu"||mode=="etalnu"||mode=="eta2lnu"||mode=="eta3lnu"||
	  mode=="etaplnu"||mode=="etaplnuE2PP"||mode=="etaplnuE3PP"||mode=="etaplnuRG") { 
    cout<<"XSLBall98::SetNormalizations  -- This model is not intended to be used with mode: "<<mode<<" !!!!"<<endl;
    cout<<"You probably want to use XSLBall01 instead... We still fix the normalizations to 1.0, but these results might be meaningless."<<endl;
    _FLATQ2Normalization=1.0; _PHSPNormalization=1.0; _ISGW2Normalization=1.0;  }
  else {
    cout<<"XSLBall98::SetNormalizations  -- Unknown mode: "<<mode<<" !!!    Normalizations set to 1..."<<endl;
    _FLATQ2Normalization=1.0; _PHSPNormalization=1.0; _ISGW2Normalization=1.0;  }
  
  return;
}




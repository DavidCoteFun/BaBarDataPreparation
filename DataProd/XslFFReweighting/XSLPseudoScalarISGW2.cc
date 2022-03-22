//*****************************************************************
//   
//   Creation: David Cote, Universite de Montreal, 12/11/03
//   see comments in the XSLEvtFFWeight.hh file or BAD #809
//
//   Calculation from EvtGenModels/EvtISGW2.cc 
#include "BaBar/BaBar.hh"

#include "XslFFReweighting/XSLPseudoScalarISGW2.hh"
#include "CLHEP/Vector/LorentzVector.h"
#include <math.h>
using std::cout;
using std::endl;
using std::string;

XSLPseudoScalarISGW2::XSLPseudoScalarISGW2( double mB, double mXu, double q2, double theta_l, string mode) : 
  XSLPseudoScalarFF(mB,mXu,q2,theta_l) 
{  SetNormalizations(mode);  Compute();  }

XSLPseudoScalarISGW2::XSLPseudoScalarISGW2( XSLKin* DecayKin, string mode) : XSLPseudoScalarFF(DecayKin) 
{  SetNormalizations(mode);  Compute(); }

XSLPseudoScalarISGW2::XSLPseudoScalarISGW2( HepLorentzVector BLab, HepLorentzVector LepLab, HepLorentzVector XuLab, string mode) : 
  XSLPseudoScalarFF( BLab, LepLab, XuLab) 
{  SetNormalizations(mode);  Compute(); }

double
XSLPseudoScalarISGW2::GetFplus(double q2)
{
  double mass=_mXu;
  double mtb, mbb;
  double msd, mx,mb,nf,nfp; 
  double msq,bx2,mbx,mtx;
  double zji,cji,w,gammaji,chiji,betaji_fppfm;
  double rfppfm,rfpmfm,f3fppfm,f3fpmfm,fppfm,fpmfm,al,ai,rcji,f3; 
  double mqm,msb,bb2,mup,mum,bbx2,tm,wt,r2,betaji_fpmfm;
  
  double t=q2;

  msb=5.2;
  msd=0.33;
  bb2=0.431*0.431;
  mbb=5.31;
  nf = 4.0;
  
  //if (theB.contains(prnt)) {
  //if (dgt==PI0||dgt==PIP||dgt==PIM||dgt==ETA||dgt==ETAPR) {
  
  msq=0.33;
  bx2=0.406*0.406;
  mbx=0.75*0.770+0.25*0.14;
  nfp = 0.0;
  //}}
  
  mtb = msb + msd;
  mtx = msq + msd;
  //  mb=EvtPDL::getNominalMass(parent);
  mb=5.279;
  mx=mass;
  
  mup=1.0/(1.0/msq+1.0/msb);
  mum=1.0/(1.0/msq-1.0/msb);
  bbx2=0.5*(bb2+bx2);
  tm=(mb-mx)*(mb-mx);
  if ( t>tm ) t=0.99*tm;
  wt=1.0+(tm-t)/(2.0*mbb*mbx);
  
  mqm = 0.1;
  
  ///////////////////////////////////////////////
  /*  r2=3.0/(4.0*msb*msq)+3*msd*msd/(2*mbb*mbx*bbx2) + 
      (16.0/(mbb*mbx*(33.0-2.0*nfp)))*
      log(EvtGetas(mqm,mqm)/EvtGetas(msq,msq));
  */
  double lqcd2 = 0.04;
  double nflav = 4;
  double As0=0.6;
  if ( mqm > 0.6 ) {
    if ( mqm < 1.85 ) {
      nflav = 3.0;}
    
    As0 = 12.0*PI / ( 33.0 - 2.0*nflav) /
      log( mqm*mqm/lqcd2);
  }
  
  double As2=0.6;
  if ( msq > 0.6 ) {
    if ( msq < 1.85 ) {
      nflav = 3.0;}
    
    As2 = 12.0*PI / ( 33.0 - 2.0*nflav) /
      log( msq*msq/lqcd2);
  }
  r2=3.0/(4.0*msb*msq)+3*msd*msd/(2*mbb*mbx*bbx2) + 
    (16.0/(mbb*mbx*(33.0-2.0*nfp)))*
    log(As0/As2);
  //////////////////////////////////////////
  
  f3 = sqrt(mtx/mtb)*pow(sqrt(bx2*bb2)/bbx2,1.5) /
    ((1.0+r2*(tm-t)/12.0)*(1.0+r2*(tm-t)/12.0));
  
  //  for w use wt def with physical masses.
  //  report(ERROR,"EvtGen") << "before w\n";
  
  w = 1.0 + (( tm - t ) / ( 2.0* mb * mx ));
  rcji = ( 1/sqrt(w*w -1 ))*log( w + sqrt( w*w -1 ));
  al = (8.0 / ( 33.0 - 2.0*nfp ))*(w*rcji -1.0 );
  ai = -1.0* ( 6.0/( 33.0 - 2.0*nf));  
  
  ////////////////////
  //  cji = pow(( EvtGetas( msb,msb ) / EvtGetas( msq,msq ) ),ai);
  double As1=0.6;
  if ( msb > 0.6 ) {
    if ( msb < 1.85 ) {
      nflav = 3.0;}
    
    As1 = 12.0*PI / ( 33.0 - 2.0*nflav) /
      log( msb*msb/lqcd2);
  }

  
  cji = pow(( As1 / As2 ),ai);
  ///////////////////////////////////
  
  zji = msq / msb;
  
  //  gammaji = EvtGetGammaji( zji );
  gammaji = 2+((2.0*zji)/(1-zji))*log(zji);
  gammaji = -1.0*gammaji;
  ///////////////////////////
  
  chiji = -1.0 - ( gammaji / ( 1- zji ));
  betaji_fppfm = gammaji - (2.0/3.0)*chiji;
  betaji_fpmfm = gammaji + (2.0/3.0)*chiji;
  
  //  rfppfm = cji *(1.0 + betaji_fppfm*EvtGetas( msq,sqrt(msb*msq) )/PI);
  double As3=0.6;
  double m3=sqrt(msb*msq);
  
  if ( m3 > 0.6 ) {
    if ( msq < 1.85 ) {
      nflav = 3.0;}
    
    As3 = 12.0*PI / ( 33.0 - 2.0*nflav) /
      log( m3*m3/lqcd2);
  }
  rfppfm = cji *(1.0 + betaji_fppfm*As3/PI);
  
  //  rfpmfm = cji *(1.0 + betaji_fpmfm*EvtGetas( msq,sqrt(msb*msq) )/PI);  
  rfpmfm = cji *(1.0 + betaji_fpmfm*As3/PI);
  ///////////////////////////////////////////////////////////
  
  f3fppfm = f3*pow(( mbb / mtb ),-0.5)*pow((mbx/mtx),0.5);
  f3fpmfm = f3*pow(( mbb / mtb ),0.5)*pow((mbx/mtx),-0.5);
  fppfm = f3fppfm* rfppfm * ( 2.0 - ( ( mtx/msq)*(1- ( (msd*msq*bb2)
						       /(2.0*mup*mtx*bbx2)))));
  fpmfm = f3fpmfm* rfpmfm * ( mtb/msq) * ( 1 - ( ( msd*msq*bb2)/
						 ( 2.0*mup*mtx*bbx2)));
  
  double fpf = (fppfm + fpmfm)/2.0;
  //*fmf = (fppfm - fpmfm)/2.0;

  return fpf;
}


void
XSLPseudoScalarISGW2::SetNormalizations(string mode)
{
  //Note: These normalization constants were determined empirically from MC samples of 1000k events (for most of the modes)
  //The quoted error is simply the statistical one, but there can be a bigger error coming from rare high weight events for
  //some region of the phase space, mainly for the Vector-ISGW2 generator.

  if(mode=="Ignore") { _FLATQ2Normalization=1.0; _PHSPNormalization=1.0; _ISGW2Normalization=1.0; }
  else if(mode=="pilnu") { 
    _FLATQ2Normalization=1.154947463;  //determined from 1000k events: 0.1% error stat
    _PHSPNormalization=1.584143515;  //determined from 1000k events: 0.1% error stat
    _ISGW2Normalization=1.0; } //1.0 by construction
  else if(mode=="pi0lnu"||mode=="pi0lnu_o") { 
    _FLATQ2Normalization=1.157747848;  //determined from 1000k events: 0.1% error stat
    _PHSPNormalization=1.591182304;  //determined from 1000k events: 0.1% error stat
    _ISGW2Normalization=1.0; }  //1.0 by construction
  else if(mode=="etalnu"||mode=="eta2lnu"||mode=="eta3lnu"||mode=="etalnu_o") { 
    _FLATQ2Normalization=0.6896565993; //determined from 1000k events: 0.1% error stat
    _PHSPNormalization=1.020871527; //determined from 1000k events: 0.1% error stat
    _ISGW2Normalization=1.0; }  //1.0 by construction
  else if(mode=="etaplnu"||mode=="etaplnuE2PP"||mode=="etaplnuE3PP"||mode=="etaplnuRG"||mode=="etaplnu_o") { 
    _FLATQ2Normalization=0.4698929596; //determined from 1000k events: 0.1% error stat
    _PHSPNormalization=0.7150845659; //determined from 1000k events: 0.1% error stat
    _ISGW2Normalization=1.0; }  //1.0 by construction
  else if(mode=="rhoClnu"||mode=="rho0lnu"||mode=="omegalnu") {
    cout<<"XSLPseudoScalarISGW2::SetNormalizations  -- This model is not intended to be used with mode: "<<mode<<" !!!!"<<endl;
    cout<<"You probably want to use XSLVectorISGW2 instead... We still fix the normalizations to 1.0, but these results might be meaningless."<<endl;
    _FLATQ2Normalization=1.0; _PHSPNormalization=1.0; _ISGW2Normalization=1.0;  }
  else {
    cout<<"XSLPseudoScalarISGW2::SetNormalizations  -- Unknown mode: "<<mode<<" !!!    Normalizations set to 1..."<<endl;
    _FLATQ2Normalization=1.0; _PHSPNormalization=1.0; _ISGW2Normalization=1.0;  }

  
  return;
}


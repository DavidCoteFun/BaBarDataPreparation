//*****************************************************************
//   
//   Creation: David Cote, Universite de Montreal, 12/11/03
//   see comments in the XSLEvtFFWeight.hh file or BAD #809
//
//   From the LCSR calculation of P. Ball - hep-ph/0412079
//   to be published in late 04/ early 05 
//
#include "BaBar/BaBar.hh"
 
#include "XslFFReweighting/XSLBall05.hh"
#include "CLHEP/Vector/LorentzVector.h"
using std::cout;
using std::endl;
using std::string;

XSLBall05::XSLBall05(double mB,double mXu,double q2,double theta_l,double theta_V,double chi,string mode,string ErrorScale) : 
  XSLVectorFF(mB,mXu,q2,theta_l,theta_V,chi) 
{ Init(mode,ErrorScale); }

XSLBall05::XSLBall05( XSLKin* DecayKin,string mode,string ErrorScale ) : XSLVectorFF(DecayKin) 
{ Init(mode,ErrorScale); }

XSLBall05::XSLBall05(HepLorentzVector BLab,HepLorentzVector LepLab,HepLorentzVector XuLab,HepLorentzVector XuDaughterLab,
		     string mode,string ErrorScale ) :  XSLVectorFF(BLab, LepLab, XuLab, XuDaughterLab) 
{ Init(mode,ErrorScale); }

void
XSLBall05::Init(string mode, string ErrorScale)
{ 
  //default values: do nothing
  _A1ScaleFunction=&XSLBall05::ReturnOne;
  _A2ScaleFunction=&XSLBall05::ReturnOne;
  _VScaleFunction=&XSLBall05::ReturnOne;

  //Scale one of the FF by some error function...
  if(ErrorScale=="A1Pos"){ _A1ScaleFunction=&XSLBall05::ErrorScaleFactorPos; } 
  else if(ErrorScale=="A1Neg"){ _A1ScaleFunction=&XSLBall05::ErrorScaleFactorNeg; } 
  else if(ErrorScale=="A2Pos"){ _A2ScaleFunction=&XSLBall05::ErrorScaleFactorPos; } 
  else if(ErrorScale=="A2Neg"){ _A2ScaleFunction=&XSLBall05::ErrorScaleFactorNeg; } 
  else if(ErrorScale=="VPos"){ _VScaleFunction=&XSLBall05::ErrorScaleFactorPos; } 
  else if(ErrorScale=="VNeg"){ _VScaleFunction=&XSLBall05::ErrorScaleFactorNeg; } 
  else if(ErrorScale!="NoError"){ cout<<"Unkown error function is XSLBall05!!"<<endl; exit(1); }

  _ErrorScale=ErrorScale; //ErrorScale doesn't need to be a class variable, but SetNormalizations needs its value 
                          //and is required to have only one argument, so this is an easy hack.
  SetNormalizations(mode);  
  SetConstants(mode); 
  Compute(); 
}


void
XSLBall05::SetConstants(string mode){

  //By default, the constants are set to B->rholnu, but they also can be set to B->omegalnu
 
  // hep-ph/0412079
  // from Table 8 
  if(mode=="omegalnu") {
    _r2_A1=-0.217;
    _mfit2_A1=37.01;
    _r1_A2=0.006;
    _r2_A2=0.192;
    _mfit2_A2=41.24;
    _r1_V=1.006;
    _r2_V=-0.713;
    _mfit2_V=37.45;
  }
  else { //rholnu
    _r2_A1=0.240;
    _mfit2_A1=37.51;
    _r1_A2=0.009;
    _r2_A2=0.212;
    _mfit2_A2=40.82;
    _r1_V=1.045;
    _r2_V=-0.721;
    _mfit2_V=38.34;
  }
  return;
}

void 
XSLBall05::GetAllFF(double *A1, double *A2, double *V)
{
  *A1=GetA1(_q2);
  *A2=GetA2(_q2);
  *V=GetV(_q2);
}


double 
XSLBall05::GetA1(double q2)
{
 // from Eqn (61)
  double a1 = _r2_A1/(1.-q2/_mfit2_A1);
  double f=(this->*_A1ScaleFunction)(q2);
  return f*a1;
}

double 
XSLBall05::GetA2(double q2)
{
  // from Eqn (60)
  double a2 = _r1_A2/(1.-q2/_mfit2_A2) + _r2_A2/pow(1.-q2/_mfit2_A2,2.);
  double f=(this->*_A2ScaleFunction)(q2);
  return f*a2;

}

double 
XSLBall05::GetV(double q2)
{ 
  const double m1 = 5.32; // B* mass
  // from Eqn (59)
  double v = _r1_V/(1.-q2/m1/m1) + _r2_V/(1.-q2/_mfit2_V);
  double f=(this->*_VScaleFunction)(q2);
  return f*v;
}

void
XSLBall05::SetNormalizations(string mode)
{
  //Note: These normalization constants were determined empirically from MC samples of 1000k events.
  //The quoted error is simply the statistical one, but there can be a bigger error coming from rare high weight events for
  //some region of the phase space, mainly for the Vector-ISGW2 generator. 
  
  _FLATQ2Normalization=1.0; _PHSPNormalization=1.0; _ISGW2Normalization=1.0; 
  
  //WARNING: Weights won't be properly normalized if you're using the PHSP generator with ErrorScale!=NoError...
  //FLATQ2 or ISGW2 Generators are fine!    :-)

  if(mode=="rhoClnu") {
    if(_ErrorScale=="NoError"){
      _FLATQ2Normalization=0.09232327364;//determined from 100k events: 0.3% error stat
      _PHSPNormalization=0.1393226273;   //determined from 100k events: 0.3% error stat
      _ISGW2Normalization=0.8447503878;  //determined from 100k events: 0.3% error stat 
                                         //but the effect of rare high weight event seem to introduce much larger uncertainties from ISGW2...
    }
    else if(_ErrorScale=="A1Pos"){
      _FLATQ2Normalization=1.0/14.10562000;
      _ISGW2Normalization=1.0/1.5742245;
    }
    else if(_ErrorScale=="A1Neg"){
      _FLATQ2Normalization=1.0/7.443169500;
      _ISGW2Normalization=1.0/0.8317134375;
    }
    else if(_ErrorScale=="A2Pos"){
      _FLATQ2Normalization=1.0/9.222988000;
      _ISGW2Normalization=1.0/1.03149775;
    }
    else if(_ErrorScale=="A2Neg"){
      _FLATQ2Normalization=1.0/11.777671;
      _ISGW2Normalization=1.0/1.313077375;
    }
    else if(_ErrorScale=="VPos"){
      _FLATQ2Normalization=1.0/10.935409;
      _ISGW2Normalization=1.0/1.219753;
    }
    else if(_ErrorScale=="VNeg"){
      _FLATQ2Normalization=1.0/9.928505;
      _ISGW2Normalization=1.0/1.110209375;
    }
    else{ cout<<"XSLBall05: ErrorScale: "<<_ErrorScale<<" not implemented for mode "<<mode<<endl; exit(1); }
  }
  else if(mode=="rho0lnu") {
    if(_ErrorScale=="NoError"){
      _FLATQ2Normalization=0.09248629758;//determined from 100k events: 0.3% error stat
      _PHSPNormalization=0.138725886;    //determined from 100k events: 0.3% error stat
      _ISGW2Normalization=0.8443531312;  //determined from 100k events: 0.3% error stat
                                         //but the effect of rare high weight event seem to introduce much larger uncertainties from ISGW2...
    }
    else if(_ErrorScale=="A1Pos"){
      _FLATQ2Normalization=1.0/14.08362064;
      _ISGW2Normalization=1.0/1.577346125;
    }
    else if(_ErrorScale=="A1Neg"){
      _FLATQ2Normalization=1.0/7.423538135;
      _ISGW2Normalization=1.0/0.8315596636;
    }
    else if(_ErrorScale=="A2Pos"){
      _FLATQ2Normalization=1.0/9.200457279;
      _ISGW2Normalization=1.0/1.031787894;
    }
    else if(_ErrorScale=="A2Neg"){
      _FLATQ2Normalization=1.0/11.75889457;
      _ISGW2Normalization=1.0/1.315351717;
    }
    else if(_ErrorScale=="VPos"){
      _FLATQ2Normalization=1.0/10.91164425;
      _ISGW2Normalization=1.0/1.220436635;
    }
    else if(_ErrorScale=="VNeg"){
      _FLATQ2Normalization=1.0/9.910329046;
      _ISGW2Normalization=1.0/1.111868493;
    }
    else{ cout<<"XSLBall05: ErrorScale: "<<_ErrorScale<<" not implemented for mode "<<mode<<endl; exit(1); }
  }
  else if(mode=="omegalnu" || mode=="omegalnu_o") {
    if(_ErrorScale=="NoError"){ 
      _FLATQ2Normalization=0.01858994868;//determined from 100k events: 0.3% error stat
      _PHSPNormalization=0.02787669746;  //determined from 100k events: 0.3% error stat
      _ISGW2Normalization=0.05927929566; //determined from 100k events: 0.3% error stat
                                         //but the effect of rare high weight event seem to introduce much larger uncertainties from ISGW2...
    }
    else{ cout<<"XSLBall05: ErrorScale: "<<_ErrorScale<<" not implemented for mode "<<mode<<endl; exit(1); }
  }
  else if(mode=="pilnu"||mode=="pi0lnu"||mode=="etalnu"||mode=="eta2lnu"||mode=="eta3lnu"||
	  mode=="etaplnu"||mode=="etaplnuE2PP"||mode=="etaplnuE3PP"||mode=="etaplnuRG") { 
    cout<<"XSLBall05::SetNormalizations  -- This model is not intended to be used with mode: "<<mode<<" !!!!"<<endl; 
    exit(1);
  }
  else if(mode!="Ignore") {
    cout<<"XSLBall05::SetNormalizations  -- Unknown mode: "<<mode<<" !!!"<<endl; 
    exit(1);
  }
  
  return;
}


double 
XSLBall05::ErrorMethodA(double q2)
{
  //The error is 10% @ q2==0, 11.5% @ q2==7, 13% @ q2==14, and the extrapolation continues linearly for higher q2.
  //see: http://babar-hn.slac.stanford.edu:5090/HyperNews/get/semi_lept_decays/380.html
  float err=0.1+0.03*q2/14.0; 
  return err;
}




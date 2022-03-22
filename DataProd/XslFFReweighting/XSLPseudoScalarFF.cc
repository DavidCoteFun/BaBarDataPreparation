//*****************************************************************
//   
//   Creation: David Cote, Universite de Montreal, 12/11/03
//   see comments in the XSLEvtFFWeight.hh file or BAD #809
#include "BaBar/BaBar.hh"

#include "XslFFReweighting/XSLPseudoScalarFF.hh"
#include "XslFFReweighting/XSLKin.hh"
#include "XslFFReweighting/XSLPseudoScalarISGW2.hh"

#include "CLHEP/Vector/LorentzVector.h"


XSLPseudoScalarFF::XSLPseudoScalarFF( double mB, double mXu, double q2, double theta_l ) :
  XSLEvtFFWeight()
{
  _q2= q2;
  _thL = theta_l;
  
  _mXu= mXu;
  _mXu2=_mXu*_mXu;
  _mB=mB;

  _pXu = sqrt(pow((_mB*_mB+_mXu2-_q2)/(2*_mB),2)-_mXu2);
}

XSLPseudoScalarFF::XSLPseudoScalarFF( HepLorentzVector BLab, HepLorentzVector LepLab, HepLorentzVector XuLab) :
  XSLEvtFFWeight()
{
  XSLKin* MyDecayKin= new XSLKin(BLab,LepLab,XuLab);
  InitWithKin( MyDecayKin );
  delete MyDecayKin;
}

XSLPseudoScalarFF::XSLPseudoScalarFF( XSLKin* DecayKin ): XSLEvtFFWeight()
{
  InitWithKin( DecayKin );
}

void
XSLPseudoScalarFF::InitWithKin( XSLKin* MyDecayKin )
{
  _q2= MyDecayKin->q2();
  _thL = MyDecayKin->theta_l();
  
  _mXu= MyDecayKin->XuLab().m();
  _mXu2=_mXu*_mXu;
  _mB=MyDecayKin->BLab().m();

  _pXu = sqrt(pow((_mB*_mB+_mXu2-_q2)/(2*_mB),2)-_mXu2);
}

XSLPseudoScalarFF::~XSLPseudoScalarFF(){}

void
XSLPseudoScalarFF::Compute()
{
  _Fplus=GetFplus();  
  return;
}


//Event Weight from specific Generators
double 
XSLPseudoScalarFF::FromISGW2ToThisModel()
{
  double w=1.0;
 
  XSLPseudoScalarFF* ISGW2 = new XSLPseudoScalarISGW2(_mB,_mXu,_q2,_thL); 
  double denom=ISGW2->Fplus2();  
  delete ISGW2;

  w= Fplus2()/denom;
  w*= _ISGW2Normalization;

  return w;
}

double 
XSLPseudoScalarFF::FromPHSPToThisModel()
{
  double p2=_pXu*_pXu;
  double ctl = cos(_thL);
  double stl2=1-ctl*ctl;

  double w=_PHSPNormalization*p2*stl2*Fplus2();
  return w;
}

double 
XSLPseudoScalarFF::FromFLATQ2ToThisModel()
{
  double p3=_pXu*_pXu*_pXu;
  double ctl = cos(_thL);
  double stl2=1-ctl*ctl;
  double w= _FLATQ2Normalization*p3*stl2*Fplus2();
  return w;
}

//Differential Decay Rate, including proper physics units
double 
XSLPseudoScalarFF::dGammadQ2dCtl()
{
  //See XSLEvtFFWeight for the assumed value of |Vub|
  double C = Vub*Vub*G_F*G_F/32.0/PI/PI/PI;

  double p3=_pXu*_pXu*_pXu;
  double ctl = cos(_thL);
  double stl2=1-ctl*ctl;

  return C*p3*stl2*Fplus2();
}

double 
XSLPseudoScalarFF::dGammadQ2()
{
  //See XSLEvtFFWeight for the assumed value of |Vub|
  double C = Vub*Vub*G_F*G_F/24.0/PI/PI/PI;
  double p3=_pXu*_pXu*_pXu;
  
  return C*p3*Fplus2();
}

double 
XSLPseudoScalarFF::dGammadQ2(double q2)
{
  //See XSLEvtFFWeight for the assumed value of |Vub|
  double C = Vub*Vub*G_F*G_F/24.0/PI/PI/PI;

  double pXu = sqrt(pow((_mB*_mB+_mXu2-q2)/(2*_mB),2)-_mXu2);
  double p3=pXu*pXu*pXu;
  
  double fp=GetFplus(q2); // f+(q2) has no physics units, I think we need

  return C*p3*fp*fp;
}

double
XSLPseudoScalarFF::Gamma()
{
  //See XSLEvtFFWeight for the assumed value of |Vub|

  //double mLep=0.1056584;       //muon
  double mLep=0.000510999;   //electron

  double q2Min=mLep*mLep;
  double q2Max=(_mB-_mXu)*(_mB-_mXu);
  double ans=0;
  int nStep=1000;
  double nStepDouble = double(nStep);

  double stepSize=(q2Max-q2Min)/nStepDouble;
  double q2=q2Min+(stepSize/2.0); //to be in the center of the bin
  while(q2<q2Max){
    ans+=dGammadQ2(q2)*stepSize*1e12; //The 1e12 is factor to express the GeV in ps-1 instead of seconds
    q2+=stepSize;
    //if(dGammadQ2(q2)>1.0) {cout<<"dGam: "<<dGammadQ2(q2)<<"    q2: "<<q2<<endl; }
  }
  
  return ans;
}





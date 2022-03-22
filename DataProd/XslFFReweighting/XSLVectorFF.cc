//*****************************************************************
//   
//   Creation: David Cote, Universite de Montreal, 12/11/03
//   see comments in the XSLEvtFFWeight.hh file or BAD #809
#include "BaBar/BaBar.hh"

#include "XslFFReweighting/XSLVectorFF.hh"
#include "XslFFReweighting/XSLVectorISGW2.hh"
#include "XslFFReweighting/XSLKin.hh"

#include <assert.h>
#include "CLHEP/Vector/LorentzVector.h"

using std::cout;
using std::endl;


XSLVectorFF::XSLVectorFF( HepLorentzVector BLab, HepLorentzVector LepLab, HepLorentzVector XuLab, HepLorentzVector XuDaughterLab) :
  XSLEvtFFWeight()
{
  XSLKin* MyDecayKin= new XSLKin(BLab,LepLab,XuLab,XuDaughterLab);
  InitWithKin(MyDecayKin);  
  delete MyDecayKin;
}

XSLVectorFF::XSLVectorFF( XSLKin* DecayKin ): XSLEvtFFWeight()
{
  //We haven't required the  constructor to accept only vector XSLKin objects,
  //but we're still gonna make sure it is the case!
  if(!DecayKin->isVector()){ 
    cout<<"An object of type XSLVectorKin is required by XSLVectorFF!!"<<endl;
    assert(0);
  }
  
  InitWithKin(DecayKin);
}

XSLVectorFF::XSLVectorFF( double mB, double mXu, double q2, double theta_l, double theta_V, double chi ) :
  XSLEvtFFWeight()
{
  _mB = mB;
  _mB2=_mB*_mB;
  _q2 = q2;
  _y=_q2/_mB2;

  _mV = mXu;
  _mV2=_mV*_mV;
  _kV = sqrt(pow((_mB2+_mV2-_q2)/(2*_mB),2)-_mV2);

  _thL=theta_l;
  _thV=theta_V;
  _chi=chi;

  _Hplus=-666;
  _Hminus=-666;
  _Hzero=-666;
  _A1=-666;
  _A2=-666;
  _V=-666;
}

XSLVectorFF::~XSLVectorFF(){};

void
XSLVectorFF::Compute()
{
  GetAllFF(&_A1,&_A2,&_V);   
  GetHelicityAmp();
  
  return;
}

void
XSLVectorFF::InitWithKin( XSLKin* MyDecayKin )
{
  _mB = MyDecayKin->BLab().m();
  _mB2=_mB*_mB;
  _q2 = MyDecayKin->q2();
  _y=_q2/_mB2;

  _mV = MyDecayKin->XuLab().m();
  _mV2=_mV*_mV;
  _kV = sqrt(pow((_mB2+_mV2-_q2)/(2*_mB),2)-_mV2);

  _thL=MyDecayKin->theta_l();
  _thV=MyDecayKin->theta_v();
  _chi=MyDecayKin->chi();

  _Hplus=-666;
  _Hminus=-666;
  _Hzero=-666;
  _A1=-666;
  _A2=-666;
  _V=-666;
}



//Event Weight from specific Generators
double 
XSLVectorFF::FromISGW2ToThisModel()
{
  double w=1.0;

  XSLVectorFF* ISGW2 = new XSLVectorISGW2(_mB, _mV, _q2,_thL,_thV,_chi); 
  double denom=ISGW2->AngularPartOfFullWeight();
  delete ISGW2;

  w= AngularPartOfFullWeight()/denom;  
  w*= _ISGW2Normalization;

  return w;
}

double 
XSLVectorFF::FromPHSPToThisModel()
{
  double w = _y*AngularPartOfFullWeight();
  w*= _PHSPNormalization;
  return w;
}

double 
XSLVectorFF::FromFLATQ2ToThisModel()
{
  double w = _kV*_y*AngularPartOfFullWeight();
  w*= _FLATQ2Normalization;
  return w;
}

double 
XSLVectorFF::AngularPartOfFullWeight()
{  
  //NOTE: we use stl instead of |stl| (and the same for th_V). 
  //This imply difference in weights for a particular event but give the same integrated distributions in the end.

  double ctl = cos(_thL);
  double stl2=1-ctl*ctl;
  //double stl = sqrt(stl2);
  double stl=sin(_thL);

  double ctv= cos(_thV);
  double ctv2= ctv*ctv;
  double stv2=1-ctv2;
  //double stv=sqrt(stv2);
  double stv=sin(_thV);

  double eta=1.0;

  double pp=(1-eta*ctl)*(1-eta*ctl)*stv2;
  double mm=(1+eta*ctl)*(1+eta*ctl)*stv2;
  double zz=4*stl2*ctv2;
  double pz=4*eta*stl*(1-eta*ctl)*stv*ctv*cos(_chi);
  double mz=4*eta*stl*(1+eta*ctl)*stv*ctv*cos(_chi);
  double pm=2*stl2*stv2*cos(2*_chi);

  double w=pp*Hplus2() + mm*Hminus2() + zz*Hzero2() - pz*Hplus()*Hzero() + mz*Hminus()*Hzero() - pm*Hplus()*Hminus();
  
  return w;
}

//Helicity Amplitudes
void XSLVectorFF::GetHelicityAmp()
{
  _Hplus = Hplus(_A1,_V,_kV);
  _Hminus = Hminus(_A1,_V,_kV);
  _Hzero = Hzero(_A1,_A2,_q2);

  return;
}

double XSLVectorFF::Hplus(double A1, double V, double kV)
{
  return (_mB+_mV)*A1-(2*_mB*kV/(_mB+_mV))*V;
}

double XSLVectorFF::Hminus(double A1, double V, double kV)
{
  return (_mB+_mV)*A1+(2*_mB*kV/(_mB+_mV))*V;
}

double XSLVectorFF::Hzero(double A1, double A2, double q2)
{
  double kV=sqrt(pow((_mB2+_mV2-q2)/(2*_mB),2)-_mV2);
  double aa=(_mB2-_mV2-q2)*(_mB+_mV);
  double bb=4*_mB2*kV*kV/(_mB+_mV);

  return (aa*A1-bb*A2)/(2*_mV*sqrt(fabs(q2)));
}

//Differential Decay Rate
double 
XSLVectorFF::dGammadQ2dCtldCtvdChi()
{
  //note: integral (1-cos)^2 dcos = integral (1+cos)^2 dcos = 8/3
  //      integral cos^2 dcos = 2/3
  //      integral sin^2 dcos = 4/3
  //      integral dchi = 2*PI
  //      One can "easily" derive the integrated versions of this formula with these! :-)


  //See XSLEvtFFWeight for the assumed value of |Vub|
  double piFac = pow((4*PI),4);
  double C = Vub*Vub*G_F*G_F*3.0/8.0/piFac;

  double w = C*_kV*_y*AngularPartOfFullWeight();
  return w;
}

double 
XSLVectorFF::dGammadQ2dCtldCtv()
{
  //See XSLEvtFFWeight for the assumed value of |Vub|
  double C = Vub*Vub*G_F*G_F*3.0/8.0/128.0/PI/PI/PI;

  double ctl = cos(_thL);
  double stl2=1-ctl*ctl;

  double ctv= cos(_thV);
  double ctv2= ctv*ctv;
  double stv2=1-ctv2;

  double eta=1.0;
  double pp=(1-eta*ctl)*(1-eta*ctl)*stv2;
  double mm=(1+eta*ctl)*(1+eta*ctl)*stv2;
  double zz=4*stl2*ctv2;

  double w=C*_kV*_y*(pp*Hplus2() + mm*Hminus2() + zz*Hzero2());
  return w;
}


double 
XSLVectorFF::dGammadQ2dCtl()
{
  //See XSLEvtFFWeight for the assumed value of |Vub|
  double C = Vub*Vub*G_F*G_F/128.0/PI/PI/PI;

  double ctl = cos(_thL);
  double stl2=1-ctl*ctl;

  double eta=1.0;
  double minus=(1+eta*ctl)*(1+eta*ctl)/2;
  double plus=(1-eta*ctl)*(1-eta*ctl)/2;

  double w=C*_kV*_y*(plus*Hplus2()+minus*Hminus2()+stl2*Hzero2());
  return w;
}

double 
XSLVectorFF::dGammadQ2dCtv()
{
  //See XSLEvtFFWeight for the assumed value of |Vub|
  double C = Vub*Vub*G_F*G_F/128.0/PI/PI/PI;

  double ctv2 = cos(_thV)*cos(_thV);
  double stv2=1-ctv2;

  double w=C*(_kV*_y)*(stv2*Hplus2()+stv2*Hminus2()+2*ctv2*Hzero2());
  return w;
}

double 
XSLVectorFF::dGammadQ2()
{
  //See XSLEvtFFWeight for the assumed value of |Vub|
  double C = Vub*Vub*G_F*G_F/96.0/PI/PI/PI;

  double w = C*_kV*_y*(Hplus2()+Hminus2()+Hzero2());
 
  return w;
}

double 
XSLVectorFF::dGammadQ2(double q2)
{
  //See XSLEvtFFWeight for the assumed value of |Vub|
  double kV=sqrt(pow((_mB2+_mV2-q2)/(2*_mB),2)-_mV2);
  double y=q2/_mB2;

  double A1=GetA1(q2);
  double A2=GetA2(q2);
  double V=GetV(q2);
  double Hp = Hplus(A1,V,kV);
  double Hm = Hminus(A1,V,kV);
  double Hz = Hzero(A1,A2,q2);

  double C = Vub*Vub*G_F*G_F/96.0/PI/PI/PI;

  double w = C*kV*y*(Hp*Hp + Hm*Hm + Hz*Hz);

  return w;
}


double
XSLVectorFF::Gamma()
{
  //See XSLEvtFFWeight for the assumed value of |Vub|

  //double mLep=0.1056584;       //muon
  double mLep=0.000510999;   //electron

  double q2Min=mLep*mLep;
  double q2Max=(_mB-_mV)*(_mB-_mV);
  double ans=0;
  int nStep=1000;
  double nStepDouble = double(nStep);

  double stepSize=(q2Max-q2Min)/nStepDouble;
  double q2=q2Min+(stepSize/2.0); //to be in the center of the bin
  while(q2<q2Max){
    ans+=dGammadQ2(q2)*stepSize*1e12; //The 1e12 is factor to express the GeV in ps-1 instead of seconds
    q2+=stepSize;
    if(dGammadQ2(q2)>1.0) {cout<<"dGam: "<<dGammadQ2(q2)<<"    q2: "<<q2<<endl; }
  }
  
  return ans;
}


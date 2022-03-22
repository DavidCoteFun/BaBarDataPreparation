//*****************************************************************
//   
//   Creation: David Cote, Universite de Montreal, 12/11/03
//   see comments in the XSLEvtFFWeight.hh file or BAD #809
#include "BaBar/BaBar.hh"

#include "XslFFReweighting/XSLBToDstrlnu_DstrToDgam_LinearQ2.hh"

#include "XslFFReweighting/XSLKin.hh"

#include "CLHEP/Vector/LorentzVector.h"
#include <assert.h>

using std::cout;
using std::endl;

XSLBToDstrlnu_DstrToDgam_LinearQ2::XSLBToDstrlnu_DstrToDgam_LinearQ2( double mB, double mDStar, double q2, double ctl, double ctv, double chi, 
								    double R1, double R2, double rho2 ) : 
  XSLEvtFFWeight()
{
  _R1=R1;
  _R2=R2;
  _rho2=rho2;

  _mB=mB;
  _q2 = q2;
  _mDStar = mDStar;
  _w = (_mB*_mB+_mDStar*_mDStar-_q2)/(2*_mB*_mDStar);

  _ctl=ctl;
  _ctv=ctv;
  _chi=chi;
}

XSLBToDstrlnu_DstrToDgam_LinearQ2::XSLBToDstrlnu_DstrToDgam_LinearQ2( HepLorentzVector BLab, HepLorentzVector LepLab, 
								    HepLorentzVector DStarLab, HepLorentzVector DStarDaugLab,
								    double R1, double R2, double rho2 ) : 
  
  XSLEvtFFWeight()
{
  _R1=R1;
  _R2=R2;
  _rho2=rho2;

  XSLKin* MyDecayKin= new XSLKin(BLab,LepLab,DStarLab,DStarDaugLab);
  InitWithKin(MyDecayKin);  
  delete MyDecayKin;
}

XSLBToDstrlnu_DstrToDgam_LinearQ2::XSLBToDstrlnu_DstrToDgam_LinearQ2( XSLKin* DecayKin, double R1, double R2, double rho2 ): XSLEvtFFWeight()
{
  _R1=R1;
  _R2=R2;
  _rho2=rho2;

  //We haven't required the  constructor to accept only vector XSLKin objects,
  //but we're still gonna make sure it is the case!
  if(!DecayKin->isVector()){ 
    cout<<"An object of type XSLVectorKin is required by XSLVectorFF!!"<<endl;
    assert(0);
  }
  
  InitWithKin(DecayKin);
}

void
XSLBToDstrlnu_DstrToDgam_LinearQ2::InitWithKin( XSLKin* MyDecayKin )
{
  _mB = MyDecayKin->BLab().m();
  _q2 = MyDecayKin->q2();
  _mDStar = MyDecayKin->XuLab().m();
  _w = MyDecayKin->w();
  
  _ctl=MyDecayKin->ctl();
  _ctv=MyDecayKin->ctv();
  _chi=MyDecayKin->chi();
  
  return;
}

XSLBToDstrlnu_DstrToDgam_LinearQ2::~XSLBToDstrlnu_DstrToDgam_LinearQ2()
{}

double 
XSLBToDstrlnu_DstrToDgam_LinearQ2::FromSP4ToThisModel()
{
  float R1=1.18,R2=0.72,rho2=0.92;
  XSLEvtFFWeight* D = new XSLBToDstrlnu_DstrToDgam_LinearQ2(_mB, _mDStar, _q2, _ctl, _ctv, _chi, R1, R2, rho2);
  double denom = D->dGammadQ2dCtldCtvdChi();
  delete D;
  
  double w=dGammadQ2dCtldCtvdChi()/denom;
  
  return w;
}

double 
XSLBToDstrlnu_DstrToDgam_LinearQ2::FromSP7ToThisModel()
{
  float R1=1.33,R2=0.92,rho2=0.77;
  XSLEvtFFWeight* D = new XSLBToDstrlnu_DstrToDgam_LinearQ2(_mB, _mDStar, _q2, _ctl, _ctv, _chi, R1, R2, rho2);
  double denom = D->dGammadQ2dCtldCtvdChi();
  delete D;
  
  double w=dGammadQ2dCtldCtvdChi()/denom;
  
  return w;
}

//Event Weight from specific Generators
double XSLBToDstrlnu_DstrToDgam_LinearQ2::FromISGW2ToThisModel()
{ cout<<"XSLBToDstrlnu_DstrToDgam_LinearQ2::FromISGW2ToThisModel not implemented!   returning 1.0"<<endl; return 1.0; }

double XSLBToDstrlnu_DstrToDgam_LinearQ2::FromPHSPToThisModel()
{ cout<<"XSLBToDstrlnu_DstrToDgam_LinearQ2::FromPHSPToThisModel not implemented!   returning 1.0"<<endl; return 1.0; }

double XSLBToDstrlnu_DstrToDgam_LinearQ2::FromFLATQ2ToThisModel()
{ cout<<"XSLBToDstrlnu_DstrToDgam_LinearQ2::FromFLATQ2ToThisModel not implemented!   returning 1.0"<<endl; return 1.0; }

//Differential Decay Rates
double XSLBToDstrlnu_DstrToDgam_LinearQ2::dGammadQ2()
{cout<<"XSLBToDstrlnu_DstrToDgam_LinearQ2::dGammadQ2 not implemented!   returning 0"<<endl; return 0;}

double XSLBToDstrlnu_DstrToDgam_LinearQ2::dGammadQ2dCtl()
{cout<<"XSLBToDstrlnu_DstrToDgam_LinearQ2::dGammadQ2dCtl not implemented!   returning 0"<<endl; return 0;}

double XSLBToDstrlnu_DstrToDgam_LinearQ2::dGammadQ2dCtv()
{cout<<"XSLBToDstrlnu_DstrToDgam_LinearQ2::dGammadQ2dCtv not implemented!   returning 0"<<endl; return 0;}

double
XSLBToDstrlnu_DstrToDgam_LinearQ2::dGammadQ2dCtldCtv()
{cout<<"XSLBToDstrlnu_DstrToDgam_LinearQ2::dGammadQ2dCtldCtv not implemented!   returning 0"<<endl; return 0;}

double XSLBToDstrlnu_DstrToDgam_LinearQ2::Gamma()
{cout<<"XSLBToDstrlnu_DstrToDgam_LinearQ2::Gamma not implemented!   returning 0"<<endl; return 0;}

double
XSLBToDstrlnu_DstrToDgam_LinearQ2::dGammadQ2dCtldCtvdChi()
{
  double w=_w;  double ctl=_ctl;  double ctd=_ctv;  double chi=_chi;
  double rate=0.0;
    
  double pDStar=_mDStar*sqrt(w*w-1);
  double opw2=(w+1)*(w+1);
  double r=_mDStar/_mB;
  double r2=r*r;
  double omr2=(1-r)*(1-r);
  
  double hA1p=(1-rho2()*(w-1));
  
  double hpfac=(1-sqrt((w-1)/(w+1))*R1());
  double hmfac=(1+sqrt((w-1)/(w+1))*R1());
  double hzfac=(1+((w-1)/(1-r))*(1-R2()));
  
  double qfac=1-2*w*r+r2;
  if(qfac<0) qfac=0;
  
  double hp=sqrt(qfac/omr2)*hpfac;
  double hm=sqrt(qfac/omr2)*hmfac;
  double hz=hzfac;
  
  double hphmTerm = hp*hm*gPM(ctl,ctd,chi);
  double hphzTerm = hp*hz*gP0(ctl,ctd,chi);
  double hmhzTerm = hm*hz*gM0(ctl,ctd,chi);
  double hp2Term  = hp*hp*gPP(ctl,ctd,chi);
  double hm2Term  = hm*hm*gMM(ctl,ctd,chi);
  double hz2Term  = hz*hz*g00(ctl,ctd,chi);
  
  rate=hp2Term+hm2Term+hz2Term+hphmTerm+hphzTerm+hmhzTerm;
  rate*=hA1p*hA1p*opw2*pDStar;
  rate/=norm();

  return rate;
}


double
XSLBToDstrlnu_DstrToDgam_LinearQ2::norm()
{
  //parameter combos
  double r1=R1();
  double r1sq=r1*r1;
  double r2=R2();
  double r2sq=r2*r2;
  double rhosq=rho2();
  double rhosqsq=rhosq*rhosq;
  
  
  //analytic normalization from maple
  double temp=
    11.90848681 
    - 7.734921707*rhosq 
    + .2849279646*r1sq
    - 6.422176045*r2 
    + 1.449766988*rhosqsq  
    + 1.204216628*r2sq
    - .962930182*rhosq*r2sq  
    + .201215453*rhosqsq*r2sq
    + 4.908824973*rhosq*r2 
    - .998749412*rhosqsq*r2
    + .02670941074*r1sq*rhosqsq
    - .1637048251*r1sq*rhosq;
  
  //DC:
  //The additional 4.20825 is added to normalize the differential decay rate to 1.0.
  //That factor was detemined
  //empirically from 1000000 events and thus has a statistical uncertainty of ~0.3%. However,
  //that additional constant has no R1,R2,rho2 dependence and is thus cancelling exactly when
  //reweighting with a ratio of linear parametrisation like is the case of FromSP4ToThisModel().
  return temp*4.20825;
 
}











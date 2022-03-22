//*****************************************************************
//   
//   Creation: David Cote, Universite de Montreal, 02/02/04
//
#include "BaBar/BaBar.hh"

#include "XslTools/XSLYAveBooster.hh"

#include "PDT/Pdt.hh"
#include "AbsEnv/AbsEnv.hh"
#include "PepEnv/PepEnv.hh"
#include "PepEnvData/PepBeams.hh"

#include "XslFFReweighting/XSLKin.hh"
#include "XslTools/XSLRotateAndBoost.hh"

YAveBooster::YAveBooster( YAveBooster & YAve )
{
  _Bchrg = YAve.Bchrg();
  _LepLab = HepLorentzVector( YAve.LepLab() );
  _XuLab = HepLorentzVector( YAve.XuLab() );
  Init();
}

YAveBooster::YAveBooster(HepLorentzVector LepLab, HepLorentzVector XuLab, int Bchrg)
{
  _Bchrg = Bchrg;
  _LepLab = HepLorentzVector(LepLab);
  _XuLab = HepLorentzVector(XuLab);
  Init();
}


void
YAveBooster::Init()
{
  //getting PepBeams parameters...
  _Ups = gblEnv->getPep()->pepBeams()->total4Momentum();

  //boosting to the Ups(4s) frame
  XSLRotateAndBoost rAndB = XSLRotateAndBoost();
  HepLorentzVector XuUps = rAndB.BoostToFrame( _XuLab, _Ups );
  HepLorentzVector LepUps = rAndB.BoostToFrame( _LepLab, _Ups );
  HepLorentzVector YUps = XuUps + LepUps;
  
  //computing Y-frame coordinate system
  _zHatY = YUps.vect().unit();
  _yHatY = (XuUps.vect().cross(_zHatY)).unit();
  _xHatY = _yHatY.cross(_zHatY);


  //B meson 
  _mB=0;
  if(_Bchrg==0) {
    static const PdtEntry* PdtB0 = Pdt::lookup(PdtLund::B0);
    _mB = PdtB0->mass(); }
  else {
    static const PdtEntry* PdtBp = Pdt::lookup(PdtLund::B_plus);
    _mB = PdtBp->mass();
  }
  
  double eBeamCM = sqrt(_Ups.mag2())/2.; 
  
  double f=1.0;
  if( eBeamCM<_mB ) {
    // 5.2895 is the mean run1-4 CM beam energy
    // see: http://babar-hn.slac.stanford.edu:5090/HyperNews/get/semi_lept_decays/177/2/1.html
    f=5.2895/eBeamCM;
  }

  //the f is for the offres case. 
  //(This is what's done in the skim, in XSLRecoYAnalyzer::cosBY, XSLRecoYAnalyzer::pBCM and YAveBooster::Init)
  double eBeam = f*eBeamCM;
  double mY = f*YUps.m();
  double pY = f*YUps.rho();
  double eY=sqrt(mY*mY+pY*pY);

  _pBCM=sqrt( eBeam*eBeam-_mB*_mB );  

  _cosBY = (2*eY*eBeam-_mB*_mB-mY*mY)/(2*_pBCM*pY);
  _OK = ( fabs(_cosBY)<=1 ) ? true : false;
  if(_OK) _sinBY = sqrt(1-_cosBY*_cosBY);
  else _sinBY = -666;

  //Initializing YAverage Frame
  double n = 4; //this is a double to accomodate HepLorentzVector::operator /= (double), but should really be an int number
  static const double PI = 3.141592654;
  _phiMin = -PI;
  _phiMax = PI;
  _dPhi = (_phiMax - _phiMin)/n;

}

YAveBooster::~YAveBooster()
{}

HepLorentzVector
YAveBooster::BYAveLab() 
{ 
  if(!_OK) { return HepLorentzVector(0,0,0,0); } 

  double wsum=0.0;		// sum of weights 
  
  //Declaring quantities of interest...
  Hep3Vector BYAveVec = Hep3Vector(0,0,0);
  
  //Here we go!
  double phi=_phiMin;
  while(phi<_phiMax-0.0001) {
    
    //compute the B with this phi
    double pBx=_pBCM*_sinBY*cos(phi);
    double pBy=_pBCM*_sinBY*sin(phi);
    double pBz=_pBCM*_cosBY;
    Hep3Vector p3B=pBx*_xHatY+pBy*_yHatY+pBz*_zHatY;
    HepLorentzVector p4B = HepLorentzVector(0,0,0,0);
    p4B.setVectM(p3B,_mB);

    //weight and sum
    double cb=p3B.cosTheta(); //Note that we're in Ups(4S) coordinates! ;-)
    double sb2=1.0-cb*cb;
    double wt=sb2;  
    wsum+=wt;

    //Compute Stuff you're interested in
    BYAveVec += wt*p3B;

    //increment phi
    phi+=_dPhi;

  }//while(phi<_phiMax-0.0001)

  //Averaging quantities of interest
  BYAveVec/=wsum;

  //Last thing: building _BYAveLab
  //We have tried three options:
  //1) Taking the momentum of the B determined by the the YAve frame and setting the mass to mB. Give worst results than standard YAve
  // technique except for q2 where we've seen a very small improvement.
  // This also ends up with a mean BYAve.rho() of _pBCM/2 and even less...
  //2) The YAve frame is only used for the DIRECTON of the B momentum, we fix mB and keep _pBCM as its magnitude
  //   like this: Hep3Vector BYAveFinalVec = _pBCM * BYAveVec.unit();
  // This give worst results than 1).
  //3) Taking the 4-momentum of the B as in YAve, fixing mB, but conserving the boostVector as well; -> Basically the same as 1) 
  //
  //Overall, none of these options is as good as standard BYAve for Kin variables

  HepLorentzVector BYAve = HepLorentzVector(0,0,0,0);
  BYAve.setVectM(BYAveVec,_mB);

  XSLRotateAndBoost rAndB = XSLRotateAndBoost();
  HepLorentzVector BYAveLabo = rAndB.BoostFromFrame( BYAve, _Ups);
  
  return BYAveLabo;
}

void
YAveBooster::ComputeKin( HepLorentzVector VDLab )
{
  if(!_OK) { _q2=-666; _thL=-666; _thV=-666; _chi=-666; return; } 

  bool IsVectorDecay=true;
  if(VDLab==HepLorentzVector(0,0,0,0)) IsVectorDecay=false;

  //Initialize and setup
  XSLRotateAndBoost rAndB = XSLRotateAndBoost();
  _q2=0;
  _thL=0;
  _thV=0;
  _chi=0;
  double wsum=0.0;		// sum of weights 

  //Here we go!
  double phi=_phiMin;
  while(phi<_phiMax-0.0001) {

    //compute the B with this phi
    double pBx=_pBCM*_sinBY*cos(phi);
    double pBy=_pBCM*_sinBY*sin(phi);
    double pBz=_pBCM*_cosBY;
    Hep3Vector p3B=pBx*_xHatY+pBy*_yHatY+pBz*_zHatY;

    //weight and sum
    double cb=p3B.cosTheta(); //Note that we're in Ups(4S) coordinates! ;-)
    double sb2=1.0-cb*cb;
    double wt=sb2;  
    wsum+=wt;

    //Compute Stuff you're interested in:
    HepLorentzVector p4B = HepLorentzVector(0,0,0,0);
    p4B.setVectM(p3B,_mB);
    HepLorentzVector BLAB = rAndB.BoostFromFrame( p4B, _Ups);
    XSLKin* Kin;
    if(IsVectorDecay)
      {
	Kin  = new XSLKin(BLAB, _LepLab, _XuLab, VDLab);
	_thV += wt*Kin->theta_v();
	_chi += wt*Kin->chi();
      }
    else Kin = new XSLKin(BLAB, _LepLab, _XuLab);
    _q2 += wt*Kin->q2();
    _thL += wt*Kin->theta_l();
    delete Kin;

    //increment phi
    phi+=_dPhi;

  }//while(phi<_phiMax-0.0001)

  //Averaging quantities of interest
  _q2/=wsum;
  _thL/=wsum;
  _thV/=wsum;
  _chi/=wsum;
}

double 
YAveBooster::delThNu( HepLorentzVector PmissLab )
{
 if(!_OK) { return -666; } 

  //Initialize and setup
  XSLRotateAndBoost rAndB = XSLRotateAndBoost();
  double delThNu=0;
  double wsum=0.0;		// sum of weights 

  //Here we go!
  double phi=_phiMin;
  while(phi<_phiMax-0.0001) {

    //compute the B with this phi
    double pBx=_pBCM*_sinBY*cos(phi);
    double pBy=_pBCM*_sinBY*sin(phi);
    double pBz=_pBCM*_cosBY;
    Hep3Vector p3B=pBx*_xHatY+pBy*_yHatY+pBz*_zHatY;

    //weight and sum
    double cb=p3B.cosTheta(); //Note that we're in Ups(4S) coordinates! ;-)
    double sb2=1.0-cb*cb;
    double wt=sb2;  
    wsum+=wt;

    //Compute Stuff you're interested in:
    HepLorentzVector p4B = HepLorentzVector(0,0,0,0);
    p4B.setVectM(p3B,_mB);
    HepLorentzVector BLAB = rAndB.BoostFromFrame( p4B, _Ups);

    delThNu += wt*BLAB.vect().angle(PmissLab.vect());

    //increment phi
    phi+=_dPhi;

  }//while(phi<_phiMax-0.0001)

  //Averaging quantities of interest
  delThNu/=wsum;

  return delThNu;
}











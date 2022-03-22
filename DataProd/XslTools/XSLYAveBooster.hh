//*****************************************************************
//   
//   Creation: David Cote, Universite de Montreal, 02/02/04
//

#ifndef XSLYAVEBOOSTER
#define XSLYAVEBOOSTER

#include "CLHEP/Vector/LorentzVector.h"
class XSLKin;

class YAveBooster {

public:

  YAveBooster( YAveBooster & YAve );
  YAveBooster(HepLorentzVector LepLab, HepLorentzVector XuLab, int Bchrg );
  ~YAveBooster();

  //Output stuff:
  HepLorentzVector BYAveLab(); 
  void ComputeKin( HepLorentzVector VDLab = HepLorentzVector(0,0,0,0) );
  double q2() { return _q2; }
  double thL() { return _thL; }
  double thV() { return _thV; }
  double chi() { return _chi; }
  double cosBY() { return _cosBY; }
  
  double delThNu( HepLorentzVector PmissLab );

  //mainly for the copy constructor...
  int Bchrg() { return _Bchrg; }
  HepLorentzVector LepLab() { return _LepLab; }
  HepLorentzVector XuLab() { return _XuLab; }

private:
  void Init();

  //For the Copy constructor:
  int _Bchrg;
  HepLorentzVector _LepLab;
  HepLorentzVector _XuLab;

  //Internal utilities
  double _pBCM;
  double _mB;
  double _cosBY;
  double _sinBY;
  double _phiMin;
  double _phiMax;
  double _dPhi;
  Hep3Vector _xHatY;
  Hep3Vector _yHatY;
  Hep3Vector _zHatY;
  HepLorentzVector _Ups;
  
  //Output Stuff
  double _q2;
  double _thL;
  double _thV;
  double _chi;

  bool _OK;

};

#endif

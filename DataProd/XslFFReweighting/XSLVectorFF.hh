//*****************************************************************
//   
//   Creation: David Cote, Universite de Montreal, 12/11/03
//   see comments in the XSLEvtFFWeight.hh file or BAD#809

#ifndef XSLVECTORFF
#define XSLVECTORFF

#include "XslFFReweighting/XSLEvtFFWeight.hh"
#include <string>

class HepLorentzVector;
class XSLKin;

class XSLVectorFF : public XSLEvtFFWeight {
public:
  XSLVectorFF( double mB, double mXu, double q2, double theta_l, double theta_V, double chi );
  XSLVectorFF( HepLorentzVector BLab, HepLorentzVector LepLab, HepLorentzVector XuLab, HepLorentzVector XuDaughterLab);
  XSLVectorFF( XSLKin* DecayKin );
  virtual ~XSLVectorFF();

  //Event Weight from specific Generators
  virtual double FromISGW2ToThisModel();
  virtual double FromPHSPToThisModel();
  virtual double FromFLATQ2ToThisModel();
  virtual double FromSP4ToThisModel(){ return FromISGW2ToThisModel(); }
  virtual double FromSP5ToThisModel(){ return FromISGW2ToThisModel(); }
  virtual double FromSP6ToThisModel(){ return FromISGW2ToThisModel(); }
  virtual double FromSP7ToThisModel(){ return FromISGW2ToThisModel(); }

  //Differential Decay Rate
  virtual double dGammadQ2();
  double dGammadQ2(double q2);
  virtual double dGammadQ2dCtl();
  virtual double dGammadQ2dCtv();
  virtual double dGammadQ2dCtldCtv();
  virtual double dGammadQ2dCtldCtvdChi();
  virtual double Gamma();


  //Helicity Amplitudes
  double Hplus(){ return _Hplus; }
  double Hplus(double A1, double V, double kV);
  double Hplus2(){ return _Hplus*_Hplus; }
  double Hminus(){ return _Hminus; }
  double Hminus(double A1, double V, double kV);
  double Hminus2(){ return _Hminus*_Hminus; }
  double Hzero(){ return _Hzero; }
  double Hzero(double A1, double A2, double q2);
  double Hzero2(){ return _Hzero*_Hzero; }

  //Form Factors
  double A1(){return _A1;}
  double A2(){return _A2;}
  double V(){return _V;}


protected:
  void InitWithKin( XSLKin* MyDecayKin );
  virtual void Compute();

  void GetHelicityAmp();
  virtual void GetAllFF(double *A1, double *A2, double *V)=0;

  virtual void SetNormalizations(std::string mode)=0;

  virtual double GetA1(double q2)=0;
  virtual double GetA2(double q2)=0;
  virtual double GetV(double q2)=0;

  double AngularPartOfFullWeight();

  double _kV;
  double _mV;
  double _mV2;
  double _mB;
  double _mB2;
  double _y;
  double _q2;
  double _thL;
  double _thV;
  double _chi;
  double _FLATQ2Normalization;
  double _PHSPNormalization;
  double _ISGW2Normalization;

private:
  double _Hplus;
  double _Hminus;
  double _Hzero;
  double _A1;
  double _A2;
  double _V;

};

#endif






//*****************************************************************
//   
//   Creation: David Cote, Universite de Montreal, 12/11/03
//   see comments in the XSLEvtFFWeight.hh file or BAD #809

#ifndef XSLPSEUDOSCALARFF
#define XSLPSEUDOSCALARFF

#include "XslFFReweighting/XSLEvtFFWeight.hh"
#include <string>

class HepLorentzVector;
class XSLKin;

class XSLPseudoScalarFF : public XSLEvtFFWeight {
public:
  XSLPseudoScalarFF( double mB, double mXu, double q2, double theta_l );
  XSLPseudoScalarFF( HepLorentzVector BLab, HepLorentzVector LepLab, HepLorentzVector XuLab);
  XSLPseudoScalarFF( XSLKin* DecayKin );
  virtual ~XSLPseudoScalarFF();

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
  virtual double dGammadQ2dCtv(){ return 0; }        //note: could be implemented in principle, but would need thV and chi
  virtual double dGammadQ2dCtldCtv(){ return 0; }    //note: could be implemented in principle, but would need thV and chi
  virtual double dGammadQ2dCtldCtvdChi(){ return 0; }//note: could be implemented in principle, but would need thV and chi
  virtual double Gamma();

  //Form Factor
  double GetFplus(){ return GetFplus(_q2); }
  double Fplus(){ return _Fplus; }
  double Fplus2(){ return _Fplus*_Fplus; }


protected:
  void InitWithKin( XSLKin* MyDecayKin );
  virtual void Compute();
  virtual double GetFplus(double q2)=0;
  virtual void SetNormalizations(std::string mode)=0;

  double _mB;
  double _mXu;
  double _mXu2;
  double _pXu;
  double _q2;
  double _thL;
  double _FLATQ2Normalization;
  double _PHSPNormalization;
  double _ISGW2Normalization;

private:
  double _Fplus;
  
};

#endif







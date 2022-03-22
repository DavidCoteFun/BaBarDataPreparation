//*****************************************************************
//   
//   Creation: David Cote, 08/24/05, Universite de Montreal
//   see comments in the XSLEvtFFWeight.hh file or BAD #809
//
//   This class can produce reweighting with the functions of:
//     -R. Hill (hep-ph/0505129)        or
//     -Becirevic-Kaidalov (Phys. Lett. B478 (2000))

#ifndef XSLGENERICFPLUS
#define XSLGENERICFPLUS

#include "XslFFReweighting/XSLPseudoScalarFF.hh"

class XSLKin;

class XSLGenericFplus  : public XSLPseudoScalarFF {
public:
  XSLGenericFplus( double mB, double mXu, double q2, double theta_l,std::string mode,double alpha,double delta=0);
  XSLGenericFplus( HepLorentzVector BLab, HepLorentzVector LepLab, HepLorentzVector XuLab,std::string mode,double alpha,double delta=0);
  XSLGenericFplus( XSLKin* DecayKin,std::string mode,double alpha,double delta=0);
  virtual ~XSLGenericFplus(){};

protected:
  virtual double GetFplus(double q2);
  virtual void SetNormalizations(std::string mode);

  double _alpha;
  double _delta;
};

#endif




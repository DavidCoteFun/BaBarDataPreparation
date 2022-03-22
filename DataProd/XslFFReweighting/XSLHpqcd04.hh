//*****************************************************************
//   
//   Creation: Jochen Dingfelder 11/15/04 SLAC
//   see comments in the XSLEvtFFWeight.hh file or BAD #809
//
//   From the paper: hep-lat/0408019

#ifndef XSLHPQCD04
#define XSLHPQCD04

#include "XslFFReweighting/XSLPseudoScalarFF.hh"

class XSLKin;

class XSLHpqcd04  : public XSLPseudoScalarFF {
public:
  XSLHpqcd04( double mB, double mXu, double q2, double theta_l,std::string mode );
  XSLHpqcd04( HepLorentzVector BLab, HepLorentzVector LepLab, HepLorentzVector XuLab,std::string mode);
  XSLHpqcd04( XSLKin* DecayKin,std::string mode );
  virtual ~XSLHpqcd04(){};

protected:
  virtual double GetFplus(double q2);
  virtual void SetNormalizations(std::string mode);
  
};

#endif




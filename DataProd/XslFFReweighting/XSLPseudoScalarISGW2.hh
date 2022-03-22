//*****************************************************************
//   
//   Creation: David Cote, Universite de Montreal, 12/11/03
//   see comments in the XSLEvtFFWeight.hh file or BAD #809
//
//   Calculation from EvtGenModels/EvtISGW2.cc 

#ifndef XSLPSEUDOSCALARISGW2
#define XSLPSEUDOSCALARISGW2

#include "XslFFReweighting/XSLPseudoScalarFF.hh"

class XSLKin;

class XSLPseudoScalarISGW2  : public XSLPseudoScalarFF {
public:
  XSLPseudoScalarISGW2( double mB, double mXu, double q2, double theta_l, std::string mode="Ignore" );
  XSLPseudoScalarISGW2( HepLorentzVector BLab, HepLorentzVector LepLab, HepLorentzVector XuLab, std::string mode="Ignore" );
  XSLPseudoScalarISGW2( XSLKin* DecayKin, std::string mode="Ignore" );
  virtual ~XSLPseudoScalarISGW2(){};

protected:
  virtual double GetFplus(double q2);
  virtual void SetNormalizations(std::string mode);
  
};

#endif




//*****************************************************************
//   
//   Creation: David Cote, Universite de Montreal, 12/11/03
//   see comments in the XSLEvtFFWeight.hh file or BAD #809
//
//
//   A Quenched Fermilab calculation  -  hep-ph/0101023v1


#ifndef XSLFNAL01
#define XSLFNAL01

#include "XslFFReweighting/XSLPseudoScalarFF.hh"

class XSLKin;

class XSLFnal01  : public XSLPseudoScalarFF {
public:
  XSLFnal01( double mB, double mXu, double q2, double theta_l, std::string mode="Ignore" );
  XSLFnal01( HepLorentzVector BLab, HepLorentzVector LepLab, HepLorentzVector XuLab, std::string mode="Ignore" );
  XSLFnal01( XSLKin* DecayKin, std::string mode="Ignore" );
  virtual ~XSLFnal01(){};

protected:
  virtual double GetFplus(double q2);
  virtual void SetNormalizations(std::string mode);
  
};

#endif




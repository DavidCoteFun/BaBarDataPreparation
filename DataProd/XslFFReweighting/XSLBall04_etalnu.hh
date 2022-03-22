//*****************************************************************
//   
//   Creation: David Cote, Universite de Montreal, 12/11/03
//   Source: hep-ph/0110115 v1
//   see comments in the XSLEvtFFWeight.hh file or BAD #809
//


#ifndef XSLBALL04_ETALNU
#define XSLBALL04_ETALNU

#include "XslFFReweighting/XSLPseudoScalarFF.hh"


class XSLKin;

class XSLBall04_etalnu  : public XSLPseudoScalarFF {
public:
  XSLBall04_etalnu( double mB, double mXu, double q2, double theta_l, std::string mode="Ignore" );
  XSLBall04_etalnu( HepLorentzVector BLab, HepLorentzVector LepLab, HepLorentzVector XuLab, std::string mode="Ignore" );
  XSLBall04_etalnu( XSLKin* DecayKin, std::string mode="Ignore" );
  virtual ~XSLBall04_etalnu(){};

protected:
  virtual double GetFplus(double q2);
  virtual void SetNormalizations(std::string mode);
  
};

#endif




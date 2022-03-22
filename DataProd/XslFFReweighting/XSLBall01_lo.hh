//*****************************************************************
//   
//   Creation: David Cote, Universite de Montreal, 12/11/03
//   Source: hep-ph/0110115 v1
//   see comments in the XSLEvtFFWeight.hh file or BAD #809
//


#ifndef XSLBALL01_LO
#define XSLBALL01_LO

#include "XslFFReweighting/XSLPseudoScalarFF.hh"


class XSLKin;

class XSLBall01_lo  : public XSLPseudoScalarFF {
public:
  XSLBall01_lo( double mB, double mXu, double q2, double theta_l, std::string mode="Ignore" );
  XSLBall01_lo( HepLorentzVector BLab, HepLorentzVector LepLab, HepLorentzVector XuLab, std::string mode="Ignore" );
  XSLBall01_lo( XSLKin* DecayKin, std::string mode="Ignore" );
  virtual ~XSLBall01_lo(){};

protected:
  virtual double GetFplus(double q2);
  virtual void SetNormalizations(std::string mode);
  
};

#endif




//*****************************************************************
//   
//   Creation: David Cote, Universite de Montreal, 12/11/03
//   Source: hep-ph/0110115 v1
//   see comments in the XSLEvtFFWeight.hh file or BAD #809
//


#ifndef XSLBALL01_HI
#define XSLBALL01_HI

#include "XslFFReweighting/XSLPseudoScalarFF.hh"


class XSLKin;

class XSLBall01_hi  : public XSLPseudoScalarFF {
public:
  XSLBall01_hi( double mB, double mXu, double q2, double theta_l, std::string mode="Ignore" );
  XSLBall01_hi( HepLorentzVector BLab, HepLorentzVector LepLab, HepLorentzVector XuLab, std::string mode="Ignore" );
  XSLBall01_hi( XSLKin* DecayKin, std::string mode="Ignore" );
  virtual ~XSLBall01_hi(){};

protected:
  virtual double GetFplus(double q2);
  virtual void SetNormalizations(std::string mode);
  
};

#endif




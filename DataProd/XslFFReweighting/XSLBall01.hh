//*****************************************************************
//   
//   Creation: David Cote, Universite de Montreal, 12/11/03
//   Source: hep-ph/0110115 v1
//   see comments in the XSLEvtFFWeight.hh file or BAD #809
//


#ifndef XSLBALL01
#define XSLBALL01

#include "XslFFReweighting/XSLPseudoScalarFF.hh"


class XSLKin;

class XSLBall01  : public XSLPseudoScalarFF {
public:
  XSLBall01( double mB, double mXu, double q2, double theta_l, std::string mode="Ignore" );
  XSLBall01( HepLorentzVector BLab, HepLorentzVector LepLab, HepLorentzVector XuLab, std::string mode="Ignore" );
  XSLBall01( XSLKin* DecayKin, std::string mode="Ignore" );
  virtual ~XSLBall01(){};

protected:
  virtual double GetFplus(double q2);
  virtual void SetNormalizations(std::string mode);
  
};

#endif




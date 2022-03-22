//*****************************************************************
//   
//   Creation: Jochen Dingfelder 11/15/04 SLAC
//   see comments in the XSLEvtFFWeight.hh file or BAD #809
//
//   From the paper: hep-lat/0409116 v1

#ifndef XSLFNAL04
#define XSLFNAL04

#include "XslFFReweighting/XSLPseudoScalarFF.hh"

class XSLKin;

class XSLFnal04  : public XSLPseudoScalarFF {
public:
  XSLFnal04( double mB, double mXu, double q2, double theta_l,std::string mode );
  XSLFnal04( HepLorentzVector BLab, HepLorentzVector LepLab, HepLorentzVector XuLab,std::string mode);
  XSLFnal04( XSLKin* DecayKin,std::string mode );
  virtual ~XSLFnal04(){};

protected:
  virtual double GetFplus(double q2);
  virtual void SetNormalizations(std::string mode);
  
};

#endif




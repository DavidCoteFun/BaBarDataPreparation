//*****************************************************************
//   
//   Creation: David Cote, Universite de Montreal, 12/11/03
//   see comments in the XSLEvtFFWeight.hh file or BAD #809
//
//   From the LCSR calculation of P. Ball -  hep-ph/9805422


#ifndef XSLBALL98
#define XSLBALL98

#include "XslFFReweighting/XSLVectorFF.hh"

class XSLKin;

class XSLBall98  : public XSLVectorFF {
public:
  XSLBall98( double mB, double mXu, double q2, double theta_l, double theta_V, double chi, std::string mode="Ignore" );
  XSLBall98(HepLorentzVector BLab, HepLorentzVector LepLab, HepLorentzVector XuLab, HepLorentzVector XuDaughterLab, std::string mode="Ignore" );
  XSLBall98( XSLKin* DecayKin, std::string mode="Ignore" );
  virtual ~XSLBall98(){};

protected:
  virtual void GetAllFF(double *A1, double *A2, double *V);
  virtual void SetNormalizations(std::string mode);
  
  virtual double GetA1(double q2);
  virtual double GetA2(double q2);
  virtual double GetV(double q2);

};

#endif




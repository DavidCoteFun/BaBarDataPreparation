//*****************************************************************
//   
//   Creation: David Cote, Universite de Montreal, 12/11/03
//   see comments in the XSLEvtFFWeight.hh file or BAD #809
//
//   Calculation from EvtGenModels/EvtISGW2.cc 

#ifndef XSLVECTORISGW2
#define XSLVECTORISGW2

#include "XslFFReweighting/XSLVectorFF.hh"

class XSLKin;

class XSLVectorISGW2  : public XSLVectorFF {
public:
  XSLVectorISGW2( double mB, double mXu, double q2, double theta_l, double theta_V, double chi, std::string mode="Ignore" );
  XSLVectorISGW2(HepLorentzVector BLab, HepLorentzVector LepLab, HepLorentzVector XuLab, HepLorentzVector XuDaughterLab, std::string mode="Ignore" );
  XSLVectorISGW2( XSLKin* DecayKin, std::string mode="Ignore" );
  virtual ~XSLVectorISGW2(){};

protected:
  virtual void GetAllFF(double *A1, double *A2, double *V);
  virtual void SetNormalizations(std::string mode);

  virtual double GetA1(double q2);
  virtual double GetA2(double q2);
  virtual double GetV(double q2);

private:
  void GetAllISGW2FF(double *A1, double *A2, double *V, double q2);

  double EvtGetGammaji ( double z );
  double EvtGetas ( double massq, double massx );
  double EvtGetas ( double mass );
  
};

#endif




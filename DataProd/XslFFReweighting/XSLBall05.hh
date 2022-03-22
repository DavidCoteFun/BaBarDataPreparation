//*****************************************************************
//   
//   Creation: David Cote, Universite de Montreal, 12/11/03
//   see comments in the XSLEvtFFWeight.hh file or BAD #809
//
//   From the LCSR calculation of P. Ball -  to be published 
//                                           in late 04/ early 05


#ifndef XSLBALL05
#define XSLBALL05

#include "XslFFReweighting/XSLVectorFF.hh"

class XSLKin;

class XSLBall05  : public XSLVectorFF {
public:
  XSLBall05(double mB,double mXu,double q2,double theta_l,double theta_V,double chi,std::string mode,std::string ErrorScale="NoError");
  XSLBall05(HepLorentzVector BLab,HepLorentzVector LepLab,HepLorentzVector XuLab,HepLorentzVector XuDaughterLab,
	    std::string mode,std::string ErrorScale="NoError");
  XSLBall05(XSLKin* DecayKin,std::string mode,std::string ErrorScale="NoError");
  virtual ~XSLBall05(){};

protected:
  virtual void GetAllFF(double *A1, double *A2, double *V);
  virtual void SetNormalizations(std::string mode);
  void SetConstants(std::string mode);
  
  virtual double GetA1(double q2);
  virtual double GetA2(double q2);
  virtual double GetV(double q2);

private:
  double _r2_A1;
  double _mfit2_A1;
  double _r1_A2;
  double _r2_A2;
  double _mfit2_A2;
  double _r1_V;
  double _r2_V;
  double _mfit2_V;

  //function pointers set in the constructor from the ErrorScale argument
  double (XSLBall05::*_A1ScaleFunction)(double);
  double (XSLBall05::*_A2ScaleFunction)(double); 
  double (XSLBall05::*_VScaleFunction)(double); 
 
  //_ScaleFunction points to either of the following three options:
  double ReturnOne(double q2){ return 1.0; }
  double ErrorScaleFactorPos(double q2){ return (1.0 + ErrorMethodA(q2)); }
  double ErrorScaleFactorNeg(double q2){ return (1.0 - ErrorMethodA(q2)); }

  //utilities
  void Init(std::string mode, std::string ErrorScale);
  double ErrorMethodA(double q2);
  std::string _ErrorScale;
};

#endif




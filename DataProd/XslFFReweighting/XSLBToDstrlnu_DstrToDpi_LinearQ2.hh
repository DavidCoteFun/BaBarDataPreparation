//*****************************************************************
//   
//   Creation: David Cote, Universite de Montreal, 5 Nov 2004 
//   from Art Snyder & Amanda Weinstein's XslFFReweighting/DStarlnuPDF class
//   see comments in the XSLEvtFFWeight.hh file or BAD#809

#ifndef XSLBTODSTRLNU_DSTRTODPI_LINEARQ2
#define XSLBTODSTRLNU_DSTRTODPI_LINEARQ2

#include "XslFFReweighting/XSLEvtFFWeight.hh"

#include <math.h>

class HepLorentzVector;
class XSLKin;

class XSLBToDstrlnu_DstrToDpi_LinearQ2 : public XSLEvtFFWeight {
public:
  XSLBToDstrlnu_DstrToDpi_LinearQ2(double mB, double mDStar, double q2, double ctl, double ctv, double chi, double R1, double R2, double rho2);
  XSLBToDstrlnu_DstrToDpi_LinearQ2(HepLorentzVector BLab, HepLorentzVector LepLab, HepLorentzVector DStarLab, HepLorentzVector DStarDaugLab, 
			    double R1, double R2, double rho2);
  XSLBToDstrlnu_DstrToDpi_LinearQ2( XSLKin* DecayKin, double R1, double R2, double rho2 );
  virtual ~XSLBToDstrlnu_DstrToDpi_LinearQ2();


  //Event Weight from specific Generators
  virtual double FromISGW2ToThisModel();
  virtual double FromPHSPToThisModel();
  virtual double FromFLATQ2ToThisModel();
  virtual double FromSP4ToThisModel();
  virtual double FromSP5ToThisModel(){ return FromSP4ToThisModel(); }
  virtual double FromSP6ToThisModel(){ return FromSP4ToThisModel(); }
  virtual double FromSP7ToThisModel();

  //Differential Decay Rate
  virtual double dGammadQ2();
  virtual double dGammadQ2dCtl();
  virtual double dGammadQ2dCtv();
  virtual double dGammadQ2dCtldCtv();
  virtual double dGammadQ2dCtldCtvdChi();
  virtual double Gamma();

private:

  void InitWithKin( XSLKin* MyDecayKin );

  double norm();

  //access 
  const double& R1()const {return _R1;}
  const double& R2()const {return _R2;}
  const double& rho2()const {return _rho2;}
 

  //angular terms
  inline double gPP(double ctl,double ctd,double chi)const
  {return (1-ctd*ctd)*(1-ctl)*(1-ctl);}

  inline double gMM(double ctl,double ctd,double chi)const
  {return (1-ctd*ctd)*(1+ctl)*(1+ctl);}

  inline double g00(double ctl,double ctd,double chi)const
  {return 4*ctd*ctd*(1-ctl*ctl);}

  inline double gPM(double ctl,double ctd,double chi)const
  {return -2*(1-ctd*ctd)*(1-ctl*ctl)*cos(2*chi);}

  inline double gP0(double ctl,double ctd,double chi)const
  { return -4*sqrt(1-ctd*ctd)*ctd*sqrt(1-ctl*ctl)*(1-ctl)*cos(chi); }

  inline double gM0(double ctl,double ctd,double chi)const 
  { return  4*sqrt(1-ctd*ctd)*ctd*sqrt(1-ctl*ctl)*(1+ctl)*cos(chi); }

  double _R1;
  double _R2;
  double _rho2;

  double _mDStar;
  double _mB;
  double _q2;
  double _w;
  double _ctl;
  double _ctv;
  double _chi;

};

#endif






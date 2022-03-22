//===========================================================================
//  File and Version Information 
//  $Id: XSLFullLattice_Param.hh,v 1.3 2005/02/01 08:03:21 cote Exp $
//  Description:
//     Definition of class implementing form-factor reweighting for a 
//  fit of the Becirevic-Kaidalov parametrization to the FNAL01, JLQCD 2001,
//  and UKCQCD99 results, as performed by J. Dingfelder.  Derived from the 
//  base class XSLPseudoScalarFF written by David Cote-Ahern
//
// Author List:
//      Amanda Weinstein           Original Author
//
// Copyright Information:
//      Copyright (C) 2004         SLAC
//
// History:
// 
//===========================================================================

#ifndef XSLFULLLATTICE_PARAM_HH
#define XSLFULLLATTICE_PARAM_HH

#include "XslFFReweighting/XSLPseudoScalarFF.hh"

class XSLKin;

class XSLFullLattice_Param  : public XSLPseudoScalarFF {
public:
  XSLFullLattice_Param( double mB, double mXu, double q2, double theta_l, std::string mode="Ignore" );
  XSLFullLattice_Param( HepLorentzVector BLab, HepLorentzVector LepLab, HepLorentzVector XuLab, std::string mode="Ignore" );
  XSLFullLattice_Param( XSLKin* DecayKin, std::string mode="Ignore" );
  virtual ~XSLFullLattice_Param(){};

protected:
  virtual double GetFplus(double q2);
  virtual void SetNormalizations(std::string mode);
  
};

#endif // XSLFULLLATTICE_PARAM_HH




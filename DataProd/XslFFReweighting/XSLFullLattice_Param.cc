//===========================================================================
//  File and Version Information 
//  $Id: XSLFullLattice_Param.cc,v 1.6 2005/02/01 08:03:21 cote Exp $
//  Description:
//     Implementation of class implementing form-factor reweighting for a 
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
#include "BaBar/BaBar.hh"

#include "XslFFReweighting/XSLFullLattice_Param.hh"
#include "CLHEP/Vector/LorentzVector.h"
using std::cout;
using std::endl;
using std::string;

XSLFullLattice_Param::XSLFullLattice_Param( double mB, double mXu, double q2, double theta_l, string mode) : 
  XSLPseudoScalarFF(mB,mXu,q2,theta_l) 
{ Compute(); }

XSLFullLattice_Param::XSLFullLattice_Param( XSLKin* DecayKin, string mode) : XSLPseudoScalarFF(DecayKin) 
{ Compute(); }

XSLFullLattice_Param::XSLFullLattice_Param( HepLorentzVector BLab, HepLorentzVector LepLab, HepLorentzVector XuLab, string mode) : 
  XSLPseudoScalarFF( BLab, LepLab, XuLab) 
{ Compute(); }

double 
XSLFullLattice_Param::GetFplus(double q2)
{
  double w = 0;
  // alpha in particular has large errors
  double c_b = 0.44;
  double alpha = 0.19;
  double  mbstar = 28.4;
  double z = q2/(mbstar*mbstar);

  double num = c_b *(1-alpha);
  double den = (1-z)*(1-alpha*z);
  w = num/den;
  return w;
}


void
XSLFullLattice_Param::SetNormalizations(string mode)
{
  //Note: These normalization constants were determined empirically from MC samples of 1000k events (for most of the modes)
  //The quoted error is simply the statistical one, but there can be a bigger error coming from rare high weight events for
  //some region of the phase space, mainly for the Vector-ISGW2 generator.

  //if(mode=="Ignore") {
   _FLATQ2Normalization=1.0; _PHSPNormalization=1.0; _ISGW2Normalization=1.0; 

   return;
}

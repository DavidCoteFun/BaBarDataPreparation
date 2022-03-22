//*****************************************************************
//
//   This class makes the exact same RotateAndBoost than the standard BetaCoreTools/BtaBooster
//   It's intended to be used in software that doesn't want/need to deal with BtaCandidates.
//
//   In short, the BoostToFrame method differs from the standard HepLorentzVector::boost 
//   by returning the boosted vector in a coordinates system where z is the boost direction.
//   This is sometimes what you really want if you're interested in absolute values of coordinates
//   such as the polar angle in the Ups(4S) frame.
//
//   See: http://www.slac.stanford.edu/BFROOT/www/Physics/Tools/BetaTools/btaBooster.html
//   
//   Creation: David Cote, Universite de Montreal, 02/10/04
//
#include "BaBar/BaBar.hh"


#include "XslTools/XSLRotateAndBoost.hh"

#include "CLHEP/Vector/LorentzRotation.h"
//#include "CLHEP/Vector/Rotation.h"       \__ These two guys are needed by this class but already included by LorentzRotation
//#include "CLHEP/Vector/LorentzVector.h"  /

XSLRotateAndBoost::XSLRotateAndBoost(){}

XSLRotateAndBoost::~XSLRotateAndBoost(){}

HepLorentzVector
XSLRotateAndBoost::BoostToFrame( HepLorentzVector VecToBoostInOriginalFrame, HepLorentzVector NewFrame )
{
  //Like BtaBooster::BoostTo
  HepLorentzRotation rAndB = GetBoostMatrix( NewFrame );
  rAndB.invert();

  HepLorentzVector BoostedVec = rAndB * VecToBoostInOriginalFrame;
  return BoostedVec;
 
}

HepLorentzVector
XSLRotateAndBoost::BoostFromFrame( HepLorentzVector VecToBoostInNewFrame, HepLorentzVector NewFrame )
{
  //Like BtaBooster::BoostFrom
  HepLorentzRotation rAndB = GetBoostMatrix( NewFrame );

  HepLorentzVector BoostedVec = rAndB * VecToBoostInNewFrame;
  return BoostedVec;
}


HepLorentzRotation
XSLRotateAndBoost::GetBoostMatrix( HepLorentzVector Frame )
{
  //-----------------------------
  //RotateAndBoost
  // the boost vector
  Hep3Vector boostVector( Frame.boostVector() );
  //
  // rotation matrix and boost
  double boost  = boostVector.mag();  
  Hep3Vector boostAlongZ( 0., 0., boost );
  HepLorentzRotation boostPart( boostAlongZ );  
  
  double alpha( boostVector.phi()   );
  double beta(  boostVector.theta() );  
  double gamma(-boostVector.phi()   );  

  HepRotation euler;  
  euler.rotateZ( gamma );  
  euler.rotateY( beta  );
  euler.rotateZ( alpha );
  
  HepLorentzRotation rotationPart( euler );
  HepLorentzRotation product( rotationPart*boostPart );

  return product;
}


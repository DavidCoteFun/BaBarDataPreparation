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

#ifndef XSLYROTATEANDBOOST
#define XSLYROTATEANDBOOST

class HepLorentzRotation;
class HepLorentzVector;

class XSLRotateAndBoost {

public:

  XSLRotateAndBoost();
  ~XSLRotateAndBoost();

  HepLorentzVector BoostToFrame( HepLorentzVector VecToBoostInOriginalFrame, HepLorentzVector NewFrame );
  HepLorentzVector BoostFromFrame( HepLorentzVector VecToBoostInNewFrame, HepLorentzVector NewFrame );
  HepLorentzRotation GetBoostMatrix( HepLorentzVector Frame );

};

#endif

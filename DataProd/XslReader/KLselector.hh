//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: KLselector.hh,v 1.3 2004/10/28 06:53:54 cote Exp $
//
// Description:
//	Class KLselector
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      PAPPAGALLO MARCO           	
//
// INPUT: Flag = 1 ---> Neutral Candidate from CalorClusterNeutral List (EMC)
//        Flag = 2 ---> Neutral Candidate from NeutralHad List (IFR)
//
//
// OUTPUT:Likelihood[0]={1.1 EMC candidate with energy < 0.2 GeV
//                       1.2 EMC candidate with a matched IFR neutral object
//                       1.3 EMC candidate without a matched IFR neutral object
//                       2.1 IFR candidate with hit layers < 2
//                       2.2 IFR candidate with a matched EMC neutral object
//                       2.3 IFR candidate without a matched EMC neutral object}
//
//        Likelihood[1]=likelihood value
//
//
// NOTE: Use BetaCoreTools V00-01-22 in case of incompatibility with BtaCandMap.hh
//
//------------------------------------------------------------------------

#ifndef KLSELECTOR_HH
#define KLSELECTOR_HH

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "Framework/AppModule.hh"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
class BtaCandidate;

//		---------------------
// 		-- Class Interface --
//		---------------------
 
class KLselector {

//--------------------
// Instance Members --
//--------------------

public:

  // Constructors
  KLselector(){};

  // Destructor
  // Operations

   static void Selector(AbsEvent* , int,  const BtaCandidate*, double [], int);

protected:

private:

   static double like2000_emc(const BtaCandidate*);
   static double like2001_emc(const BtaCandidate*);
   static double like2002_emc(const BtaCandidate*);
   static double like2000_ifr(const BtaCandidate*);
   static double like2001_ifr(const BtaCandidate*);
   static double like2002_ifr(const BtaCandidate*);

};

#endif








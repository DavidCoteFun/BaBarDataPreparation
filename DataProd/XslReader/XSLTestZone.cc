//------------------------------------------------------------------------------//
//                                                                              //
//  Analysis module for the B -> pi/pi0/eta/eta'/rhoC/rho0/omega l nu analyses  //
//                                                                              //
//     Sylvie Brunet  2002-2004 Universite de Montreal                          //
//     David Cote     2003-2004 Universite de Montreal                          //
//     Benoit Viaud        2004 Universite de Montreal                          //
//                                                                              //
//------------------------------------------------------------------------------//
#include "BaBar/BaBar.hh"

#include "XslReader/XSLReader.hh"

#include "HepTuple/Tuple.h"
#include "AbsEvent/getTmpAList.hh"
#include "PDT/Pdt.hh"
#include "CLHEP/Alist/AIterator.h"

///////////////

void
XSLReader::TestZoneOne( AbsEvent* anEvent )
{


  BtaCandidate* TestCand;

  HepAList<BtaCandidate>* EtapRhoGam;
  getTmpAList (anEvent, EtapRhoGam,  IfdStrKey(HepString("etaPrgDefault")) );

  HTValOrderedVector<float> mRG,pRG;
  HepAListIterator<BtaCandidate> iter_EtapRhoGam(*EtapRhoGam);  
  while ( 0 != ( TestCand = iter_EtapRhoGam()) ) 
    {
      mRG.append(TestCand->mass());
      pRG.append(TestCand->p());
    }

  _TestZoneOneNtuple->column("NbRG",EtapRhoGam->length(),0 ,"indRG",HTRange<int>(0,1000));
  _TestZoneOneNtuple->column("mRG",mRG,"NbRG",0, "indRG");
  _TestZoneOneNtuple->column("pRG",pRG,"NbRG",0, "indRG");



  HepAList<BtaCandidate>* EtapEtaPiPi;
  getTmpAList (anEvent, EtapEtaPiPi,  IfdStrKey(HepString("etaPeppDefault")) );

  HTValOrderedVector<float> mEPP,pEPP;
  //HTValOrderedVector<int> overlap;  
  //BtaCandidate* eta=0,*fille;
  //int UIDplus=-666,UIDminus=-666,same=0;
  HepAListIterator<BtaCandidate> iter_EtapEtaPiPi(*EtapEtaPiPi);  
  while ( 0 != ( TestCand = iter_EtapEtaPiPi()) ) 
    {
      mEPP.append(TestCand->mass());
      pEPP.append(TestCand->p());

      //Get the info...
      //int lund = TestCand->pdtEntry()->lundId();
      //if(lund==PdtLund::eta) eta = new BtaCandidate(*TestCand);
      //else if(lund==PdtLund::pi_plus) UIDplus = (int)TestCand->uid();
      //else if(lund==PdtLund::pi_minus) UIDminus = (int)TestCand->uid();
    }

  //Something's wrong with that crap!
  /*
  //...then use it!  
  iter_EtapEtaPiPi.rewind();
  while ( 0 != ( TestCand = iter_EtapEtaPiPi()) ) 
    {
      int n = eta->nDaughters();
      if(eta!=0 && n!=0)
	{
	  HepAListIterator<BtaCandidate>iter_fille = eta->daughterIterator();
	  while ( 0 != (fille =iter_fille()) )   
	    {
	      int lund2 = fille->pdtEntry()->lundId();
	      if(lund2==PdtLund::pi_plus && UIDplus==(int)fille->uid()) same=1;
	      else if(lund2==PdtLund::pi_minus && UIDminus==(int)fille->uid()) same=1;	  
	    }
	}
      overlap.append(same);
    }
  
  delete eta;
  */

  _TestZoneOneNtuple->column("NbEPP",EtapEtaPiPi->length(),0 ,"indEPP",HTRange<int>(0,1000));
  _TestZoneOneNtuple->column("mEPP",mEPP,"NbEPP",0, "indEPP");
  _TestZoneOneNtuple->column("pEPP",pEPP,"NbEPP",0, "indEPP");
  //_TestZoneOneNtuple->column("overlap",overlap,"NbEPP",0, "indEPP");
  
  _TestZoneOneNtuple->dumpData();

  return;
}

#include "EidData/EidEventTriplet.hh"
#include "EidData/EidCondKeyTriplet.hh"
#include "BdbTime/BdbTime.hh"
#include "AbsEvent/AbsEventID.hh"
using std::cout;
using std::endl;

void
XSLReader::TestZoneTwo( AbsEvent* anEvent )
{

  if ( _lookMCTruth.value()==true ) 
    {

      HepAList< AbsEventID >* eidList( 0 );
      getTmpAList( anEvent , eidList , IfdStrKey("AbsEventIDList") );
      BdbTime bdbFGTime;
      BdbTime bdbBGTime;
      EidCondKeyTriplet bkgroundTriplet;
      EidCondKeyTriplet foregroundTriplet;
      
      foregroundTriplet = ((*eidList)[ 0 ])->condKeyTriplet();
      bdbFGTime = foregroundTriplet.key();
      std::string fString = bdbFGTime.asString("%c", BdbTime::Local);
      cout << "Foreground conditions time is " << fString << endl;
      
      bkgroundTriplet = ((*eidList)[ 1 ])->condKeyTriplet();
      bdbBGTime = bkgroundTriplet.key();
      std::string bString = bdbBGTime.asString("%c", BdbTime::Local);
      cout << "Background conditions time is " << bString << endl;
    }

  AbsEventID* eventID = Ifd<AbsEventID>::get( anEvent , "AbsEventID" );
  if (eventID != 0) 
    {
      BdbTime myTime = eventID->condKeyTriplet().key();
      std::string myString = myTime.asString("%c", BdbTime::Local);
      cout << "my conditions time is " << myString << endl;
    }
  
  return;
}











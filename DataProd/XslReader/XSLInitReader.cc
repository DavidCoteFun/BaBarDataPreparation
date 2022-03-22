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

#include "AbsEnv/AbsEnv.hh"
#include "AbsEvent/AbsEvent.hh"
#include "AbsEvent/getTmpAList.hh"
#include "ErrLogger/ErrLog.hh"
#include "HepTuple/Tuple.h"
#include "HepTuple/TupleManager.h"
#include "GenEnv/GenEnv.hh"
#include "UsrData/UsrEventBlock.hh"
#include "PDT/Pdt.hh"
#include "BetaMC/BtaMcAssocGHit.hh"
using std::cout;
using std::endl;

///////////////

XSLReader::XSLReader( const char* const theName, 
					  const char* const theDescription )
  : AppModule( theName, theDescription )
  //StrKey 
  , _eventInfoList("eventInfoList", this, "Default")
  , _InputYList("inputYList",this,"XSLBtoXulnuSkimmedYlist") 
  , _InputTrkList("InputTrkList",this,"GoodTracksVeryLoose")
  , _InputPhotonList("InputPhotonList",this,"CalorNeutral")
  , _EmcDigiList("EmcDigiList",this,"Default")
  , _InputIfrList("InputIfrList",this,"NeutralHadrons")
  , _TruthMap("TruthMap", this, "GHit") 
  , _TruthList("TruthList", this, "MCTruth") 
  , _eventData("eventData", this, "XSLBtoXulnuEventData") 
  , _signalModeString("signalModeString",this,"pilnu")
  , _readingMode("readingMode",this,"Run4") //BypassUsrEventData,Run3,Run4

  //double (cuts)
  , _r2allMAX("r2allMAX",this,1)//no cut:>=1
  , _ishMIN("ishMIN",this,0)//no cut:0, cut:1
  , _ptMaxGhost("ptMaxGhost",this,0.35) //no cut: 100 
  , _nDchMinGhost("nDchMinGhost",this,1) //no cut: -1
  , _nSvtMinLooper("nSvtMinLooper",this,1) //no cut: -1
  , _cosThetaMaxLooper("cosThetaMaxLooper",this,0.2) //no cut: 1000
  , _ptMaxLooper("ptMaxLooper",this,0.25) //no cut: 1000
  , _dpMaxToMatch("dpMaxToMatch",this,0.2) //no cut: 1000

  //bool
  , _analyzePilnu("AnalyzePilnu",this,true)
  , _analyzePi0lnu("AnalyzePi0lnu",this,true)
  , _analyzeEtalnu("AnalyzeEtalnu",this,true)
  , _analyzeEtaplnu("AnalyzeEtaplnu",this,true)
  , _analyzeRhoClnu("AnalyzeRhoClnu",this,true)
  , _analyzeRho0lnu("AnalyzeRho0lnu",this,true)
  , _analyzeOmegalnu("AnalyzeOmegalnu",this,true)
  , _analyzeGammalnu("AnalyzeGammalnu",this,true)
  , _MyVerbose("MyVerbose", this, false)
  , _lookMCTruth("lookMCTruth", this, false)
  , _dumpEventInfo("dumpEventInfo",this,true)
  , _dumpTrkInfo("dumpTrkInfo",this,true)
  , _dumpPhotonInfo("dumpPhotonInfo",this,true)
  , _dumpPi0Info("dumpPi0Info",this,true)
  , _dumpCoupleInfo("dumpCoupleInfo",this,true)
  , _sigOnly("sigOnly",this,true)
  , _goToTestZone("goToTestZone",this,false)
  , _applyCuts("applyCuts",this,false)

  //UsrData
  , _decayMode("decayMode")
{
  //StrKey
  commands()->append( & _InputYList );
  commands()->append( & _TruthMap );
  commands()->append( & _TruthList );
  commands()->append( & _eventData );
  commands()->append( & _InputTrkList );
  commands()->append( & _InputPhotonList );
  commands()->append( & _EmcDigiList );
  commands()->append( & _InputIfrList );
  commands()->append( & _signalModeString);
  commands()->append( & _readingMode );

  //double
  commands()->append( & _r2allMAX );
  commands()->append( & _ishMIN );
  commands()->append( & _nDchMinGhost );
  commands()->append( & _ptMaxGhost );
  commands()->append( & _nSvtMinLooper );
  commands()->append( & _cosThetaMaxLooper );
  commands()->append( & _ptMaxLooper );
  commands()->append( & _dpMaxToMatch );

  //bool
  commands()->append( & _analyzePilnu );
  commands()->append( & _analyzePi0lnu );
  commands()->append( & _analyzeEtalnu );
  commands()->append( & _analyzeEtaplnu );
  commands()->append( & _analyzeRhoClnu );
  commands()->append( & _analyzeRho0lnu );
  commands()->append( & _analyzeOmegalnu );
  commands()->append( & _analyzeGammalnu );
  commands()->append( & _MyVerbose );
  commands()->append( & _lookMCTruth );
  commands()->append( & _dumpEventInfo );
  commands()->append( & _dumpTrkInfo );
  commands()->append( & _dumpPhotonInfo );
  commands()->append( & _dumpPi0Info );
  commands()->append( & _dumpCoupleInfo );
  commands()->append( & _sigOnly );
  commands()->append( & _goToTestZone );
  commands()->append( & _applyCuts );
}

XSLReader::~XSLReader() {}


AppResult
XSLReader::beginJob( AbsEvent* anEvent )
{
  ErrMsg(routine)<<"begin reader Job"<<endmsg; 

  HepTupleManager* manager = gblEnv->getGen()->ntupleManager();
  assert(manager != 0);

  //ntuples
  _MainNtuple = manager->ntuple("MainNtuple",1);
  assert( NULL != _MainNtuple);
  _JobNtuple = manager->ntuple("JobNtuple",2);
  assert( NULL != _JobNtuple);
  _TestZoneOneNtuple = manager->ntuple("TestZoneOneNtuple",10);
  assert( NULL != _TestZoneOneNtuple);
  _piOnlyNtuple = manager->ntuple("piOnlyNtuple",101);
  assert( NULL != _piOnlyNtuple);
  _pi0OnlyNtuple = manager->ntuple("pi0OnlyNtuple",102);
  assert( NULL != _pi0OnlyNtuple);
  _etaOnlyNtuple = manager->ntuple("etaOnlyNtuple",103);
  assert( NULL != _etaOnlyNtuple);
  _rho0OnlyNtuple = manager->ntuple("rho0OnlyNtuple",104);
  assert( NULL != _rho0OnlyNtuple);
  _rhoCOnlyNtuple = manager->ntuple("rhoCOnlyNtuple",105);
  assert( NULL != _rhoCOnlyNtuple);
  _omegaOnlyNtuple = manager->ntuple("omegaOnlyNtuple",106);
  assert( NULL != _omegaOnlyNtuple);


  _sigModeString=_signalModeString.value();
  if(_sigModeString=="pilnu") { _sigModeInt=1; }
  else if(_sigModeString=="pi0lnu") { _sigModeInt=2; }
  else if(_sigModeString=="etalnu") { _sigModeInt=3; }
  else if(_sigModeString=="etaplnu") { _sigModeInt=4; }
  else if(_sigModeString=="rhoClnu") { _sigModeInt=5; }
  else if(_sigModeString=="rho0lnu") { _sigModeInt=6; }
  else if(_sigModeString=="omegalnu") { _sigModeInt=7; }
  else _sigModeInt=0;


  _nbEvtProcessed=0;

  return AppResult::OK;
}


bool
XSLReader::InitializeEvent( AbsEvent* anEvent )
{
  
  ///////First get the had-l (Y) list...
  getTmpAList(anEvent,_YList,_InputYList.value());
  if (_YList == 0)   ErrMsg(fatal)<<"No YList!!"<<endmsg;
  
  int len=_YList->length();
  if( _MyVerbose.value())  cout<<"len: "<<len<<endl;
  if(len==0 && _applyCuts.value() ) return false; 

  _isPi = _analyzePilnu.value();
  _isPi0 = _analyzePi0lnu.value();
  _isEta = _analyzeEtalnu.value();
  _isEtap = _analyzeEtaplnu.value();
  _isRhoC = _analyzeRhoClnu.value();
  _isRho0 = _analyzeRho0lnu.value();
  _isOmega = _analyzeOmegalnu.value();
  _isGamma = _analyzeGammalnu.value(); 

  //Then look at the (UsrEventData) tag bit info...
  if(_readingMode.value()!="BypassUsrEventData")
    {
      //Event cut: reading back the bit mask...
     UsrEventBlock * eventData = UsrIfd<UsrEventBlock>::get( anEvent, _eventData.value() );
     
     bool found = false;
     if ( eventData != 0 ) found = eventData->get( _decayMode );    
     assert(found);
     int i=_decayMode;
          
     if(_readingMode.value()=="Run3")
       {
	 int a=i&1;
	 if( _isPi && (a!=1) ) _isPi= false;
	 a=i&2;
	 if( _isPi0 && (a!=2) ) _isPi0= false;
	 a=i&4;
	 if( _isEta && (a!=4) ) _isEta= false;
	 a=i&8;
	 if( _isRhoC && (a!=8) ) _isRhoC= false;
	 a=i&16;
	 if( _isRho0 && (a!=16) ) _isRho0= false;
	 a=i&32;
	 if( _isOmega && (a!=32) ) _isOmega= false;
	 a=i&64;
	 if( _isGamma && (a!=64) ) _isGamma= false;
       }
     if(_readingMode.value()=="Run4")
       {
	 int a=i&1;
	 if( _isPi && (a!=1) ) _isPi= false;
	 a=i&2;
	 if( _isPi0 && (a!=2) ) _isPi0= false;
	 a=i&4;
	 if( _isEta && (a!=4) ) _isEta= false;
	 a=i&8;
	 if( _isEtap && (a!=8) ) _isEtap= false;
	 a=i&16;
	 if( _isRhoC && (a!=16) ) _isRhoC= false;
	 a=i&32;
	 if( _isRho0 && (a!=32) ) _isRho0= false;
	 a=i&64;
	 if( _isOmega && (a!=64) ) _isOmega= false;
	 a=i&128;
	 if( _isGamma && (a!=128) ) _isGamma= false;	 
       }

     if( _MyVerbose.value())
       {
	 cout<<"found "<<found<<endl;
	 cout<<"decayMode"<<i<<endl;
       }
     if(!(_isPi||_isPi0||_isEta||_isEtap||_isRhoC||_isRho0||_isOmega||_isGamma)) { return false; } //This event won't be dropped in the ntuple!
    }  

  //Get event's basic lists...
  getTmpAList(anEvent,BasicTrkList,_InputTrkList.value());
  getTmpAList(anEvent,BasicBumpList,_InputPhotonList.value());
  getTmpAList(anEvent,BasicIfrList,_InputIfrList.value());

  ////////////Get MC Truth map...
  truthMapWithCon=0;
  if ( _lookMCTruth.value()==true )
    {
      //The main truth list
      getTmpAList (anEvent, _MCTruthList, _TruthList.value());

      //Emiss composition stuff...
      _nBadRecoed=0;
      _BadRecoed=HepLorentzVector(0,0,0,0);

      _dpMax=_dpMaxToMatch.value();

      getTmpAList( anEvent,_MatchedTruthList,IfdStrKey(HepString("MatchedTruthList")) );

      //MC-Reco map
      truthMapWithCon = (BtaMcAssocGHit*)Ifd< BtaMcAssoc >::get( anEvent, "GHit" );
      if (!truthMapWithCon) ErrMsg(fatal) << "No GHit matcher found" << endmsg;

    }//if lookMC


  return true;
}



AppResult
XSLReader::endJob( AbsEvent* anEvent )
{

 _JobNtuple->column("NbProcessedEvt", _nbEvtProcessed );
 _JobNtuple->dumpData();

  return AppResult::OK;
}

//////////////






//------------------------------------------------------------------------------//
//                                                                              //
//  Analysis module for the B -> pi/pi0/eta/eta'/rhoC/rho0/omega l nu analyses  //
//                                                                              //
//     Sylvie Brunet  2002-2004 Universite de Montreal                          //
//     David Cote     2003-2004 Universite de Montreal                          //
//     Benoit Viaud        2004 Universite de Montreal                          //
//                                                                              //
//------------------------------------------------------------------------------//


#ifndef XSLREADER_HH
#define XSLREADER_HH

//----------------------
// Base Class Headers --
//----------------------
#include "Framework/AppModule.hh"


//------------------------------------
// Collaborating Class Declarations --
//------------------------------------

class HepHistogram;
class HepTuple;
class AbsEvent;
class BtaMcAssocGHit;
class EventInfo;
class XSLMCTruthAnalyzer;
class XSLRecoYAnalyzer;

#include "AbsParm/AbsParmIfdStrKey.hh"
#include "Framework/AbsParmBool.hh"
#include "Framework/AbsParmDouble.hh"
#include "Framework/AbsParmString.hh"
#include "Beta/BtaCandidate.hh"
#include "CLHEP/Alist/AList.h"
#include "CLHEP/Alist/AListBase.h"
#include "HepTuple/HTValOrderedVector.h"
#include "UsrData/UsrVariable.hh"

//		---------------------
// 		-- Class Interface --
//		---------------------
 

class XSLReader : public AppModule {

//--------------------
// Instance Members --
//--------------------

public:

  // Constructors
  XSLReader( const char* const theName, const char* const theDescription );
  
  // Destructor
  virtual ~XSLReader( );
  
  // Operations
  virtual AppResult           beginJob( AbsEvent* anEvent );
  virtual AppResult           event   ( AbsEvent* anEvent );
  virtual AppResult           endJob  ( AbsEvent* anEvent );
  

  bool acceptedEvent();
  bool GetEventInfoFromReco( AbsEvent* anEvent );
  bool InitializeEvent( AbsEvent* anEvent );
  void TestZoneOne( AbsEvent* anEvent );
  void TestZoneTwo( AbsEvent* anEvent );
  void FillMainNtuple( AbsEvent* anEvent );
  void FillNtupleWithEventVariables( HepTuple* ); 
  void FillNtupleWithBumpListBlock( AbsEvent* , HepTuple* );
  BtaCandidate* GetNearestTrk( AbsEvent* anEvent, BtaCandidate* c, double &distToClosestTrk );
  int NearestTrkIndBasedOnAlpha(BtaCandidate* neut);
  void FillNtupleWithIfrListBlock( AbsEvent* anEvent, HepTuple* ntuple );
  void FillNtupleWithTrkListBlock( AbsEvent* anEvent, HepTuple* ntuple );
  XSLMCTruthAnalyzer* FillNtupleWithMCTruthEventVariables( AbsEvent* anEvent, HepTuple* ); 
  void FillNtupleWithMCTruthBlockInfo(AbsEvent* anEvent, HepTuple* ntuple);
  void FillNtupleWithCoupleVariables( string mode,  HepAList<BtaCandidate>* List, HepTuple* ntuple, AbsEvent* anEvent, XSLMCTruthAnalyzer* EvtTruth);
  void Unused();

  int GetBumpIndex( BtaCandidate candid );
  int GetTrkIndex( BtaCandidate candid );
  int GetBremPhotonIndex( BtaCandidate lepton );
  bool CandInList( BtaCandidate* cand, HepAList<BtaCandidate>* list );
  bool CandIsDaughterInList( BtaCandidate* cand, HepAList<BtaCandidate>* list );
  int CandIndexInList( BtaCandidate* cand, HepAList<BtaCandidate>* list );
  void GetPidBit( AbsEvent* anEvent, BtaCandidate* c, int *bita, int *bitb );
  float EncodeValAndErrInFloat(float val, float err);

  float IsSignal( XSLRecoYAnalyzer* recoY, XSLMCTruthAnalyzer* EvtTruth ); 
  int IsSignalOld( BtaCandidate* recoY ); 
  BtaCandidate* TruthMatch( BtaCandidate* rec, float &Consistency, float &OpeningAngle, float &delP,
			    int &nMatch, bool anyLund=true, bool noQualityCut=false );
  bool TruthMatchTheseCand( BtaCandidate* rec, BtaCandidate* mc, float &Consistency, 
			    bool mcIsMother=false, bool anyLund=true, bool noQualityCut=false );
  BtaCandidate* TruthMatchMerged( BtaCandidate* rec, float &OpeningAngle, float &delP, int &nFound, bool anyLund=true, bool noQualityCut=false );
  bool TruthMatchTheseMergedCand( BtaCandidate* rec, BtaCandidate* mc, bool mcIsMother=false, bool anyLund=true, bool noQualityCut=false );
  bool IndicesOfTruthMatched( BtaCandidate* rec, int ind[4], float cons[4], int &nMatch);
  float DeltaEBAD53(HepLorentzVector B);
  float mESBAD53(Hep3Vector B);

private:
  AbsParmIfdStrKey _eventInfoList;
  AbsParmIfdStrKey _InputYList;
  AbsParmIfdStrKey _InputTrkList;
  AbsParmIfdStrKey _InputPhotonList;
  AbsParmIfdStrKey _EmcDigiList;
  AbsParmIfdStrKey _InputIfrList;
  AbsParmIfdStrKey _TruthMap;    
  AbsParmIfdStrKey _TruthList;    
  AbsParmIfdStrKey _eventData;
  AbsParmString _signalModeString;
  AbsParmString _readingMode;

  AbsParmDouble      _r2allMAX;
  AbsParmDouble        _ishMIN;
  AbsParmDouble    _ptMaxGhost;
  AbsParmDouble  _nDchMinGhost;
  AbsParmDouble _nSvtMinLooper;
  AbsParmDouble _cosThetaMaxLooper;
  AbsParmDouble _ptMaxLooper;
  AbsParmDouble _dpMaxToMatch;

  AbsParmBool _analyzePilnu;
  AbsParmBool _analyzePi0lnu;
  AbsParmBool _analyzeEtalnu;
  AbsParmBool _analyzeEtaplnu;
  AbsParmBool _analyzeRhoClnu;
  AbsParmBool _analyzeRho0lnu;
  AbsParmBool _analyzeOmegalnu;
  AbsParmBool _analyzeGammalnu;
  AbsParmBool _MyVerbose;
  AbsParmBool _lookMCTruth;
  AbsParmBool _dumpEventInfo;
  AbsParmBool _dumpTrkInfo;
  AbsParmBool _dumpPhotonInfo;
  AbsParmBool _dumpPi0Info;
  AbsParmBool _dumpCoupleInfo;
  AbsParmBool _sigOnly;
  AbsParmBool _goToTestZone;
  AbsParmBool _applyCuts;

  HepAList<BtaCandidate>* _YList;
  HepAList<BtaCandidate>* BasicTrkList;
  HepAList<BtaCandidate>* BasicBumpList;
  HepAList<BtaCandidate>* BasicIfrList;
  HepAList<BtaCandidate>* _MCTruthList;
  HepAList<BtaCandidate>* _MatchedTruthList;

  HepTuple*        _MainNtuple;
  HepTuple*        _JobNtuple;
  HepTuple*        _TestZoneOneNtuple;
  HepTuple*        _piOnlyNtuple;
  HepTuple*        _pi0OnlyNtuple;
  HepTuple*        _etaOnlyNtuple;
  HepTuple*        _rho0OnlyNtuple;
  HepTuple*        _rhoCOnlyNtuple;
  HepTuple*        _omegaOnlyNtuple;

  UsrVariable<int> _decayMode;  

  BtaMcAssocGHit* truthMapWithCon;

  bool _isPi;
  bool _isPi0;
  bool _isEta;
  bool _isEtap;
  bool _isRhoC;
  bool _isRho0;
  bool _isOmega;
  bool _isGamma;
  bool _trueSigFound;

  int _nbEvtProcessed;

  //From CoupleStuff
  int _nDaughters;

  //From MCTruth
  int _nBadRecoed;
  HepLorentzVector _BadRecoed;
  double _dpMax;
  int _sigModeInt;
  string _sigModeString;

  //From EventStuff
  int _runNb;
  int _modulus;
  long int _date; //year/month/day as a long int
  int _platform;
  int _partition;
  int _upperID;
  int _lowerID;
  int _L3OutEmc;
  int _L3OutDch;
  EventInfo* _eventInfo;
  float _eTotCal;
  float _s;
  double _EbeamCM;
  float _sphericityall;	//Sphericity of the Event CMS
  float _r2;			//2nd Fox-Wolfram Moment with some tracks/neutrals selection
  float _r2all;		//2nd Fox-Wolfram Moment with everything
  int _ismultihadint;		//# Charged Tracks > 2 et R2 < 0.98
  int _isBCMH;
  HepLorentzVector _UpsLab;
  HepLorentzVector _PtotLab;
  HepLorentzVector _PtotUps;
  HepLorentzVector _PmissLab;
  HepLorentzVector _PmissUps;
  BbrPointErr _beamSpot;
  Hep3Vector _Thrust_Axis;
  double _Thrust_Mag;
  double _firstLepMax;
  double _secondLepMax;

};



#endif








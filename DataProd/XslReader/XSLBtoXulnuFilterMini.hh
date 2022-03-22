//--------------------------------------------------------------------------
//                                                                        //
//   Filter module for the B -> pi/pi0/eta/rho/rho0/omega l nu analyses   //
//                                                                        //
//     Sylvie Brunet  2003 Universite de Montreal (BaBar)                 //
//     David Cote     2003 Universite de Montreal                         //
//                                                                        //
//--------------------------------------------------------------------------

#ifndef XSLBTOXULNUFILTERMINI_HH
#define XSLBTOXULNUFILTERMINI_HH

//----------------------
// Base Class Headers --
//----------------------
#include "Framework/AppModule.hh"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------

class HepHistogram;
class HepTuple;

class AbsParmIfdStrKey;
class AbsParmDouble;
class AbsParmBool;
class AbsParmString;

//Necessaire d'inclure BtaCandidate pour pouvoir avoir des objets de la classe
//tels que HepLorentzVector, HepPoint, etc.
#include "Beta/BtaCandidate.hh"
class AbsEvent;

class BtaMcAssoc;
class EventInfo;
class BtaBooster;


//		---------------------
// 		-- Class Interface --
//		---------------------
 

class XSLBtoXulnuFilterMini : public AppModule {

//--------------------
// Instance Members --
//--------------------

public:

    // Constructors
    XSLBtoXulnuFilterMini( const char* const theName, const char* const theDescription );

    // Destructor
    virtual ~XSLBtoXulnuFilterMini( );

    // Operations
    virtual AppResult           beginJob( AbsEvent* anEvent );
    virtual AppResult           event   ( AbsEvent* anEvent );
    virtual AppResult           endJob  ( AbsEvent* anEvent );
      
  
protected:
  double mB0;
  double mB;

  AbsParmIfdStrKey*  _eventInfoList;
  AbsParmIfdStrKey*  _InputLeptonList;
  AbsParmIfdStrKey*  _InputPiList;
  AbsParmIfdStrKey*  _InputPi0List;
  AbsParmIfdStrKey*  _InputEtaList;
  AbsParmIfdStrKey*  _InputRho0List;
  AbsParmIfdStrKey*  _InputRhoCList;
  AbsParmIfdStrKey*  _InputOmegaList;
  AbsParmIfdStrKey*  _InputGammaList;
  AbsParmIfdStrKey*  _OutputYList;

  HepAList<BtaCandidate>* InputLeptonList;
  HepAList<BtaCandidate>* InputPiList;
  HepAList<BtaCandidate>* InputPi0List;
  HepAList<BtaCandidate>* InputEtaList;
  HepAList<BtaCandidate>* InputRho0List;
  HepAList<BtaCandidate>* InputRhoCList;
  HepAList<BtaCandidate>* InputOmegaList;
  HepAList<BtaCandidate>* InputGammaList;

  HepAList<BtaCandidate> NewLepList;
  HepAList<BtaCandidate> NewPiList;
  HepAList<BtaCandidate> NewPi0List;
  HepAList<BtaCandidate> NewEtaList;
  HepAList<BtaCandidate> NewRho0List;
  HepAList<BtaCandidate> NewRhoCList;
  HepAList<BtaCandidate> NewOmegaList;
  HepAList<BtaCandidate> NewGammaList;

  HepAList<BtaCandidate>* _eTightLH;
  HepAList<BtaCandidate>* _piLooseLH;
  HepAList<BtaCandidate>* _muTight;
  HepAList<BtaCandidate>* _muVeryTight;


private:

  void MakeHadronsAndLeptonsList( AbsEvent* anEvent );
  bool SurvivedISH( AbsEvent* anEvent );
  bool SurvivedR2All( AbsEvent* anEvent );
  bool pionPID( BtaCandidate* candid );
  bool electronPID( BtaCandidate* candid );
  bool muonPID( BtaCandidate* candid );
  virtual AppResult QuitEvent();

  //evt cuts
  AbsParmDouble*    _r2allMAX;
  AbsParmDouble*    _ishMIN;

  //trk building cuts
  AbsParmDouble*  _pLABmin_electron;
  AbsParmDouble*  _pLABmin_muon;
  AbsParmDouble*  _thetaLABmin_electron;
  AbsParmDouble*  _thetaLABmin_muon;
  AbsParmDouble*  _thetaLABmax_electron;
  AbsParmDouble*  _thetaLABmax_muon;
  AbsParmDouble*  _pLABmax_pi0;
  AbsParmDouble*  _pLABmax_eta;
  AbsParmDouble*  _pLABmax_rhoC;
  AbsParmDouble*  _pLABmax_rho0;
  AbsParmDouble*  _pLABmax_omega;
  AbsParmDouble*  _pLABmax_gamma;

  //Y selection
  AbsParmDouble*  _cosBYmin;
  AbsParmDouble*  _cosBYmax;
  AbsParmDouble*  _pStar2D_slopePseudoScalar;
  AbsParmDouble*  _pStar2D_slopePseudoVector;
  AbsParmDouble*  _pStar2D_sumPseudoScalar;
  AbsParmDouble*  _pStar2D_sumPseudoVector;

  //bool
  AbsParmBool* _analyzeChargedPi;
  AbsParmBool* _analyzeNeutralPi;
  AbsParmBool* _analyzeEta;
  AbsParmBool* _analyzeChargedRho;
  AbsParmBool* _analyzeNeutralRho;
  AbsParmBool* _analyzeOmega;
  AbsParmBool* _analyzeGamma;
  AbsParmBool* _Verbose;
  AbsParmBool* _doPionPidYourself;

  BtaBooster* _CMS;  //created and deleted in event (after SurvivedEventCuts() )

  //For Verbose stuff...
  float  _TotalNbEvents;
  float  _NbAfterISHCuts;
  float  _NbAfterR2AllCuts;
  float  _NbAfterOneLeptonTight;
  float  _NbAfterEventCuts;
  float  _NbAllModesAfterCoupleCuts;
  float  _NbPiAfterCoupleCuts;
  float  _NbPiAfterLeptonSelCuts;
  float  _NbPiAfterCosBY;
  float  _NbPiAfterpStar2D;
  float  _NbPi0AfterCoupleCuts;
  float  _NbPi0AfterLeptonSelCuts;
  float  _NbPi0AfterCosBY;
  float  _NbPi0AfterpStar2D;
  float  _NbEtaAfterCoupleCuts;
  float  _NbEtaAfterLeptonSelCuts;
  float  _NbEtaAfterCosBY;
  float  _NbEtaAfterpStar2D;
  float  _NbRhoCAfterCoupleCuts;
  float  _NbRhoCAfterLeptonSelCuts;
  float  _NbRhoCAfterCosBY;
  float  _NbRhoCAfterpStar2D;
  float  _NbRho0AfterCoupleCuts;
  float  _NbRho0AfterLeptonSelCuts;
  float  _NbRho0AfterCosBY;
  float  _NbRho0AfterpStar2D;
  float  _NbOmeAfterCoupleCuts;
  float  _NbOmeAfterLeptonSelCuts;
  float  _NbOmeAfterCosBY;
  float  _NbOmeAfterpStar2D;
  float  _NbGammaAfterCoupleCuts;
  float  _NbGammaAfterLeptonSelCuts;
  float  _NbGammaAfterCosBY;
  float _NbCoupleAllModesAfterAllCuts;
  float _NbCouplePiAfterAllCuts;
  float _NbCouplePi0AfterAllCuts;
  float _NbCoupleEtaAfterAllCuts;
  float _NbCoupleRhoCAfterAllCuts;
  float _NbCoupleRho0AfterAllCuts;
  float _NbCoupleOmegaAfterAllCuts;
  float _NbCoupleGammaAfterAllCuts;
};




#endif








//--------------------------------------------------------------------------
//                                                                        //
//   Filter module for the B -> pi/pi0/eta/rho/rho0/omega l nu analyses   //
//                                                                        //
//     Sylvie Brunet  2003 Universite de Montreal (BaBar)                 //
//     David Cote     2003 Universite de Montreal                         //
//                                                                        //
//--------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

//This module is build in an non-optimized way in order to reproduce the preliminary results presented here:
//http://www.slac.stanford.edu/~brunet/Presentations/Vub/SKIM/Aug4_2003/cut_prod_details_4August2003.pdf
//and obtained with an older code also written by us (BruCo).
//Further optimized version of this module will be in directory 5SkimAnal14 and higher.


#include "XslReader/XSLBtoXulnuFilterMini.hh"

///////////////
#include "CLHEP/Alist/AList.h"
#include "CLHEP/Alist/AIterator.h"
#include "AbsParm/AbsParmIfdStrKey.hh"
#include "AbsEnv/AbsEnv.hh"
#include "AbsEvent/AbsEvent.hh"
#include "AbsEvent/getTmpAList.hh"
#include "ErrLogger/ErrLog.hh"
#include "Framework/AbsParmBool.hh"
#include "Framework/AbsParmDouble.hh"
#include "AbsEventTag/AbsEventTag.hh"
#include "BetaCoreTools/BtaOpMakeTree.hh"
#include "Beta/EventInfo.hh"
#include "BetaCoreTools/BtaBooster.hh"
using std::cout;
using std::endl;


XSLBtoXulnuFilterMini::XSLBtoXulnuFilterMini( const char* const theName, 
					  const char* const theDescription )
  : AppModule( theName, theDescription )
  , mB0(5.2794), mB(5.2790)

  //I/O lists
  , _eventInfoList( new AbsParmIfdStrKey( "eventInfoList", this, "Default" ) )
  , _InputLeptonList(new AbsParmIfdStrKey("inputLeptonList",this,"eLHBremGTL"))
  , _InputPiList(new AbsParmIfdStrKey("inputPiList",this,"GoodTracksLoose"))
  , _InputPi0List(new AbsParmIfdStrKey("inputPi0List",this,"pi0AllVeryLoose")) 
  , _InputEtaList(new AbsParmIfdStrKey("inputEtaList",this,"etaAllSpecial3Pi"))
  , _InputRho0List(new AbsParmIfdStrKey("inputRho0List",this,"rho0VeryLooseMassAndPid"))
  , _InputRhoCList(new AbsParmIfdStrKey("inputRhoCList",this,"rhoCVeryLooseMassAndPid"))  
  , _InputOmegaList(new AbsParmIfdStrKey("inputOmegaList",this,"omegaVeryLooseAndPid"))
  , _InputGammaList(new AbsParmIfdStrKey("inputGammaList",this,"GoodPhotonDefault"))
  , _OutputYList(new AbsParmIfdStrKey("outputYList",this,"XSLBtoXulnuSkimmedYlist")) 

  //event parameters
  , _r2allMAX(new AbsParmDouble("r2allMAX",this,0.5))//no cut:>=1
  , _ishMIN(new AbsParmDouble("ishMIN",this,1))//no cut:0, cut:1
  //The one lepton tight requirement in the evt is fulfilled automatically by the taking tight electron and muon lists

  //list building parameters
  , _pLABmin_electron(new AbsParmDouble("pLABmin_electron",this,0.5))//no cut:0
  , _pLABmin_muon(new AbsParmDouble("pLABmin_muon",this,1))//no cut:0
  , _thetaLABmin_electron(new AbsParmDouble("thetaLABmin_electron",this,0.41))//no cut:0
  , _thetaLABmin_muon(new AbsParmDouble("thetaLABmin_muon",this,0.41))//no cut:0
  , _thetaLABmax_electron(new AbsParmDouble("thetaLABmax_electron",this, 2.54))//no cut: 3.14
  , _thetaLABmax_muon(new AbsParmDouble("thetaLABmax_muon",this, 2.54))//no cut: 3.14
  , _pLABmax_pi0(new AbsParmDouble("pLABmax_pi0",this,10))//no cut:100
  , _pLABmax_eta(new AbsParmDouble("pLABmax_eta",this,10))//no cut:100
  , _pLABmax_rhoC(new AbsParmDouble("pLABmax_rhoC",this,10))//no cut:100
  , _pLABmax_rho0(new AbsParmDouble("pLABmax_rho0",this,10))//no cut:100
  , _pLABmax_omega(new AbsParmDouble("pLABmax_omega",this,10))//no cut:100
  , _pLABmax_gamma(new AbsParmDouble("pLABmax_gamma",this,10))//no cut:100

  //Y selection
  //The cosBY cut is the same for all hadrons so we don't duplicate parameters
  //The cosBY cut is also the same for all leptons, since we're using BremRecovered electrons, so we don't duplicate parameters
  , _cosBYmin(new AbsParmDouble("cosBYmin",this,-1.5))//no cut:-5000
  , _cosBYmax(new AbsParmDouble("cosBYmax",this,1.3))//no cut:+1000
  , _pStar2D_slopePseudoScalar(new AbsParmDouble("pStar2D_slopePseudoScalar",this,1))
  , _pStar2D_slopePseudoVector(new AbsParmDouble("pStar2D_slopePseudoVector",this,0.696969697)) // 0.6969697=23/33
  , _pStar2D_sumPseudoScalar(new AbsParmDouble("pStar2D_sumPseudoScalar",this,2.6))//no cut: 1000
  , _pStar2D_sumPseudoVector(new AbsParmDouble("pStar2D_sumPseudoVector",this,2.3))//no cut: 1000

  //modes to be studied
  , _analyzeChargedPi(new AbsParmBool("AnalyzeChargedPi",this,true))
  , _analyzeNeutralPi(new AbsParmBool("AnalyzeNeutralPi",this,true))
  , _analyzeEta(new AbsParmBool("AnalyzeEta",this,true))
  , _analyzeChargedRho(new AbsParmBool("AnalyzeChargedRho",this,true))
  , _analyzeNeutralRho(new AbsParmBool("AnalyzeNeutralRho",this,true))
  , _analyzeOmega(new AbsParmBool("AnalyzeOmega",this,true))
  , _analyzeGamma(new AbsParmBool("AnalyzeGamma",this,true))

  //bool
  , _Verbose(new AbsParmBool("Verbose",this,false))
  , _doPionPidYourself(new AbsParmBool("doPionPidYourself",this,false))
{
  //I/O lists
  commands()->append( _eventInfoList);
  commands()->append( _InputLeptonList);
  commands()->append( _InputPiList);
  commands()->append( _InputPi0List);
  commands()->append( _InputEtaList);
  commands()->append( _InputRho0List);
  commands()->append( _InputRhoCList);
  commands()->append( _InputOmegaList);
  commands()->append( _InputGammaList);
  commands()->append( _OutputYList);

  //evt cuts
  commands()->append( _r2allMAX );
  commands()->append( _ishMIN );

  //trk building cuts
  commands()->append( _pLABmin_electron );
  commands()->append( _pLABmin_muon );
  commands()->append( _thetaLABmin_electron );
  commands()->append( _thetaLABmin_muon );
  commands()->append( _thetaLABmax_electron );
  commands()->append( _thetaLABmax_muon );
  commands()->append( _pLABmax_pi0 );
  commands()->append( _pLABmax_eta );
  commands()->append( _pLABmax_rhoC );
  commands()->append( _pLABmax_rho0 );
  commands()->append( _pLABmax_omega );
  commands()->append( _pLABmax_gamma );

  //Y selection
  commands()->append( _cosBYmin );
  commands()->append( _cosBYmax );
  commands()->append( _pStar2D_slopePseudoScalar );
  commands()->append( _pStar2D_slopePseudoVector );
  commands()->append( _pStar2D_sumPseudoScalar );
  commands()->append( _pStar2D_sumPseudoVector );

  //bool
  commands()->append( _analyzeChargedPi );
  commands()->append( _analyzeNeutralPi );
  commands()->append( _analyzeEta );
  commands()->append( _analyzeChargedRho );
  commands()->append( _analyzeNeutralRho );
  commands()->append( _analyzeOmega );
  commands()->append( _analyzeGamma );
  commands()->append( _Verbose );
  commands()->append( _doPionPidYourself );
}

XSLBtoXulnuFilterMini::~XSLBtoXulnuFilterMini()
{
  //I/O lists
  delete  _eventInfoList;
  delete  _InputLeptonList;
  delete  _InputPiList;
  delete  _InputPi0List;
  delete  _InputEtaList;
  delete  _InputRho0List;
  delete  _InputRhoCList;
  delete  _InputOmegaList;
  delete  _InputGammaList;
  delete  _OutputYList;

  //evt
  delete  _r2allMAX;
  delete  _ishMIN;

  //trk building
  delete _pLABmin_electron;
  delete _pLABmin_muon;
  delete _thetaLABmin_electron;
  delete _thetaLABmin_muon;
  delete _thetaLABmax_electron;
  delete _thetaLABmax_muon;
  delete _pLABmax_pi0;
  delete _pLABmax_eta;
  delete _pLABmax_rhoC;
  delete _pLABmax_rho0;
  delete _pLABmax_omega;
  delete _pLABmax_gamma;

  //Y selection
  delete _cosBYmin;
  delete _cosBYmax;
  delete _pStar2D_slopePseudoScalar;
  delete _pStar2D_slopePseudoVector;
  delete _pStar2D_sumPseudoScalar;
  delete _pStar2D_sumPseudoVector;

  //bool
  delete _analyzeChargedPi;
  delete _analyzeNeutralPi;
  delete _analyzeEta;
  delete _analyzeChargedRho;
  delete _analyzeNeutralRho;
  delete _analyzeOmega;
  delete _analyzeGamma;
  delete _Verbose;
  delete _doPionPidYourself;
}


AppResult
XSLBtoXulnuFilterMini::beginJob( AbsEvent* anEvent )
{
  ErrMsg(routine)<<"begin filter Job"<<endmsg; 

  _TotalNbEvents=0;
  _NbAfterISHCuts=0;
  _NbAfterR2AllCuts=0;
  _NbAfterOneLeptonTight=0;
  _NbAfterEventCuts=0;
  _NbAllModesAfterCoupleCuts=0;
  _NbPiAfterCoupleCuts=0;
  _NbPiAfterLeptonSelCuts=0;
  _NbPiAfterCosBY=0;
  _NbPiAfterpStar2D=0;
  _NbPi0AfterCoupleCuts=0;
  _NbPi0AfterLeptonSelCuts=0;
  _NbPi0AfterCosBY=0;
  _NbPi0AfterpStar2D=0;
  _NbEtaAfterCoupleCuts=0;
  _NbEtaAfterLeptonSelCuts=0;
  _NbEtaAfterCosBY=0;
  _NbEtaAfterpStar2D=0;
  _NbRhoCAfterCoupleCuts=0;
  _NbRhoCAfterLeptonSelCuts=0;
  _NbRhoCAfterCosBY=0;
  _NbRhoCAfterpStar2D=0;
  _NbRho0AfterCoupleCuts=0;
  _NbRho0AfterLeptonSelCuts=0;
  _NbRho0AfterCosBY=0;
  _NbRho0AfterpStar2D=0;
  _NbOmeAfterCoupleCuts=0;
  _NbOmeAfterLeptonSelCuts=0;
  _NbOmeAfterCosBY=0;
  _NbOmeAfterpStar2D=0;
  _NbGammaAfterCoupleCuts=0;
  _NbGammaAfterLeptonSelCuts=0;
  _NbGammaAfterCosBY=0;
  _NbCoupleAllModesAfterAllCuts=0;
  _NbCouplePiAfterAllCuts=0;
  _NbCouplePi0AfterAllCuts=0;
  _NbCoupleEtaAfterAllCuts=0;
  _NbCoupleRhoCAfterAllCuts=0;
  _NbCoupleRho0AfterAllCuts=0;
  _NbCoupleOmegaAfterAllCuts=0;
  _NbCoupleGammaAfterAllCuts=0;
  return AppResult::OK;
}


AppResult
XSLBtoXulnuFilterMini::event( AbsEvent* anEvent )
{

  //First thing: create the output list
  //The getTmpAList function makes the list available for all AppModule's children within an event
  //A zero length YList will tell to the reader to reject the event
  HepAList< BtaCandidate >* outYList;
  getTmpAList( anEvent, outYList, _OutputYList->value() );
  assert( outYList!=0 );

  _TotalNbEvents+=1;

  //The _CMS booster object is not created at this point so we cannot use QuitEvent()!
  if (SurvivedISH( anEvent ))  _NbAfterISHCuts+=1; 
  else return AppResult::OK; 
  if(SurvivedR2All( anEvent ))  _NbAfterR2AllCuts+=1; 
  else return AppResult::OK;

  // 4-Vector P(e+)+P(e-) CM Upsilon(4S)
  HepAList< EventInfo >* infoList = Ifd<HepAList< EventInfo > >::get(anEvent, _eventInfoList->value());
  if (infoList == 0){ ErrMsg(fatal) << "Could not locate event info list of name ";}
  EventInfo* eventInfo = infoList->first();  
  //HepLorentzVector UpsP4 = eventInfo->cmFrame();
  HepLorentzVector UpsP4(-0.1102553491300846,0.0,5.8772159072547439,12.101750056479895); //this way there is no problem for cosBY with OffPeak
  double s = UpsP4.mag2(); 
  BbrPointErr beamSpot=eventInfo->beamSpot();

  //Create the Ups(4S) frame booster
  _CMS = new BtaBooster( UpsP4 );

  //Making lists
  MakeHadronsAndLeptonsList( anEvent );

  if(NewLepList.length()!=0) 
    { 
      _NbAfterEventCuts+=1;  
      _NbAfterOneLeptonTight+=1;    
    }
  else return QuitEvent(); 

  //////////////////////////////////////
  //  We now build the YList made out of all Y couples passing their mode-specific criteria
  /////////////////////////////////////
  BtaCandidate* lepton;
  BtaCandidate* hadron;

  HepAListIterator<BtaCandidate> iter_lep(NewLepList);  
  HepAListIterator<BtaCandidate> iter_pi(NewPiList);  
  HepAListIterator<BtaCandidate> iter_pi0(NewPi0List);  
  HepAListIterator<BtaCandidate> iter_eta(NewEtaList);  
  HepAListIterator<BtaCandidate> iter_rho0(NewRho0List);  
  HepAListIterator<BtaCandidate> iter_rhoC(NewRhoCList);  
  HepAListIterator<BtaCandidate> iter_ome(NewOmegaList);  
  HepAListIterator<BtaCandidate> iter_gam(NewGammaList);  


  int piLevel=0,pi0Level=0,etaLevel=0,rhoCLevel=0,rho0Level=0,omeLevel=0,gamLevel=0;
  
  while ( 0 != ( lepton = iter_lep() ) ) {
   
    ///////////////////////////////////////////////////////////////////////    
    if( _analyzeChargedPi->value()==true ) {
      iter_pi.rewind();

      while ( 0 != ( hadron = iter_pi() ) ) {	
	

	if( !(hadron->overlaps(*lepton)) && lepton->charge()!=hadron->charge() ) 
	  { 
	    if( piLevel==0 ) { piLevel=1;  _NbPiAfterLeptonSelCuts+=1; }
	    
	    BtaCandidate boostedLep = _CMS->boostTo( *lepton );
	    BtaCandidate boostedHad = _CMS->boostTo( *hadron );
	    double sum = boostedLep.p() + _pStar2D_slopePseudoScalar->value()*boostedHad.p();

	    HepLorentzVector boostedP4Y = boostedLep.p4() + boostedHad.p4();
	    HepLorentzVector P4Y = lepton->p4() + hadron->p4();
	    
	    double numerat = sqrt(s)*boostedP4Y.e()-(mB0*mB0)- P4Y.mag2();  
	    double denomin = 2*sqrt( fabs(s/4-(mB0*mB0)) )*boostedP4Y.vect().mag(); 
	    //the abs is for the offres case. It has no effect otherwise!
	    double cosBY = numerat/denomin;

	    if( cosBY>=_cosBYmin->value() &&  cosBY<=_cosBYmax->value() ) 
	      {
		if( piLevel==1 ) { piLevel=2;  _NbPiAfterCosBY+=1; }

		if( sum >= _pStar2D_sumPseudoScalar->value() )		
		  {
		    if( piLevel==2 ) { piLevel=3;  _NbPiAfterpStar2D+=1;  _NbPiAfterCoupleCuts+=1; }
		    
		    BtaOpMakeTree comb;
		    BtaCandidate* Y = comb.create(*hadron,*lepton);
		    outYList->append( new BtaCandidate(*Y) );
		    _NbCouplePiAfterAllCuts+=1;
		    delete Y;
		  }//pStar2D
	      }//cosBY
	  }//opposite charge and no overlap
      }//while ( 0 != ( hadron = iter_pi() ) )
    }//if(_chargedPi)

    ///////////////////////////////////////////////////////////////////////
    if( _analyzeNeutralPi->value()==true ) {
      iter_pi0.rewind();

      while ( 0 != ( hadron = iter_pi0() ) ) {
		
	if( true )
	  {
	    if( pi0Level==0 ) { pi0Level=1;  _NbPi0AfterLeptonSelCuts+=1; }
	    
	    BtaCandidate boostedLep = _CMS->boostTo( *lepton );
	    BtaCandidate boostedHad = _CMS->boostTo( *hadron );
	    double sum = boostedLep.p() + _pStar2D_slopePseudoScalar->value()*boostedHad.p();
	    
	    HepLorentzVector boostedP4Y = boostedLep.p4() + boostedHad.p4();
	    HepLorentzVector P4Y = lepton->p4() + hadron->p4();
	    
	    double numerat = sqrt(s)*boostedP4Y.e()-(mB0*mB0)- P4Y.mag2();  
	    double denomin = 2*sqrt( fabs(s/4-(mB0*mB0)) )*boostedP4Y.vect().mag(); 
	    //the abs is for the offres case. It has no effect otherwise!
	    double cosBY = numerat/denomin;
	    
	    if( cosBY>=_cosBYmin->value() &&  cosBY<=_cosBYmax->value() ) 
	      {
		if( pi0Level==1 ) { pi0Level=2;  _NbPi0AfterCosBY+=1; }
		
		if( sum >= _pStar2D_sumPseudoScalar->value() )		
		  {
		    if( pi0Level==2 ) { pi0Level=3;  _NbPi0AfterpStar2D+=1;   _NbPi0AfterCoupleCuts+=1; }
			
		    BtaOpMakeTree comb;
		    BtaCandidate* Y = comb.create(*hadron,*lepton);
		    outYList->append( new BtaCandidate(*Y) );
		    _NbCouplePi0AfterAllCuts+=1;
		    delete Y;
		  }//pStar2D
	      }//cosBY
	  }//lepton selection
      }//while ( 0 != ( hadron = iter_pi0() ) )
    }//if(_analyzeNeutralPi)
    
    ///////////////////////////////////////////////////////////////////////////////
    if( _analyzeEta->value()==true ) {
      iter_eta.rewind();

      while ( 0 != ( hadron = iter_eta() ) ) {	
	
	//We don't want Eta's daughters to be the lepton!
	bool HadNotMixedToLepton = true;
	BtaCandidate* temp;
	HepAListIterator<BtaCandidate> fille = hadron->daughterIterator();
	while (  0 != ( temp = fille() ) )
	  {  if( temp->uid()==lepton->uid() ) HadNotMixedToLepton = false; }
	
	
	if( HadNotMixedToLepton )
	  {
	    if( etaLevel==0 ) { etaLevel=1;   _NbEtaAfterLeptonSelCuts+=1; }
	    
	    BtaCandidate boostedLep = _CMS->boostTo( *lepton );
	    BtaCandidate boostedHad = _CMS->boostTo( *hadron );
	    double sum = boostedLep.p() + _pStar2D_slopePseudoScalar->value()*boostedHad.p();
	    
	    HepLorentzVector boostedP4Y = boostedLep.p4() + boostedHad.p4();
	    HepLorentzVector P4Y = lepton->p4() + hadron->p4();
	    
	    double numerat = sqrt(s)*boostedP4Y.e()-(mB0*mB0)- P4Y.mag2();  
	    double denomin = 2*sqrt( fabs(s/4-(mB0*mB0)) )*boostedP4Y.vect().mag(); 
	    //the abs is for the offres case. It has no effect otherwise!
	    double cosBY = numerat/denomin;

	    if( cosBY>=_cosBYmin->value() &&  cosBY<=_cosBYmax->value() ) 
	      {
		if( etaLevel==1 ) { etaLevel=2;  _NbEtaAfterCosBY+=1; }

		if( sum >= _pStar2D_sumPseudoScalar->value() )		
		  {
		    if( etaLevel==2 ) { etaLevel=3;  _NbEtaAfterpStar2D+=1;  _NbEtaAfterCoupleCuts+=1; }
		    
			BtaOpMakeTree comb;
			BtaCandidate* Y = comb.create(*hadron,*lepton);
			outYList->append( new BtaCandidate(*Y) );
			_NbCoupleEtaAfterAllCuts+=1;
			delete Y;
		  }//pStar2D
	      }//cosBY
	  }//if HadNotMixedToLepton
      }//while ( 0 != ( hadron = iter_eta() ) )
    }//if(_analyzeEta)


    ///////////////////////////////////////////////////////////////////////////////
    if( _analyzeChargedRho->value()==true ) {
      iter_rhoC.rewind();

      while ( 0 != ( hadron = iter_rhoC() ) ) {	
	
	//We don't want Rho's daughters to be the lepton!
	bool HadNotMixedToLepton = true;
	BtaCandidate* temp;
	HepAListIterator<BtaCandidate> fille = hadron->daughterIterator();
	while (  0 != ( temp = fille() ) )
	  {  if( temp->uid()==lepton->uid() ) HadNotMixedToLepton = false; }

		
	if( lepton->charge()!=hadron->charge() && HadNotMixedToLepton )
	  {
	    if( rhoCLevel==0 ) { rhoCLevel=1;  _NbRhoCAfterLeptonSelCuts+=1; }

	    BtaCandidate boostedLep = _CMS->boostTo( *lepton );
	    BtaCandidate boostedHad = _CMS->boostTo( *hadron );
	    double sum = boostedLep.p() + _pStar2D_slopePseudoVector->value()*boostedHad.p();

	    HepLorentzVector boostedP4Y = boostedLep.p4() + boostedHad.p4();
	    HepLorentzVector P4Y = lepton->p4() + hadron->p4();
	    
	    double numerat = sqrt(s)*boostedP4Y.e()-(mB0*mB0)- P4Y.mag2();  
	    double denomin = 2*sqrt( fabs(s/4-(mB0*mB0)) )*boostedP4Y.vect().mag(); 
	    //the abs is for the offres case. It has no effect otherwise!
	    double cosBY = numerat/denomin;

	    if( cosBY>=_cosBYmin->value() &&  cosBY<=_cosBYmax->value() ) 
	      {
		if( rhoCLevel==1 ) { rhoCLevel=2;  _NbRhoCAfterCosBY+=1; }

		if( sum >= _pStar2D_sumPseudoVector->value() )		
		  {
		    if( rhoCLevel==2 ) { rhoCLevel=3;  _NbRhoCAfterpStar2D+=1;   _NbRhoCAfterCoupleCuts+=1; }
		    
		    BtaOpMakeTree comb;
		    BtaCandidate* Y = comb.create(*hadron,*lepton);
		    outYList->append( new BtaCandidate(*Y) );
		    _NbCoupleRhoCAfterAllCuts+=1;
		    delete Y;
		  }//pStar2D
	      }//cosBY
	  }//if opposite charges and HadNotMixedToLepton
      }//while ( 0 != ( hadron = iter_rhoC() ) )
    }//if(_analyzeChargedRho)

    ///////////////////////////////////////////////////////////////////////////////
    if( _analyzeNeutralRho->value()==true ) {
      iter_rho0.rewind();

      while ( 0 != ( hadron = iter_rho0() ) ) {	
	
	//We don't want Rho's daughters to be the lepton!
	bool HadNotMixedToLepton = true;
	BtaCandidate* temp;
	HepAListIterator<BtaCandidate> fille = hadron->daughterIterator();
	while (  0 != ( temp = fille() ) )
	  {  if( temp->uid()==lepton->uid() ) HadNotMixedToLepton = false; }

	
	if( HadNotMixedToLepton )
	  {
	    if( rho0Level==0 ) { rho0Level=1;  _NbRho0AfterLeptonSelCuts+=1; }

	    BtaCandidate boostedLep = _CMS->boostTo( *lepton );
	    BtaCandidate boostedHad = _CMS->boostTo( *hadron );
	    double sum = boostedLep.p() + _pStar2D_slopePseudoVector->value()*boostedHad.p();

	    HepLorentzVector boostedP4Y = boostedLep.p4() + boostedHad.p4();
	    HepLorentzVector P4Y = lepton->p4() + hadron->p4();
	    
	    double numerat = sqrt(s)*boostedP4Y.e()-(mB0*mB0)- P4Y.mag2();  
	    double denomin = 2*sqrt( fabs(s/4-(mB0*mB0)) )*boostedP4Y.vect().mag(); 
	    //the abs is for the offres case. It has no effect otherwise!
	    double cosBY = numerat/denomin;

	    if( cosBY>=_cosBYmin->value() &&  cosBY<=_cosBYmax->value() ) 
	      {
		if( rho0Level==1 ) { rho0Level=2;  _NbRho0AfterCosBY+=1; }

		if( sum >= _pStar2D_sumPseudoVector->value() )		
		  {
		    if( rho0Level==2 ) { rho0Level=3;  _NbRho0AfterpStar2D+=1;  _NbRho0AfterCoupleCuts+=1; }
		    
		    BtaOpMakeTree comb;
		    BtaCandidate* Y = comb.create(*hadron,*lepton);
		    outYList->append( new BtaCandidate(*Y) );
		    _NbCoupleRho0AfterAllCuts+=1;
		    delete Y;
		  }//pStar2D
	      }//cosBY
	  }//if HadNotMixedToLepton
      }//while ( 0 != ( hadron = iter_rho0() ) )
    }//if(_analyzeNeutralRho)

    ///////////////////////////////////////////////////////////////////////////////
    if( _analyzeOmega->value()==true ) {
      iter_ome.rewind();

      while ( 0 != ( hadron = iter_ome() ) ) {	
	
	//We don't want Omega's daughters to be the lepton!
	bool HadNotMixedToLepton = true;
	BtaCandidate* temp;
	HepAListIterator<BtaCandidate> fille = hadron->daughterIterator();
	while (  0 != ( temp = fille() ) )
	  {  if( temp->uid()==lepton->uid() ) HadNotMixedToLepton = false; }

	
	if( HadNotMixedToLepton )
	  {
	    if( omeLevel==0 ) { omeLevel=1;  _NbOmeAfterLeptonSelCuts+=1; }

	    BtaCandidate boostedLep = _CMS->boostTo( *lepton );
	    BtaCandidate boostedHad = _CMS->boostTo( *hadron );
	    double sum = boostedLep.p() + _pStar2D_slopePseudoVector->value()*boostedHad.p();

	    HepLorentzVector boostedP4Y = boostedLep.p4() + boostedHad.p4();
	    HepLorentzVector P4Y = lepton->p4() + hadron->p4();
	    
	    double numerat = sqrt(s)*boostedP4Y.e()-(mB0*mB0)- P4Y.mag2();  
	    double denomin = 2*sqrt( fabs(s/4-(mB0*mB0)) )*boostedP4Y.vect().mag(); 
	    //the abs is for the offres case. It has no effect otherwise!
	    double cosBY = numerat/denomin;

	    if( cosBY>=_cosBYmin->value() &&  cosBY<=_cosBYmax->value() ) 
	      {
		if( omeLevel==1 ) { omeLevel=2;  _NbOmeAfterCosBY+=1; }

		if( sum >= _pStar2D_sumPseudoVector->value() )		
		  {
		    if( omeLevel==2 ) { omeLevel=3;  _NbOmeAfterpStar2D+=1;   _NbOmeAfterCoupleCuts+=1; }
		    
		    BtaOpMakeTree comb;
		    BtaCandidate* Y = comb.create(*hadron,*lepton);
		    outYList->append( new BtaCandidate(*Y) );
		    _NbCoupleOmegaAfterAllCuts+=1;
		    delete Y;
		  }//pStar2D
	      }//cosBY
	  }//if HadNotMixedToLepton
      }//while ( 0 != ( hadron = iter_ome() ) )
    }//if(_analyzeOmega)

    ///////////////////////////////////////////////////////////////////////////////
    if( _analyzeGamma->value()==true ) {
      iter_gam.rewind();
      
      while ( 0 != ( hadron = iter_gam() ) ) {	
	
	  if( gamLevel==0 ) { gamLevel=1;  _NbGammaAfterLeptonSelCuts+=1; }

	  BtaCandidate boostedLep = _CMS->boostTo( *lepton );
	  BtaCandidate boostedHad = _CMS->boostTo( *hadron );

	  HepLorentzVector boostedP4Y = boostedLep.p4() + boostedHad.p4();
	  HepLorentzVector P4Y = lepton->p4() + hadron->p4();
	  
	  double numerat = sqrt(s)*boostedP4Y.e()-(mB0*mB0)- P4Y.mag2();  
	  double denomin = 2*sqrt( fabs(s/4-(mB0*mB0)) )*boostedP4Y.vect().mag(); 
	  //the abs is for the offres case. It has no effect otherwise!
	  double cosBY = numerat/denomin;

	  if( cosBY>=_cosBYmin->value() &&  cosBY<=_cosBYmax->value() ) 
	    {
	      if( gamLevel==1 ) { gamLevel=2;  _NbGammaAfterCosBY+=1; _NbGammaAfterCoupleCuts+=1; }
	      
	      BtaOpMakeTree comb;
	      BtaCandidate* Y = comb.create(*hadron,*lepton);
	      outYList->append( new BtaCandidate(*Y) );
	      _NbCoupleGammaAfterAllCuts+=1;
	      delete Y;
	    }//cosBY
      }//while ( 0 != ( hadron = iter_gam() ) )
    }//if(_analyzeGamma)
    
  }//while ( 0 != ( lepton = iter_lep() ) )

  if( outYList->length()!=0 ) { _NbAllModesAfterCoupleCuts+=1; _NbCoupleAllModesAfterAllCuts+=outYList->length(); }

  //To be made persistent:
  bool isPi=false,isPi0=false,isEta=false,isRhoC=false,isRho0=false,isOmega=false,isGam=false;
  if(piLevel==3) isPi=true;
  if(pi0Level==3) isPi0=true;
  if(etaLevel==3) isEta=true;
  if(rhoCLevel==3) isRhoC=true;
  if(rho0Level==3) isRho0=true;
  if(omeLevel==3) isOmega=true;
  if(gamLevel==2) isGam=true;

  return QuitEvent(); 

}


bool
XSLBtoXulnuFilterMini::SurvivedISH(AbsEvent* anEvent)
{
  bool evtOK=false;
  float ish=-1;

  AbsEventTag* tag = Ifd<AbsEventTag>::get( anEvent );
  if(tag==0) 
    { 
      cout<<"NOTHING IN AbsEvent!"<<endl;
    }
  else
    {
      bool ismultihad;
      if(tag->getBool(ismultihad,"BGFMultiHadron"))	{
	if(ismultihad) ish=1;
	else ish=0;
      }
    }
  
  if( ish>= _ishMIN->value() ) evtOK=true;
  
  return evtOK;
}

bool
XSLBtoXulnuFilterMini::SurvivedR2All(AbsEvent* anEvent)
{
  bool evtOK=false;
  float r2all=2;

  AbsEventTag* tag = Ifd<AbsEventTag>::get( anEvent );
  if(tag==0) 
    { 
      cout<<"NOTHING IN AbsEvent!"<<endl;
    }
  else
    {
      if(tag->getFloat(r2all,"R2All"));
    }  
  
  if( r2all<=_r2allMAX->value() ) evtOK=true;

  return evtOK;
}


void
XSLBtoXulnuFilterMini::MakeHadronsAndLeptonsList(AbsEvent* anEvent)
{
  //We apply here additional kinematic and quality cuts to the various lists
  //Note that the base list are already somewhat customized in various CompositionSequences modules
  BtaCandidate* candid;  


  //PID Lists   //PidLHElectrons, piLHLoose, cutBased muSelector
  getTmpAList (anEvent, _eTightLH,  IfdStrKey(HepString("PidLHElectrons")) );
  if (_eTightLH == 0)  ErrMsg(fatal) << "Could not locate list eTightLH" << endmsg;

  getTmpAList (anEvent, _piLooseLH,  IfdStrKey(HepString("piLHLoose")) );
  if (_piLooseLH == 0)  ErrMsg(fatal) << "Could not locate list piLHLoose" << endmsg;

  getTmpAList (anEvent, _muTight,  IfdStrKey(HepString("muMicroTight")) );
  if (_muTight == 0)  ErrMsg(fatal) << "Could not locate list muMicroTight" << endmsg;
  getTmpAList (anEvent, _muVeryTight,  IfdStrKey(HepString("muMicroVeryTight")) );
  if (_muVeryTight == 0)  ErrMsg(fatal) << "Could not locate list muMicroVeryTight" << endmsg;


  ///////////////
  //  LEPTONS  //
  ///////////////
  getTmpAList(anEvent,InputLeptonList,_InputLeptonList->value());
  if (InputLeptonList == 0)   ErrMsg(fatal)<<"No leptonList!!"<<endmsg;
  HepAListIterator<BtaCandidate> iter_InputLepton(*InputLeptonList);  
  while ( 0 != ( candid = iter_InputLepton()) ) 
    {
      bool ePID = electronPID( candid );
      bool muPID = muonPID( candid );

      bool eCut=false;
      bool muCut=false;
      if(candid->p() >= _pLABmin_electron->value() && ePID ) eCut=true;
      if(candid->p() >= _pLABmin_muon->value() && muPID ) muCut=true;

      if( ( eCut || muCut)	
	  && candid->p()<10  
	  && candid->p3().theta() >= _thetaLABmin_electron->value()
	  && candid->p3().theta() <= _thetaLABmax_electron->value() )  
	NewLepList.append( candid );
    }


  ///////////////
  //  HADRONS  //
  ///////////////

  //Put Hadron Selection cuts here! 

  if(_analyzeChargedPi->value()==true)
    {
      getTmpAList(anEvent,InputPiList,_InputPiList->value());
      if (InputPiList == 0)   ErrMsg(fatal)<<"No piList!!"<<endmsg;
      HepAListIterator<BtaCandidate> iter_InputPi(*InputPiList);  
      while ( 0 != ( candid = iter_InputPi()) ) 
	{
	  bool daughtersOK=false;
	  if( _doPionPidYourself->value()==false) daughtersOK=true;
	  else daughtersOK = pionPID(candid);

	  if( daughtersOK && candid->nDaughters()==0 ) NewPiList.append( candid ); 
	}
    }

  if(_analyzeNeutralPi->value()==true)
    {
      getTmpAList(anEvent,InputPi0List,_InputPi0List->value());
      if (InputPi0List == 0)   ErrMsg(fatal)<<"No pi0List!!"<<endmsg;
      HepAListIterator<BtaCandidate> iter_InputPi0(*InputPi0List);  
      while ( 0 != ( candid = iter_InputPi0()) ) 
	{
	  if( candid->p() <= _pLABmax_pi0->value() ) NewPi0List.append( candid );
	}
    }

  if(_analyzeEta->value()==true)
    {
      getTmpAList(anEvent,InputEtaList,_InputEtaList->value());
      if (InputEtaList == 0)   ErrMsg(fatal)<<"No etaList!!"<<endmsg;
      HepAListIterator<BtaCandidate> iter_InputEta(*InputEtaList);  
      while ( 0 != ( candid = iter_InputEta()) ) 
	{
	  bool daughtersOK=false;
	  if(candid->nDaughters()==2 || _doPionPidYourself->value()==false ) daughtersOK=true; //eta-> gam gam mode
	  else if( candid->nDaughters()==3 && _doPionPidYourself->value()==true ) //eta->pi+ pi- pi0 mode
	    {
	      BtaCandidate* cand;
	      int g=0;
	      HepAListIterator<BtaCandidate>iter_fille = candid->daughterIterator();
	      //We ask the two charged pions daughters to pass the piLH selector
	      while ( 0 != (cand =iter_fille()) && g>=0 ) 
		{ if(cand->charge()!=0 && !pionPID(cand) ) g=-1; }
	      if(g==0) daughtersOK=true;
	    }

	  if( daughtersOK && candid->p()<=_pLABmax_eta->value() ) NewEtaList.append( candid );
	}
    }

  if(_analyzeNeutralRho->value()==true)
    {
      getTmpAList(anEvent,InputRho0List,_InputRho0List->value());
      if (InputRho0List == 0)   ErrMsg(fatal)<<"No rho0List!!"<<endmsg;
      HepAListIterator<BtaCandidate> iter_InputRho0(*InputRho0List);  
      while ( 0 != ( candid = iter_InputRho0()) ) 
	{
	  bool daughtersOK=false;

	  if(_doPionPidYourself->value()==false) daughtersOK=true;
	  else 
	    {
	      BtaCandidate* cand;
	      int g=0;
	      HepAListIterator<BtaCandidate>iter_fille = candid->daughterIterator();
	      //We ask the two charged pions daughters to pass the piLH selector
	      while ( 0 != (cand =iter_fille()) && g>=0 ) 
		{ if(cand->charge()!=0 && !pionPID(cand) ) g=-1; }
	      if(g==0) daughtersOK=true;
	    }

	  if( daughtersOK && candid->p() <= _pLABmax_rho0->value() ) NewRho0List.append( candid );
	}
    }

  if(_analyzeChargedRho->value()==true)
    {
      getTmpAList(anEvent,InputRhoCList,_InputRhoCList->value());
      if (InputRhoCList == 0)   ErrMsg(fatal)<<"No rhoCList!!"<<endmsg;
      HepAListIterator<BtaCandidate> iter_InputRhoC(*InputRhoCList);  
      while ( 0 != ( candid = iter_InputRhoC()) ) 
	{
	  bool daughtersOK=false;

	  if(_doPionPidYourself->value()==false) daughtersOK=true;
	  else 
	    {
	      BtaCandidate* cand;
	      int g=0;
	      HepAListIterator<BtaCandidate>iter_fille = candid->daughterIterator();
	      //We ask the two charged pions daughters to pass the piLH selector
	      while ( 0 != (cand =iter_fille()) && g>=0 ) 
		{ if(cand->charge()!=0 && !pionPID(cand) ) g=-1; }
	      if(g==0) daughtersOK=true;
	    }

	  if( daughtersOK && candid->p()<=_pLABmax_rhoC->value() ) NewRhoCList.append( candid );
	}
    }

  if(_analyzeOmega->value()==true)
    {
      getTmpAList(anEvent,InputOmegaList,_InputOmegaList->value());
      if (InputOmegaList == 0)   ErrMsg(fatal)<<"No omegaList!!"<<endmsg;
      HepAListIterator<BtaCandidate> iter_InputOmega(*InputOmegaList);  
      while ( 0 != ( candid = iter_InputOmega()) ) 
	{
	  bool daughtersOK=false;

	  if(_doPionPidYourself->value()==false) daughtersOK=true;
	  else 
	    {
	      BtaCandidate* cand;
	      int g=0;
	      HepAListIterator<BtaCandidate>iter_fille = candid->daughterIterator();
	      //We ask the two charged pions daughters to pass the piLH selector
	      while ( 0 != (cand =iter_fille()) && g>=0 ) 
		{ if(cand->charge()!=0 && !pionPID(cand) ) g=-1; }
	      if(g==0) daughtersOK=true;
	    }

	  if( daughtersOK && candid->p()<=_pLABmax_omega->value() ) NewOmegaList.append( candid );
	}
    }

  if(_analyzeGamma->value()==true)
    {
      getTmpAList(anEvent,InputGammaList,_InputGammaList->value());
      if (InputGammaList == 0)   ErrMsg(fatal)<<"No gammaList!!"<<endmsg;
      HepAListIterator<BtaCandidate> iter_InputGamma(*InputGammaList);  
      while ( 0 != ( candid = iter_InputGamma()) ) 
	{
	  if( candid->p() <= _pLABmax_gamma->value() ) NewGammaList.append( candid );
	}
    } // if _analyzeGamma
  
}



AppResult
XSLBtoXulnuFilterMini::endJob( AbsEvent* anEvent )
{

  if(_Verbose->value())
    {
      cout<<"Results from B->Xu l nu XSL skim"<<endl;
      cout<<endl;
      cout<<"SUMMARY:"<<endl;
      cout<<"Total number of events before any cut:              "<<_TotalNbEvents<<endl;
      cout<<"Events passing Event cuts:                          "<<_NbAfterEventCuts<<" --> "<<_NbAfterEventCuts/_TotalNbEvents*100<<"%"<<endl;
      cout<<"Events passing Couple cuts (All Modes):             "<<_NbAllModesAfterCoupleCuts<<" --> "<<_NbAllModesAfterCoupleCuts/_TotalNbEvents*100<<"%"<<endl;
      cout<<"Events passing Couple cuts (Pi Mode):               "<<_NbPiAfterCoupleCuts<<" --> "<<_NbPiAfterCoupleCuts/_TotalNbEvents*100<<"%"<<endl;
      cout<<"Events passing Couple cuts (Pi0 Mode):              "<<_NbPi0AfterCoupleCuts<<" --> "<<_NbPi0AfterCoupleCuts/_TotalNbEvents*100<<"%"<<endl;
      cout<<"Events passing Couple cuts (Eta Mode):              "<<_NbEtaAfterCoupleCuts<<" --> "<<_NbEtaAfterCoupleCuts/_TotalNbEvents*100<<"%"<<endl;
      cout<<"Events passing Couple cuts (Rho0 Mode):             "<<_NbRho0AfterCoupleCuts<<" --> "<<_NbRho0AfterCoupleCuts/_TotalNbEvents*100<<"%"<<endl;
      cout<<"Events passing Couple cuts (RhoC Mode):             "<<_NbRhoCAfterCoupleCuts<<" --> "<<_NbRhoCAfterCoupleCuts/_TotalNbEvents*100<<"%"<<endl;
      cout<<"Events passing Couple cuts (Omega Mode):            "<<_NbOmeAfterCoupleCuts<<" --> "<<_NbOmeAfterCoupleCuts/_TotalNbEvents*100<<"%"<<endl;
      cout<<"Events passing Couple cuts (Gamma Mode):            "<<_NbGammaAfterCoupleCuts<<" --> "<<_NbGammaAfterCoupleCuts/_TotalNbEvents*100<<"%"<<endl;
      cout<<endl;
      cout<<"BtaCandidates INFO:"<<endl;
      cout<<"Nb of BtaCandidates/accepted event after all cuts (All Modes):  "<<_NbCoupleAllModesAfterAllCuts<<" --> "<<_NbCoupleAllModesAfterAllCuts/_NbAllModesAfterCoupleCuts<<endl;
      cout<<"Nb of BtaCandidates/accepted event after all cuts (pi Mode):    "<<_NbCouplePiAfterAllCuts<<" --> "<<_NbCouplePiAfterAllCuts/_NbPiAfterCoupleCuts<<endl;
      cout<<"Nb of BtaCandidates/accepted event after all cuts (pi0 Mode):   "<<_NbCouplePi0AfterAllCuts<<" --> "<<_NbCouplePi0AfterAllCuts/_NbPi0AfterCoupleCuts<<endl;
      cout<<"Nb of BtaCandidates/accepted event after all cuts (eta Mode):   "<<_NbCoupleEtaAfterAllCuts<<" --> "<<_NbCoupleEtaAfterAllCuts/_NbEtaAfterCoupleCuts<<endl;
      cout<<"Nb of BtaCandidates/accepted event after all cuts (rhoC Mode):  "<<_NbCoupleRhoCAfterAllCuts<<" --> "<<_NbCoupleRhoCAfterAllCuts/_NbRhoCAfterCoupleCuts<<endl;
      cout<<"Nb of BtaCandidates/accepted event after all cuts (rho0 Mode):  "<<_NbCoupleRho0AfterAllCuts<<" --> "<<_NbCoupleRho0AfterAllCuts/_NbRho0AfterCoupleCuts<<endl;
      cout<<"Nb of BtaCandidates/accepted event after all cuts (omega Mode): "<<_NbCoupleOmegaAfterAllCuts<<" --> "<<_NbCoupleOmegaAfterAllCuts/_NbOmeAfterCoupleCuts<<endl;
      cout<<"Nb of BtaCandidates/accepted event after all cuts (gamma Mode): "<<_NbCoupleGammaAfterAllCuts<<" --> "<<_NbCoupleGammaAfterAllCuts/_NbGammaAfterCoupleCuts<<endl;
      cout<<endl;
      cout<<"CUT BY CUT RESULTS:"<<endl;
      cout<<"Events passing ISH:                                 "<<_NbAfterISHCuts<<" --> "<<_NbAfterISHCuts/_TotalNbEvents*100<<"%"<<endl;
      cout<<"Events passing R2All:                               "<<_NbAfterR2AllCuts<<" --> "<<_NbAfterR2AllCuts/_TotalNbEvents*100<<"%"<<endl;
      cout<<"Events passing >=1 lepton tight:                    "<<_NbAfterOneLeptonTight<<" --> "<<_NbAfterOneLeptonTight/_TotalNbEvents*100<<"%"<<endl;
      cout<<endl;
      cout<<"Events passing Couple lepton Sel cuts (Pi Mode):    "<<_NbPiAfterLeptonSelCuts<<" --> "<<_NbPiAfterLeptonSelCuts/_TotalNbEvents*100<<"%"<<endl;
      cout<<"Events passing Couple cosBY cut (Pi Mode):          "<<_NbPiAfterCosBY<<" --> "<<_NbPiAfterCosBY/_TotalNbEvents*100<<"%"<<endl;
      cout<<"Events passing Couple pStar2D cut (Pi Mode):        "<<_NbPiAfterpStar2D<<" --> "<<_NbPiAfterpStar2D/_TotalNbEvents*100<<"%"<<endl;
      cout<<endl;
      cout<<"Events passing Couple lepton Sel cuts (Pi0 Mode):   "<<_NbPi0AfterLeptonSelCuts<<" --> "<<_NbPi0AfterLeptonSelCuts/_TotalNbEvents*100<<"%"<<endl;
      cout<<"Events passing Couple cosBY cut (Pi0 Mode):         "<<_NbPi0AfterCosBY<<" --> "<<_NbPi0AfterCosBY/_TotalNbEvents*100<<"%"<<endl;
      cout<<"Events passing Couple pStar2D cut (Pi0 Mode):       "<<_NbPi0AfterpStar2D<<" --> "<<_NbPi0AfterpStar2D/_TotalNbEvents*100<<"%"<<endl;
      cout<<endl;
      cout<<"Events passing Couple lepton Sel cuts (Eta Mode):   "<<_NbEtaAfterLeptonSelCuts<<" --> "<<_NbEtaAfterLeptonSelCuts/_TotalNbEvents*100<<"%"<<endl;
      cout<<"Events passing Couple cosBY cut (Eta Mode):         "<<_NbEtaAfterCosBY<<" --> "<<_NbEtaAfterCosBY/_TotalNbEvents*100<<"%"<<endl;
      cout<<"Events passing Couple pStar2D cut (Eta Mode):       "<<_NbEtaAfterpStar2D<<" --> "<<_NbEtaAfterpStar2D/_TotalNbEvents*100<<"%"<<endl;
      cout<<endl;
      cout<<"Events passing Couple lepton Sel cuts (Rho0 Mode):  "<<_NbRho0AfterLeptonSelCuts<<" --> "<<_NbRho0AfterLeptonSelCuts/_TotalNbEvents*100<<"%"<<endl;
      cout<<"Events passing Couple cosBY cut (Rho0 Mode):        "<<_NbRho0AfterCosBY<<" --> "<<_NbRho0AfterCosBY/_TotalNbEvents*100<<"%"<<endl;
      cout<<"Events passing Couple pStar2D cut (Rho0 Mode):      "<<_NbRho0AfterpStar2D<<" --> "<<_NbRho0AfterpStar2D/_TotalNbEvents*100<<"%"<<endl;
      cout<<endl;
      cout<<"Events passing Couple lepton Sel cuts (RhoC Mode):  "<<_NbRhoCAfterLeptonSelCuts<<" --> "<<_NbRhoCAfterLeptonSelCuts/_TotalNbEvents*100<<"%"<<endl;
      cout<<"Events passing Couple cosBY cut (RhoC Mode):        "<<_NbRhoCAfterCosBY<<" --> "<<_NbRhoCAfterCosBY/_TotalNbEvents*100<<"%"<<endl;
      cout<<"Events passing Couple pStar2D cut (RhoC Mode):      "<<_NbRhoCAfterpStar2D<<" --> "<<_NbRhoCAfterpStar2D/_TotalNbEvents*100<<"%"<<endl;
      cout<<endl;
      cout<<"Events passing Couple lepton Sel cuts (Omega Mode): "<<_NbOmeAfterLeptonSelCuts<<" --> "<<_NbOmeAfterLeptonSelCuts/_TotalNbEvents*100<<"%"<<endl;
      cout<<"Events passing Couple cosBY cut (Omega Mode):       "<<_NbOmeAfterCosBY<<" --> "<<_NbOmeAfterCosBY/_TotalNbEvents*100<<"%"<<endl;
      cout<<"Events passing Couple pStar2D cut (Omega Mode):     "<<_NbOmeAfterpStar2D<<" --> "<<_NbOmeAfterpStar2D/_TotalNbEvents*100<<"%"<<endl;
      cout<<endl;
      cout<<"Events passing Couple lepton Sel cuts (Gamma Mode): "<<_NbGammaAfterLeptonSelCuts<<" --> "<<_NbGammaAfterLeptonSelCuts/_TotalNbEvents*100<<"%"<<endl;
      cout<<"Events passing Couple cosBY cut (Gamma Mode):       "<<_NbGammaAfterCosBY<<" --> "<<_NbGammaAfterCosBY/_TotalNbEvents*100<<"%"<<endl;
      cout<<"No pStar2D cut for the Gamma Mode"<<endl;
      cout<<endl;
      cout<<"C't'a peu pres ca pour a soir... Byebye... Byebye..."<<endl;
    }
  
  return AppResult::OK;

}

//////////////

bool XSLBtoXulnuFilterMini::pionPID( BtaCandidate* candid )
{
  BtaCandidate* PIDcand;
  
  HepAListIterator<BtaCandidate> piL(*_piLooseLH);
  while( 0 != (PIDcand = piL() ) ){
    if( PIDcand->overlaps(*candid) )
      { 
	return true;
      }
  }
  
  return false;  
}

bool XSLBtoXulnuFilterMini::electronPID( BtaCandidate* candid )
{
  BtaCandidate* PIDcand;
  
  HepAListIterator<BtaCandidate> eT(*_eTightLH);
  while(0 != ( PIDcand = eT() ) ){
    if( PIDcand->overlaps(*candid) )
      { 
	return true;
      }
  }
  
  return false;  
}

bool XSLBtoXulnuFilterMini::muonPID( BtaCandidate* candid )
{
  BtaCandidate* PIDcand;
  
  HepAListIterator<BtaCandidate> muT(*_muTight);
  while(0 != ( PIDcand = muT()) ){
    if( PIDcand->overlaps(*candid) )
      { 
	return true;
      }
  }

  HepAListIterator<BtaCandidate> muVT(*_muVeryTight);
  while(0 != ( PIDcand = muVT()) ){
    if( PIDcand->overlaps(*candid) )
      { 
	return true;
      }
  }

    
  return false;
}


AppResult
XSLBtoXulnuFilterMini::QuitEvent()
{
  //CLEANING the list before quiting the event...
    
  NewLepList.removeAll();
  NewPiList.removeAll();
  NewPi0List.removeAll();
  NewEtaList.removeAll();
  NewOmegaList.removeAll();
  NewRho0List.removeAll();
  NewRhoCList.removeAll();
  NewGammaList.removeAll();
  
  HepAListDeleteAll( NewLepList );
  HepAListDeleteAll( NewPiList );
  HepAListDeleteAll( NewPi0List );
  HepAListDeleteAll( NewEtaList );
  HepAListDeleteAll( NewOmegaList );
  HepAListDeleteAll( NewRho0List );
  HepAListDeleteAll( NewRhoCList );
  HepAListDeleteAll( NewGammaList );
  
  delete _CMS;

  return AppResult::OK;

}

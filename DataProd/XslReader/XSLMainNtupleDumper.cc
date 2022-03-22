//------------------------------------------------------------------------------//
//                                                                              //
//  Analysis module for the B -> pi/pi0/eta/eta'/rhoC/rho0/omega l nu analyses  //
//                                                                              //
//      Sylvie Brunet  2002-2004 Universite de Montreal                         //
//      David Cote     2003-2004 Universite de Montreal                         //
//      Benoit Viaud        2004 Universite de Montreal                         //
//                                                                              //
//------------------------------------------------------------------------------//
#include "BaBar/BaBar.hh"

#include "XslReader/XSLReader.hh"

#include "Beta/EventInfo.hh"
#include "AbsEvent/AbsEvent.hh"
#include "AbsEvent/getTmpAList.hh"
#include "ErrLogger/ErrLog.hh"
#include "HepTuple/Tuple.h"
#include "HepTuple/HTValOrderedVector.h"
#include "BetaCoreTools/BtaBooster.hh"
#include "PDT/Pdt.hh"
#include "BetaCoreTools/BtaThrust.hh"

//Fit
#include "ProbTools/probab.hh"

//XSL code
#include "XslTools/XSLRotateAndBoost.hh"
#include "XslTools/XSLMCTruthAnalyzer.hh"
#include "XslTools/XSLRecoYAnalyzer.hh"
#include "XslTools/XSLYAveBooster.hh"
#include "XslFFReweighting/XSLKin.hh"

using std::cout;
using std::endl;

///////////////


void 
XSLReader::FillMainNtuple( AbsEvent* anEvent )
{
  BtaCandidate* Y;  
  HepAListIterator<BtaCandidate> iter_YList(*_YList);  
  HepAList<BtaCandidate> *PiList,*Pi0List,*Eta2List,*Eta3List,*EtapE2PPList,*EtapE3PPList,*EtapRGList,*RhoCList,*Rho0List,*OmegaList,*GammaList;
  getTmpAList (anEvent, PiList, IfdStrKey(HepString("PiList")) );
  getTmpAList (anEvent, Pi0List, IfdStrKey(HepString("Pi0List")) );
  getTmpAList (anEvent, Eta2List, IfdStrKey(HepString("Eta2List")) );
  getTmpAList (anEvent, Eta3List, IfdStrKey(HepString("Eta3List")) );
  getTmpAList (anEvent, EtapE2PPList, IfdStrKey(HepString("EtapE2PPList")) );
  getTmpAList (anEvent, EtapE3PPList, IfdStrKey(HepString("EtapE3PPList")) );
  getTmpAList (anEvent, EtapRGList, IfdStrKey(HepString("EtapRGList")) );
  getTmpAList (anEvent, RhoCList, IfdStrKey(HepString("RhoCList")) );
  getTmpAList (anEvent, Rho0List, IfdStrKey(HepString("Rho0List")) );
  getTmpAList (anEvent, OmegaList, IfdStrKey(HepString("OmegaList")) );
  getTmpAList (anEvent, GammaList, IfdStrKey(HepString("GammaList")) );

  //Fill the different mode lists before sending to the ntuple fillers
  while ( 0 != ( Y = iter_YList()) ) 
    {
      BtaCandidate* fille,*fille2;
      HepAListIterator<BtaCandidate>iter_fille = Y->daughterIterator();
      int d=0;
      while ( 0 != (fille =iter_fille()) && d==0 )   
	{  
	  //The first daughter is always the hadron by construction
	  //WARNING! this is the case only because the Y are created like this: BtaCandidate* Y = comb.create(*hadron,*lepton);
	  //The hadron would however be the second daughter if Y would be created like: BtaCandidate* Y = comb.create(*lepton,*hadron);
	  int lund = abs( (int)fille->pdtEntry()->lundId() );
	  if( lund == PdtLund::pi_plus && d==0 ) PiList->append( new BtaCandidate(*Y) ); 
	  else if( lund == PdtLund::pi0 && d==0 ) Pi0List->append( new BtaCandidate(*Y) ); 
	  else if( lund == PdtLund::eta && fille->nDaughters()==2 && d==0 ) Eta2List->append( new BtaCandidate(*Y) ); 
	  else if( lund == PdtLund::eta && fille->nDaughters()==3 && d==0 ) Eta3List->append( new BtaCandidate(*Y) ); 
	  else if( lund == PdtLund::eta_prime && fille->nDaughters()==3 && d==0 ) 
	    {
	      HepAListIterator<BtaCandidate>iter_fille2 = fille->daughterIterator();
	      bool notFound=true;
	      while ( 0 != (fille2 =iter_fille2()) && notFound )   
		{
		  if(fille2->pdtEntry()->lundId()==PdtLund::eta && fille2->nDaughters()==2) 
		    { EtapE2PPList->append( new BtaCandidate(*Y) ); notFound=false; }
		  else if(fille2->pdtEntry()->lundId()==PdtLund::eta && fille2->nDaughters()==3) 
		    { EtapE3PPList->append( new BtaCandidate(*Y) ); notFound=false; }
		}
	    }
	  else if( lund == PdtLund::eta_prime && fille->nDaughters()==2 && d==0 ) EtapRGList->append( new BtaCandidate(*Y) ); 
	  else if( lund == PdtLund::rho_plus && d==0 ) RhoCList->append( new BtaCandidate(*Y) ); 
	  else if( lund == PdtLund::rho0 && d==0 ) Rho0List->append( new BtaCandidate(*Y) ); 
	  else if( lund == PdtLund::omega && d==0 ) OmegaList->append( new BtaCandidate(*Y) ); 
	  else if( lund == PdtLund::gamma && d==0 ) GammaList->append( new BtaCandidate(*Y) ); 
	  else if(d==0) cout<<"Problem in XSLReader!!! (unknown hadron)"<<endl;

	  if(d>0) cout<<"Problem in XSLReader!!! (d>0)"<<endl;
	  d+=1;
	}
    }//while ( 0 != ( candid = iter_YList()) ) 


  //Fill the ntuple!
  FillNtupleWithEventVariables( _MainNtuple ); 
  FillNtupleWithBumpListBlock( anEvent, _MainNtuple );
  FillNtupleWithIfrListBlock( anEvent, _MainNtuple );
  FillNtupleWithTrkListBlock( anEvent, _MainNtuple );
  if( _lookMCTruth.value()==true ) FillNtupleWithMCTruthBlockInfo(anEvent, _MainNtuple);

  //FillNtupleWithMCTruthEventVariables has to be executed AFTER Trk/Bump list blocks because it needs _nBadRecoed and _MatchedTruthList
  //but BEFORE CoupleVariables which needs the XSLMCTruthAnalyzer
  XSLMCTruthAnalyzer* EvtTruth;
  if( _lookMCTruth.value()==true ){ EvtTruth = FillNtupleWithMCTruthEventVariables( anEvent, _MainNtuple ); }
  else EvtTruth = new XSLMCTruthAnalyzer();



  if(_analyzePilnu.value()==true) FillNtupleWithCoupleVariables( "pilnu", PiList, _MainNtuple, anEvent, EvtTruth );
  if(_analyzePi0lnu.value()==true) FillNtupleWithCoupleVariables( "pi0lnu", Pi0List, _MainNtuple, anEvent, EvtTruth );
  if(_analyzeEtalnu.value()==true) FillNtupleWithCoupleVariables( "eta2lnu", Eta2List, _MainNtuple, anEvent, EvtTruth );
  if(_analyzeEtalnu.value()==true) FillNtupleWithCoupleVariables( "eta3lnu", Eta3List, _MainNtuple, anEvent, EvtTruth );
  if(_analyzeEtaplnu.value()==true) FillNtupleWithCoupleVariables( "etaplnuE2PP", EtapE2PPList, _MainNtuple, anEvent, EvtTruth );
  if(_analyzeEtaplnu.value()==true) FillNtupleWithCoupleVariables( "etaplnuE3PP", EtapE3PPList, _MainNtuple, anEvent, EvtTruth );
  //if(_analyzeEtaplnu.value()==true) FillNtupleWithCoupleVariables( "etaplnuRG", EtapRGList, _MainNtuple, anEvent, EvtTruth ); 
  if(_analyzeRhoClnu.value()==true) FillNtupleWithCoupleVariables( "rhoClnu", RhoCList, _MainNtuple, anEvent, EvtTruth );
  if(_analyzeRho0lnu.value()==true) FillNtupleWithCoupleVariables( "rho0lnu", Rho0List, _MainNtuple, anEvent, EvtTruth );
  if(_analyzeOmegalnu.value()==true) FillNtupleWithCoupleVariables( "omegalnu", OmegaList, _MainNtuple, anEvent, EvtTruth );
  if(_analyzeGammalnu.value()==true) _MainNtuple->column("Gammalnu_Nb",GammaList->length());

  _MainNtuple->dumpData();

  delete EvtTruth;
  return;
}


void
XSLReader::FillNtupleWithCoupleVariables( string mode,  HepAList<BtaCandidate>* List , HepTuple* ntuple, 
					  AbsEvent* anEvent, XSLMCTruthAnalyzer* EvtTruth )
{
  
  //Declare couple block vectors here:
  //MC
  HTValOrderedVector<int> oldSigMC;
  HTValOrderedVector<float> SigMC;

  //Kin
  HTValOrderedVector<float> q2Ups,thLUps,thVUps,chiUps,q2Ave,thLAve,thVAve,chiAve,q2Neutrino;

  //Kin with TreeFitter
  HTValOrderedVector<float> q2UpsYTRConstr,q2AveYTRConstr,thVAveYTRConstr,thLAveYTRConstr,chiAveYTRConstr;

  //Lepton
  HTValOrderedVector<int> indexLep,BremGamInd; 
  HTValOrderedVector<float> pLepLab,thLepLab,phiLepLab,pLepUps,thLepUps,phiLepUps;

  //Xu
  HTValOrderedVector<float> pXuLab,thXuLab,phiXuLab,pXuUps,thXuUps,phiXuUps;
  HTValOrderedVector<float> pTXuLab,pTXuUps,mXu;
  HTValOrderedVector<int> NDaughters;

  //Y
  HTValOrderedVector<float> cosBY,opAngL0,helAngL2,helAngL3,mXuLep2,mm12,mm23,probChi2YTree,probChi2YFast,probChi2TreeToCutOn,cosBYFit;
  HTValOrderedVector<float> posXY,posZ,posTh;

  //Y+Pmiss
  HTValOrderedVector<float> delE,delE2,mES,mESfitted,delThNu,delThNuAve,delThNuAveFit;
  HTValOrderedVector<float> delEFitted,cosBXu,fittedcosBXu;

  //continuum rejection
  HTValOrderedVector<float> L0,L1,L2,L3,cosTT,cosYZ,cosTZ,thrustROE;  

  //Pi mode
  HTValOrderedVector<int> indexPi;  

  //Family of the Composite Xu's (all other modes)
  //fille1: Gam1/pi+ , fille2: Gam2/pi-, X0: pi0,eta,rho0
  HTValOrderedVector<int>  f1Ind,f2Ind,pf1Ind,pf2Ind,ppfGam1Ind,ppfGam2Ind; 

  HTValOrderedVector<float>  p_fX0Lab, pT_fX0Lab,th_fX0Lab,phi_fX0Lab,m_fX0;
  HTValOrderedVector<float>  p_fX0Ups, pT_fX0Ups,th_fX0Ups,phi_fX0Ups;

  HTValOrderedVector<float>  p_pfPi0Lab, pT_pfPi0Lab,th_pfPi0Lab,phi_pfPi0Lab,m_pfPi0;
  HTValOrderedVector<float>  p_pfPi0Ups, pT_pfPi0Ups,th_pfPi0Ups,phi_pfPi0Ups;


  ///////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////

  BtaBooster* CMS = new BtaBooster(_UpsLab); //deleted at the end of the method


  //Starting the Y couples loop
  BtaCandidate* Y;  
  HepAListIterator<BtaCandidate> iter_List(*List);  
  while ( 0 != ( Y = iter_List()) ) 
    {            
      XSLRecoYAnalyzer* recoY = new XSLRecoYAnalyzer();
      recoY->Init(Y);
      if(mode!=recoY->mode()) { 
	cout<<"Problem with mode in FillNtupleWithCoupleVariables!!"<<endl;
	cout<<"mode: "<<mode<<"        recoY->mode(): "<<recoY->mode()<<endl;
      }
      //else cout<<"recoMode: "<<recoY->mode()<<endl;

      //Lep
      BtaCandidate LepLab = recoY->LepLab();
      BtaCandidate LepUps = recoY->LepUps();

      pLepLab.append((float)LepLab.p());
      thLepLab.append((float)LepLab.p3().theta());
      phiLepLab.append((float)LepLab.p3().phi());
      pLepUps.append((float)LepUps.p());
      thLepUps.append((float)LepUps.p3().theta());
      phiLepUps.append((float)LepUps.p3().phi());

      int LepBremGamInd=GetBremPhotonIndex(LepLab);
      BremGamInd.append(LepBremGamInd);
      int lepind=GetTrkIndex(LepLab);
      indexLep.append(lepind);

      //Xu
      BtaCandidate XuLab = recoY->XuLab();
      BtaCandidate XuUps = recoY->XuUps();

      int Bchrg = (int)(XuLab.charge() + LepLab.charge());

      pXuLab.append((float)XuLab.p());
      //pTXuLab.append((float)XuLab.pt());
      thXuLab.append((float)XuLab.p3().theta());
      phiXuLab.append((float)XuLab.p3().phi());
      pXuUps.append((float)XuUps.p());
      //pTXuUps.append((float)XuUps.pt());
      thXuUps.append((float)XuUps.p3().theta());
      phiXuUps.append((float)XuUps.p3().phi());
      mXu.append((float)XuLab.mass());

      if(mode!="pilnu"){ NDaughters.append(XuLab.nDaughters()); }

      
      //Y
      cosBY.append((float)recoY->cosBY());
      opAngL0.append((float)recoY->openAngleL0());
      mXuLep2.append((float)recoY->mXuLep2());
      HepLorentzVector BUps = XuUps.p4()+LepUps.p4()+_PmissUps;
      double cBXu=BUps.vect().unit()*XuUps.p3().unit();
      cosBXu.append((float)cBXu);
      
      
      //Kin
      BtaCandidate VDLab = recoY->VDLab();
      HepLorentzVector BLab = 0.5*_UpsLab;
      XSLKin* UpsFrame;
      if(mode=="pilnu") UpsFrame = new XSLKin(BLab,LepLab.p4(),XuLab.p4());
      else UpsFrame = new XSLKin(BLab,LepLab.p4(),XuLab.p4(),VDLab.p4());

      q2Ups.append((float)UpsFrame->q2());
      thLUps.append((float)UpsFrame->theta_l());
      if(mode!="pilnu") 
	{
	  thVUps.append((float)UpsFrame->theta_v());
	  chiUps.append((float)UpsFrame->chi());
	}
      delete UpsFrame;
      
      //YAve
      YAveBooster YAve(LepLab.p4(),XuLab.p4(),Bchrg);
      YAve.ComputeKin(VDLab.p4());
      q2Ave.append((float) YAve.q2() );
      thLAve.append((float) YAve.thL() );
      if(mode!="pilnu") 
	{
	  thVAve.append((float) YAve.thV() );
	  chiAve.append((float) YAve.chi() );
	}
      delThNuAve.append((float) YAve.delThNu( _PmissLab ) );      


      HepLorentzVector NuLab(_PmissLab.x(),_PmissLab.y(),_PmissLab.z(),_PmissLab.rho()); //NuLab is used below for deltaE/mES
      HepLorentzVector WLab=NuLab+LepLab.p4();
      q2Neutrino.append((float)WLab.m2());

      //Y + Pmiss     
      delThNu.append((float)recoY->delThetaNu(_PmissUps));
      //delE.append((float)recoY->deltaE(_PmissUps));
      //delE2.append((float)recoY->deltaE2(_PmissUps));
      //mES.append((float)recoY->mES(_PmissUps));

      //Computing deltaE,mES as in BAD 53.      
      HepLorentzVector B1=NuLab+LepLab.p4()+XuLab.p4();
      HepLorentzVector B2=_PmissLab+LepLab.p4()+XuLab.p4();
      delE.append(DeltaEBAD53(B1));
      delE2.append(DeltaEBAD53(B2));
      mES.append(mESBAD53(B1.vect()));


      /////////////////////////////////////
      //Tree Fitter
      HepLorentzVector in = HepLorentzVector(666,666,666,666);
      BtaCandidate fittedXu(in);
      BtaCandidate fittedLep(in);
      BtaCandidate f1FitY(in),f2FitY(in),f0FitY(in);
      BtaCandidate fittedY = recoY->fittedY_Tree(_eventInfo,fittedXu,fittedLep,f1FitY,f2FitY,f0FitY);
	  
      int nDofYTree=0;
      double chi2YTree=0,probYTree=-666,cosBYF=-666,delThNuAveF=-666;
      double q2UpsYT=-666,q2AveYT=-666,thVAveYT=-666,thLAveYT=-666,chiAveYT=-666;
      float vtxPosXY=-666,vtxPosZ=-666,vtxPosTh=-666;

      if(fittedY.p()!=0) {
	if(fittedY.decayVtx()!=0)
	  {
	    nDofYTree=fittedY.decayVtx()->nDof();
	    chi2YTree=fittedY.decayVtx()->chiSquared();
	    probYTree=probab(nDofYTree,chi2YTree);	    

	    vtxPosTh=fittedY.decayVtx()->point().theta();
	    double tmpX=fittedY.decayVtx()->point().x();
	    double tmpY=fittedY.decayVtx()->point().y();
	    double tmpZ=fittedY.decayVtx()->point().z();

	    double evtX=0,evtY=0,evtZ=0;
	    if(_eventInfo->primaryVtx()!=0){
	      evtX= _eventInfo->primaryVtx()->point().x();
	      evtY= _eventInfo->primaryVtx()->point().y();
	      evtZ= _eventInfo->primaryVtx()->point().z();  }
	    else if(_eventInfo->beamSpot()!=0){
	      evtX=_eventInfo->beamSpot().x();
	      evtY=_eventInfo->beamSpot().y();
	      evtZ=_eventInfo->beamSpot().z();
	    }
	    else cout<<"XSLMainNtupleDumper: evt(X,Y,Z) not found and set to (0,0,0)"<<endl;
	      
	    double dX = evtX-tmpX;
	    double dY = evtY-tmpY;
	    vtxPosXY = sqrt(dX*dX + dY*dY);
	    vtxPosZ = evtZ-tmpZ;

	    //Kin after fit
	    BtaCandidate YTRVDLabFit = recoY->VDLabFromFittedDaughters(f1FitY,f2FitY,f0FitY);

	    XSLKin* UpsFrameFit;
	    if(mode=="pilnu") UpsFrameFit = new XSLKin(BLab,fittedLep.p4(),fittedXu.p4());
	    else UpsFrameFit = new XSLKin(BLab,fittedLep.p4(),fittedXu.p4(),YTRVDLabFit.p4());
	    q2UpsYT = UpsFrameFit->q2();
	    delete UpsFrameFit;

	    //YAve
	    YAveBooster YAveFit(fittedLep.p4(),fittedXu.p4(),Bchrg);
	    YAveFit.ComputeKin(YTRVDLabFit.p4());
	    q2AveYT=YAveFit.q2();
	    thLAveYT=YAveFit.thL();
	    thVAveYT=YAveFit.thV();
	    chiAveYT=YAveFit.chi();
	    cosBYF=YAveFit.cosBY();
	    delThNuAveF=YAveFit.delThNu( _PmissLab );
	  }
      }
      probChi2YTree.append((float)probYTree);
      cosBYFit.append((float)cosBYF);
      delThNuAveFit.append((float)delThNuAveF);      
      posXY.append(vtxPosXY);
      posZ.append(vtxPosZ);
      posTh.append(vtxPosTh);

      if(mode=="eta3lnu"||mode=="etaplnuE2PP"||mode=="etaplnuE3PP"||mode=="etaplnuRG"||mode=="omegalnu")
	{ probChi2TreeToCutOn.append( (float)recoY->ProbChi2TreeToCutOn(_eventInfo) ); }

      q2UpsYTRConstr.append((float)q2UpsYT);
      q2AveYTRConstr.append((float)q2AveYT);
      thLAveYTRConstr.append((float)thLAveYT);
      if(mode!="pilnu") 
	{
	  thVAveYTRConstr.append((float)thVAveYT);
	  chiAveYTRConstr.append((float)chiAveYT);
	}      
      
      ////////////////////////////
      //Fast Vtx
      if(mode!="pi0lnu"&&mode!="eta2lnu") probChi2YFast.append( (float)recoY->ProbChi2Y_Fast() );	  	
      
      
      //////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////
      //Continuum
      HepAList<BtaCandidate> YList; //The YList is filled with trks and photons only
      BtaCandidate* filleL1,*filleL2,*filleL3,*filleL4; //L1 is Xu, L2 is eta'->eta, L3 is eta->pi0, L4 is pi0->gamma      
      HepAListIterator<BtaCandidate> YL1 = Y->daughterIterator();
      while ( 0 != (filleL1 = YL1()) )   
	{
	  if(filleL1->nDaughters()==0) YList.append(filleL1);
	  else {
	    HepAListIterator<BtaCandidate> YL2 = filleL1->daughterIterator();
	    while ( 0 != (filleL2 = YL2()) )
	      {   
		if(filleL2->nDaughters()==0) YList.append(filleL2);
		else {
		  HepAListIterator<BtaCandidate> YL3 = filleL2->daughterIterator();
		  while ( 0 != (filleL3 = YL3()) )
		    {   
		      if(filleL3->nDaughters()==0) YList.append(filleL3);
		      else {
			HepAListIterator<BtaCandidate> YL4 = filleL3->daughterIterator();
			while ( 0 != (filleL4 = YL4()) )
			  {   
			    if(filleL4->nDaughters()==0) YList.append(filleL4);
			    else ErrMsg(fatal) << "Problem in XSLReader::FillNtupleWithCoupleVariables !!"<<endl;
			  }//while level4
		      }//else level4
		    }//while level 3
		}//else Level3
	      }//while Level2
	  }//else Level2
	}//Level1
      
      
      HepAList<BtaCandidate> ListAll; 
      ListAll.append( *BasicTrkList );
      ListAll.append( *BasicBumpList );
      HepAListIterator<BtaCandidate> iterAll( ListAll );

      HepAList<BtaCandidate> ROEList;  //Rest of event 
      BtaCandidate *CandAllL1,*CandAllL2, *CandY;
      while ( 0 != ( CandAllL1 = iterAll()) ) 
	{
	  if(CandAllL1->nDaughters()==0)
	    {
	      bool StillNotFound=true;
	      HepAListIterator<BtaCandidate> iterY( YList );
	      while ( 0 != ( CandY = iterY()) && StillNotFound ) 
		{
		  if(CandY->uid()==CandAllL1->uid()) StillNotFound=false;
		}//
	      if(StillNotFound) ROEList.append(CandAllL1);
	    }// if(CandAllL1->nDaughters()==0)
	  else
	    {
	      HepAListIterator<BtaCandidate> AllL2 = CandAllL1->daughterIterator();
	      while ( 0 != (CandAllL2 = AllL2()) )
		{   
		  bool StillNotFound=true;
		  HepAListIterator<BtaCandidate> iterY( YList );
		  while ( 0 != ( CandY = iterY()) && StillNotFound ) 
		    {
		      if(CandY->uid()==CandAllL2->uid()) StillNotFound=false;
		    }//
		  if(StillNotFound) ROEList.append(CandAllL2);
		}	      
	    }//else Level2 (this happens with Brem reco'ed electrons)
	}// while ( 0 != ( CandAllL1 = iterAll()) ) 



      ////////////////////////////////////////
      //  Thrust (input list should NOT be in the CMS frame)
      BtaThrust YTH(YList, *_eventInfo);		    
      Hep3Vector Yaxis = YTH.thrust_axis();
      
      BtaThrust RoeThrust(ROEList, *_eventInfo);		    
      thrustROE.append((float)RoeThrust.thrust());      
      Hep3Vector RoeAxis = RoeThrust.thrust_axis();

      cosTT.append((float) fabs( Yaxis.dot(RoeAxis) ) );


      //////////////////////////////
      //  ANGLES
      Hep3Vector Ydirection = recoY->YUps().p3().unit();
      Hep3Vector Z = _UpsLab.vect().unit();
      
      cosYZ.append((float)fabs( Z.dot(Ydirection) ));
      cosTZ.append((float)fabs( Z.dot(Yaxis) ));

      //////////////////////////
      //  The monomials! 
      double mom=0,cos=0,cos2=0;
      double Lzero=0,Lone=0,Ltwo=0,Lthree=0;

      BtaCandidate* roeCand;
      HepAListIterator<BtaCandidate> iterRoe(ROEList);
      while( 0 != (roeCand = iterRoe()) ){

	BtaCandidate boostedCand = CMS->boostTo(*roeCand);

	//monomials
	mom= boostedCand.p();
	cos = Yaxis.dot(boostedCand.p3().unit());
	cos = fabs(cos);
	cos2 = cos*cos;
	
	Lzero += mom; 
	Lone += mom*cos;
	Ltwo += mom*cos2;
	Lthree += mom*cos2*cos;
      }

      L0.append((float)Lzero);
      L1.append((float)Lone);
      L2.append((float)Ltwo);
      L3.append((float)Lthree);

      //end of continuum part
      ////////////////////////////////////////////////////////////////////////////////      
      ////////////////////////////////////////////////////////////////////////////////      
      
      
      if( mode == "pilnu" ) 
	{ 
	  indexPi.append(GetTrkIndex(XuLab)); 
	}

      else if(mode=="pi0lnu"||mode=="eta2lnu")
        {
	  int phoind2 = -666;
	  if(XuLab.nDaughters()!=0) //composite pi0/eta
	    {
	      //pi0/eta -> Gam1 Gam2
	      //Gam1
	      BtaCandidate Pho1=recoY->fille1Lab();
	      int phoind1=GetBumpIndex(Pho1);
	      f1Ind.append(phoind1);
	      
	      //Gam2
	      BtaCandidate Pho2=recoY->fille2Lab();
	      phoind2=GetBumpIndex(Pho2);
	      f2Ind.append(phoind2);
	    }
	  else //merged pi0/eta!
	    {
	      int phoind1=GetBumpIndex(XuLab);
	      f1Ind.append(phoind1);
	      f2Ind.append(phoind2);	      
	    }
        }//mode == "pi0lnu"
      
      else if( mode == "rho0lnu")
	{      
	  //rho0-> pi+ pi-

	  //pi+
	  BtaCandidate Pip=recoY->fille1Lab();
	  int pipind=GetTrkIndex(Pip);
	  f1Ind.append(pipind);

	  //pi-
	  BtaCandidate Pim=recoY->fille2Lab();
	  int pimind=GetTrkIndex(Pim);	  
	  f2Ind.append(pimind);          
	  
	}//mode == "rho0lnu"
      
      else if( mode == "rhoClnu")
	{   
	  helAngL2.append((float)recoY->HelicityAngleL2());

	  //rho+/- -> pi+/- pi0
	  BtaCandidate Pi;
	  if( XuLab.charge()>0)
	    {
	      Pi= recoY->fille1Lab() ;  
	    }else{
	      Pi= recoY->fille2Lab() ;  
	    }//(XuUps.charge()) == 1
	  
	  //pi +/-
	  int piind=GetTrkIndex(Pi);
	  f1Ind.append(piind);   //in this particular case, fPipInd stands for pi+ or pi- (but the name will be fPiInd in the ntuple)
	  
	  //pi0 (lab)
	  BtaCandidate Pi0Lab=recoY->filleX0Lab();
	  p_fX0Lab.append((float)Pi0Lab.p());
	  //pT_fX0Lab.append((float)Pi0Lab.pt());
	  th_fX0Lab.append((float)Pi0Lab.p3().theta());
	  phi_fX0Lab.append((float)Pi0Lab.p3().phi());
	  m_fX0.append((float)Pi0Lab.mass());           
	  
	  //pi0 (ups)
	  BtaCandidate Pi0Ups=recoY->filleX0Ups();
	  p_fX0Ups.append((float)Pi0Ups.p());
	  //pT_fX0Ups.append((float)Pi0Ups.pt());
	  th_fX0Ups.append((float)Pi0Ups.p3().theta());
	  phi_fX0Ups.append((float)Pi0Ups.p3().phi());

	  int phoind2 = -666;
	  if(Pi0Lab.nDaughters()!=0) //composite pi0
	    {
	      //pi0 -> Gam1 Gam2
	      //Gam1
	      BtaCandidate Pho1=recoY->pfille1Lab();
	      int phoind1=GetBumpIndex(Pho1);
	      pf1Ind.append(phoind1);
	      
	      //Gam2
	      BtaCandidate Pho2=recoY->pfille2Lab();
	      phoind2=GetBumpIndex(Pho2);
	      pf2Ind.append(phoind2);
	    }
	  else //merged pi0!
	    {
	      int phoind1=GetBumpIndex(Pi0Lab);
	      pf1Ind.append(phoind1);
	      pf2Ind.append(phoind2);	      
	    }
	  
	}// mode == "rhoClnu"
      
      else if(mode=="omegalnu"||mode=="eta3lnu"||mode=="etaplnuE2PP"||mode=="etaplnuE3PP")
	{   
	  helAngL2.append((float)recoY->HelicityAngleL2());
	  double m12=-666,m23=-666;
	  recoY->DalitzL1( m12, m23 );
	  mm12.append((float)m12);
	  mm23.append((float)m23);

	  //omega/eta/eta' -> pi+ pi- pi0/eta
	  //pi+
	  BtaCandidate Pip=recoY->fille1Lab();
	  int pipind=GetTrkIndex(Pip);
	  f1Ind.append(pipind);

	  //pi-
	  BtaCandidate Pim=recoY->fille2Lab();
	  int pimind=GetTrkIndex(Pim);	  
	  f2Ind.append(pimind);          

	  //pi0 (lab)
	  BtaCandidate X0Lab=recoY->filleX0Lab();
	  p_fX0Lab.append((float)X0Lab.p());
	  //pT_fX0Lab.append((float)X0Lab.pt());
	  th_fX0Lab.append((float)X0Lab.p3().theta());
	  phi_fX0Lab.append((float)X0Lab.p3().phi());
	  m_fX0.append((float)X0Lab.mass());           
	  
	  //pi0 (ups)
	  BtaCandidate X0Ups=recoY->filleX0Ups();
	  p_fX0Ups.append((float)X0Ups.p());
	  //pT_fX0Ups.append((float)X0Ups.pt());
	  th_fX0Ups.append((float)X0Ups.p3().theta());
	  phi_fX0Ups.append((float)X0Ups.p3().phi());

	  int phoind2 = -666;
	  if(X0Lab.nDaughters()==2) //composite pi0/eta ->mode=="omegalnu"||mode=="eta3lnu"||mode=="etaplnuE2PP"
	    {
	      //pi0/eta -> Gam1 Gam2
	      //Gam1
	      BtaCandidate Pho1=recoY->pfille1Lab();
	      int phoind1=GetBumpIndex(Pho1);
	      pf1Ind.append(phoind1);
	      
	      //Gam2
	      BtaCandidate Pho2=recoY->pfille2Lab();
	      phoind2=GetBumpIndex(Pho2);
	      pf2Ind.append(phoind2);
	    }
	  else if(X0Lab.nDaughters()==0) //merged pi0/eta! ->mode=="omegalnu"||mode=="eta3lnu"||mode=="etaplnuE2PP"
	    {
	      int phoind1=GetBumpIndex(X0Lab);
	      pf1Ind.append(phoind1);
	      pf2Ind.append(phoind2);	      
	    }	  
	  else if(X0Lab.nDaughters()==3)// mode == "etaplnuE3PP"
	    {   
	      /*
	      BtaCandidate fittedEta = recoY->fittedEtaL1(_eventInfo);
	      int nDofEta=0;
	      double chi2Eta=0,probEta=-666;
	      if(fittedEta!=0) {
		if(fittedEta.decayVtx()!=0)
		  {
		    nDofEta=fittedEta.decayVtx()->nDof();
		    chi2Eta=fittedEta.decayVtx()->chiSquared();
		    probEta=probab(nDofEta,chi2Eta);
		  }
	      }
	      probChi2EtaL1.append((float)probEta);
	      */

	      helAngL3.append((float)recoY->HelicityAngleL3());

	      //eta -> pi+ pi- pi0
	      //pi+
	      BtaCandidate Pip2=recoY->pfille1Lab();
	      int pipind2=GetTrkIndex(Pip2);
	      pf1Ind.append(pipind2);
	      
	      //pi-
	      BtaCandidate Pim2=recoY->pfille2Lab();
	      int pimind2=GetTrkIndex(Pim2);	  
	      pf2Ind.append(pimind2);          
	      
	      //pi0 (lab)
	      BtaCandidate Pi0Lab=recoY->pfillePi0Lab();
	      p_pfPi0Lab.append((float)Pi0Lab.p());
	      //pT_pfPi0Lab.append((float)Pi0Lab.pt());
	      th_pfPi0Lab.append((float)Pi0Lab.p3().theta());
	      phi_pfPi0Lab.append((float)Pi0Lab.p3().phi());
	      m_pfPi0.append((float)Pi0Lab.mass());           
	      
	      //pi0 (ups)
	      BtaCandidate Pi0Ups=recoY->pfillePi0Ups();
	      p_pfPi0Ups.append((float)Pi0Ups.p());
	      //pT_pfPi0Ups.append((float)Pi0Ups.pt());
	      th_pfPi0Ups.append((float)Pi0Ups.p3().theta());
	      phi_pfPi0Ups.append((float)Pi0Ups.p3().phi());

	      int phoind2 = -666;
	      if(Pi0Lab.nDaughters()!=0) //composite pi0
		{
		  //pi0 -> Gam1 Gam2
		  //Gam1
		  BtaCandidate Pho1=recoY->ppfilleGam1Lab();
		  int phoind1=GetBumpIndex(Pho1);
		  ppfGam1Ind.append(phoind1);
		  
		  //Gam2
		  BtaCandidate Pho2=recoY->ppfilleGam2Lab();
		  phoind2=GetBumpIndex(Pho2);
		  ppfGam2Ind.append(phoind2);
		}
	      else //merged pi0!
		{
		  int phoind1=GetBumpIndex(Pi0Lab);
		  ppfGam1Ind.append(phoind1);
		  ppfGam2Ind.append(phoind2);	      
		}
	      
	    }//else if(X0Lab.nDaughters()==3)// mode == "etaplnuE3PP"
	} //else if(mode=="omegalnu"||mode=="eta3lnu"||mode=="etaplnuE2PP"||mode=="etaplnuE3PP")
      else if( mode == "etaplnuRG")
	{   
	  helAngL2.append((float)recoY->HelicityAngleL2());

	  //etap ->rho0 Gam1
	  //Gam1
	  BtaCandidate Pho1=recoY->fille1Lab();
	  int phoind1=GetBumpIndex(Pho1);
	  f1Ind.append(phoind1);
	  
	  //rho0 (lab)
	  BtaCandidate Rho0Lab=recoY->filleX0Lab();
	  p_fX0Lab.append((float)Rho0Lab.p());
	  //pT_fX0Lab.append((float)Rho0Lab.pt());
	  th_fX0Lab.append((float)Rho0Lab.p3().theta());
	  phi_fX0Lab.append((float)Rho0Lab.p3().phi());
	  m_fX0.append((float)Rho0Lab.mass());           
	      
	  //rho0 (ups)
	  BtaCandidate Rho0Ups=recoY->filleX0Ups();
	  p_fX0Ups.append((float)Rho0Ups.p());
	  //pT_fX0Ups.append((float)Rho0Ups.pt());
	  th_fX0Ups.append((float)Rho0Ups.p3().theta());
	  phi_fX0Ups.append((float)Rho0Ups.p3().phi());
	  
	  //rho0-> pi+ pi-
	  //pi+
	  BtaCandidate Pip=recoY->pfille1Lab();
	  int pipind=GetTrkIndex(Pip);
	  pf1Ind.append(pipind);

	  //pi-
	  BtaCandidate Pim=recoY->pfille2Lab();
	  int pimind=GetTrkIndex(Pim);	  
	  pf2Ind.append(pimind);          

	}// mode == "etaplnuRG"        
      
      
      ///////////////////////////////////////////////////////////////////////////////////////      
      //MCTruth part
      
      if ( _lookMCTruth.value()==true ) { 
	oldSigMC.append( IsSignalOld(Y) );
	SigMC.append( IsSignal(recoY,EvtTruth) );
	
	int typeB1=EvtTruth->typeB1();
	int typeB2=EvtTruth->typeB2();
	if((typeB1>=1&&typeB1<=77)||(typeB2>=1&&typeB2<=77))
	  {
	    HepLorentzVector p4SigNuLab(0,0,0,0);
	    if(typeB2>=1&&typeB2<=77) p4SigNuLab = EvtTruth->Nu2_Lab().p4();  //we start with B2 because it's the good one for signal collections
	    else if(typeB1>=1&&typeB1<=77) p4SigNuLab = EvtTruth->Nu1_Lab().p4();

	    ////////////////////////
	    //boosting to the Ups(4s) frame
	    HepLorentzVector Ups;
	    if(_sigModeString=="continuum") Ups=_UpsLab;
	    else Ups = EvtTruth->Ups4S().p4();  //pure MC truth
	    XSLRotateAndBoost rAndB = XSLRotateAndBoost();	  
	    HepLorentzVector p4SigNuUps = rAndB.BoostToFrame( p4SigNuLab, Ups );
	  }
      }
      
      delete recoY;

    }//while ( 0 != ( candid = iter_List()) ) 


  ntuple->column(mode+"_Nb",List->length(),0 ,"cou"+mode,HTRange<int>(0,2000));

  if(mode!="pilnu"){ ntuple->column(mode+"_nDaugXu",NDaughters, mode+"_Nb",0, "cou"+mode); }
  ntuple->column(mode+"_pXuLab",pXuLab, mode+"_Nb",0, "cou"+mode);
  //ntuple->column(mode+"_pTXuLab",pTXuLab, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_thXuLab",thXuLab, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_phiXuLab",phiXuLab, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_pXuUps",pXuUps, mode+"_Nb",0, "cou"+mode);
  //ntuple->column(mode+"_pTXuUps",pTXuUps, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_thXuUps",thXuUps, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_phiXuUps",phiXuUps, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_mXu",mXu, mode+"_Nb",0, "cou"+mode);
  if( mode == "pilnu" ) ntuple->column(mode+"_XuTrkInd",indexPi, mode+"_Nb",0, "cou"+mode);

  ntuple->column(mode+"_pLepLab",pLepLab, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_thLepLab",thLepLab, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_phiLepLab",phiLepLab, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_pLepUps",pLepUps, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_thLepUps",thLepUps, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_phiLepUps",phiLepUps, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_indLep",indexLep, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_LepBremGamInd",BremGamInd, mode+"_Nb",0, "cou"+mode);
  
  ntuple->column(mode+"_cosBY",cosBY, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_opAngL0",opAngL0, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_mXuLep2",mXuLep2, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_deltaE",delE, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_deltaE2",delE2, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_mES",mES, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_cosBXu",cosBXu, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_delThNu",delThNu, mode+"_Nb",0, "cou"+mode);

  ntuple->column(mode+"_L0",L0, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_L1",L1, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_L2",L2, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_L3",L3, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_cosTT",cosTT, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_cosYZ",cosYZ, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_cosTZ",cosTZ, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_thrustROE",thrustROE, mode+"_Nb",0, "cou"+mode);

  ntuple->column(mode+"_q2Nu",q2Neutrino, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_q2Ups",q2Ups, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_q2UpsFit",q2UpsYTRConstr, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_q2Ave",q2Ave, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_q2AveFit",q2AveYTRConstr, mode+"_Nb",0, "cou"+mode);

  ntuple->column(mode+"_thLUps",thLUps, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_thLAve",thLAve, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_thLAveFit",thLAveYTRConstr, mode+"_Nb",0, "cou"+mode);

  if(mode!="pilnu")
    {     
      ntuple->column(mode+"_thVUps",thVUps, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_thVAve",thVAve, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_thVAveFit",thVAveYTRConstr, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_chiUps",chiUps, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_chiAve",chiAve, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_chiAveFit",chiAveYTRConstr, mode+"_Nb",0, "cou"+mode);
    }
  
  ntuple->column(mode+"_delThNuAve",delThNuAve, mode+"_Nb",0, "cou"+mode);  
  ntuple->column(mode+"_delThNuAveFit",delThNuAveFit, mode+"_Nb",0, "cou"+mode);  
  ntuple->column(mode+"_cosBYFit",cosBYFit, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_probChi2YTree",probChi2YTree, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_posXY",posXY, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_posZ",posZ, mode+"_Nb",0, "cou"+mode);
  ntuple->column(mode+"_posTh",posTh, mode+"_Nb",0, "cou"+mode);

  if(mode=="eta3lnu"||mode=="etaplnuE2PP"||mode=="etaplnuE3PP"||mode=="etaplnuRG"||mode=="omegalnu")
    {  ntuple->column(mode+"_probChi2TreeCut",probChi2TreeToCutOn, mode+"_Nb",0, "cou"+mode); }

  if(mode!="pi0lnu"&&mode!="eta2lnu")
    {
      ntuple->column(mode+"_probChi2YFast",probChi2YFast, mode+"_Nb",0, "cou"+mode);
    }  
  
  if(mode=="pi0lnu"||mode=="eta2lnu")
    {
      ntuple->column(mode+"_f1EmcInd",f1Ind, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_f2EmcInd",f2Ind, mode+"_Nb",0, "cou"+mode);
      
    }//mode == "pi0lnu" || "eta2lnu"
    
  if(mode=="rho0lnu")
    {
      ntuple->column(mode+"_f1TrkInd",f1Ind, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_f2TrkInd",f2Ind, mode+"_Nb",0, "cou"+mode);
      
    }//mode == "rho0lnu"

  if( mode == "rhoClnu")
    {   
      ntuple->column(mode+"_f1TrkInd",f1Ind, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_pLabfX0",p_fX0Lab, mode+"_Nb",0, "cou"+mode);
      //ntuple->column(mode+"_pTLabfX0",pT_fX0Lab, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_thLabfX0",th_fX0Lab, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_phiLabfX0",phi_fX0Lab, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_mfX0",m_fX0, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_pUpsfX0",p_fX0Ups, mode+"_Nb",0, "cou"+mode);
      //ntuple->column(mode+"_pTUpsfX0",pT_fX0Ups, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_thUpsfX0",th_fX0Ups, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_phiUpsfX0",phi_fX0Ups, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_pf1EmcInd",pf1Ind, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_pf2EmcInd",pf2Ind, mode+"_Nb",0, "cou"+mode);

      ntuple->column(mode+"_helAngL2",helAngL2, mode+"_Nb",0, "cou"+mode);
      
    }// mode == "rhoClnu"
  
  if(mode=="omegalnu"||mode=="eta3lnu"||mode=="etaplnuE2PP"||mode=="etaplnuE3PP")
    {   
      ntuple->column(mode+"_f1TrkInd",f1Ind, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_f2TrkInd",f2Ind, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_pLabfX0",p_fX0Lab, mode+"_Nb",0, "cou"+mode);
      //ntuple->column(mode+"_pTfX0Lab",pT_fX0Lab, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_thLabfX0",th_fX0Lab, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_phiLabfX0",phi_fX0Lab, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_mfX0",m_fX0, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_pUpsfX0",p_fX0Ups, mode+"_Nb",0, "cou"+mode);
      //ntuple->column(mode+"_pTUpsfX0",pT_fX0Ups, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_thUpsfX0",th_fX0Ups, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_phiUpsfX0",phi_fX0Ups, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_helAngL2",helAngL2, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_m12",mm12, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_m23",mm23, mode+"_Nb",0, "cou"+mode);

      if(mode=="etaplnuE2PP"||mode=="omegalnu"||mode=="eta3lnu")
	{
	  ntuple->column(mode+"_pf1EmcInd",pf1Ind, mode+"_Nb",0, "cou"+mode);
	  ntuple->column(mode+"_pf2EmcInd",pf2Ind, mode+"_Nb",0, "cou"+mode);
	}
      
      if( mode == "etaplnuE3PP")
	{   
	  ntuple->column(mode+"_pf1TrkInd",pf1Ind, mode+"_Nb",0, "cou"+mode);
	  ntuple->column(mode+"_pf2TrkInd",pf2Ind, mode+"_Nb",0, "cou"+mode);
	  ntuple->column(mode+"_pLabpfX0",p_pfPi0Lab, mode+"_Nb",0, "cou"+mode);
	  //ntuple->column(mode+"_pTLabpfX0",pT_pfPi0Lab, mode+"_Nb",0, "cou"+mode);
	  ntuple->column(mode+"_thLabpfX0",th_pfPi0Lab, mode+"_Nb",0, "cou"+mode);
	  ntuple->column(mode+"_phiLabpfX0",phi_pfPi0Lab, mode+"_Nb",0, "cou"+mode);
	  ntuple->column(mode+"_mpfX0",m_pfPi0, mode+"_Nb",0, "cou"+mode);
	  ntuple->column(mode+"_pUpspfX0",p_pfPi0Ups, mode+"_Nb",0, "cou"+mode);
	  //ntuple->column(mode+"_pTUpspfX0",pT_pfPi0Ups, mode+"_Nb",0, "cou"+mode);
	  ntuple->column(mode+"_thUpspfX0",th_pfPi0Ups, mode+"_Nb",0, "cou"+mode);
	  ntuple->column(mode+"_phiUpspfX0",phi_pfPi0Ups, mode+"_Nb",0, "cou"+mode);
	  ntuple->column(mode+"_ppf1EmcInd",ppfGam1Ind, mode+"_Nb",0, "cou"+mode);
	  ntuple->column(mode+"_ppf2EmcInd",ppfGam2Ind, mode+"_Nb",0, "cou"+mode);

	  ntuple->column(mode+"_helAngL3",helAngL3, mode+"_Nb",0, "cou"+mode);
	  //ntuple->column(mode+"_probChi2EtaL1",probChi2EtaL1, mode+"_Nb",0, "cou"+mode);
	  
	}// mode == "etaplnuE3PP"   
      
    }//(mode=="omegalnu"||mode=="eta3lnu"||mode=="etaplnuE2PP"||mode=="etaplnuE3PP")
     
  if( mode == "etaplnuRG")
    {   
      ntuple->column(mode+"_f1EmcInd",f1Ind, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_pLabfX0",p_fX0Lab, mode+"_Nb",0, "cou"+mode);
      //ntuple->column(mode+"_pTLabfX0",pT_fX0Lab, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_thLabfX0",th_fX0Lab, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_phiLabfX0",phi_fX0Lab, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_mfX0",m_fX0, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_pUpsfX0",p_fX0Ups, mode+"_Nb",0, "cou"+mode);
      //ntuple->column(mode+"_pTUpsfX0",pT_fX0Ups, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_thUpsfX0",th_fX0Ups, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_phiUpsfX0",phi_fX0Ups, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_pf1TrkInd",pf1Ind, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_pf2TrkInd",pf2Ind, mode+"_Nb",0, "cou"+mode);

      ntuple->column(mode+"_helAngL2",helAngL2, mode+"_Nb",0, "cou"+mode);
    }// mode == "etaplnuRG"        
  
  
  if ( _lookMCTruth.value()==true ) 
    {
      ntuple->column(mode+"_oldSigMC",oldSigMC, mode+"_Nb",0, "cou"+mode);
      ntuple->column(mode+"_SigMC",SigMC, mode+"_Nb",0, "cou"+mode);
    }

  delete CMS;
  
  return;
}



float
XSLReader::DeltaEBAD53(HepLorentzVector B){
  return (2.*_UpsLab.dot(B) - _s)/2./sqrt(_s);
}

float
XSLReader::mESBAD53(Hep3Vector B){
  //One difference with BAD 53, the value of mES is systematically shifted to max=5.29
  float prod = B.dot(_UpsLab.vect());
  float Ezero = _UpsLab.e();
  float arg = (_s/2 + prod)*(_s/2 + prod)/Ezero/Ezero - B.mag2();
  float mES = (arg>0) ? sqrt(arg) : -sqrt(-arg);
  
  float offset=_EbeamCM-5.29;
  return mES-offset;
}







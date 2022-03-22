//-----------------------------------------------------------------------------//
//                                                                             //
//  Analysis module for the B -> pi/pi0/eta/eta'/rho/rho0/omega l nu analyses  //
//                                                                             //
//     Sylvie Brunet  2002-2004 Universite de Montreal                         //
//     David Cote     2003-2004 Universite de Montreal                         //
//     Benoit Viaud        2004 Universite de Montreal                         //
//                                                                             //
//-----------------------------------------------------------------------------//
#include "BaBar/BaBar.hh"

#include "XslReader/XSLReader.hh"

#include "BetaMC/BtaMcAssocGHit.hh"
#include "HepTuple/HTValOrderedVector.h"
#include "HepTuple/Tuple.h"
#include "AbsEvent/getTmpAList.hh"
//#include "BetaCoreTools/BtaPrintTree.hh"

//XSL code
#include "XslTools/XSLRotateAndBoost.hh"
#include "XslTools/XSLMCTruthAnalyzer.hh"
#include "XslTools/XSLYAveBooster.hh"
#include "XslFFReweighting/XSLKin.hh"
#include "XslFFReweighting/XSLVectorISGW2.hh"
#include "XslFFReweighting/XSLPseudoScalarISGW2.hh"
#include "XslFFReweighting/XSLBall98.hh"
#include "XslFFReweighting/XSLBall01.hh"
#include "XslFFReweighting/XSLBall04_pilnu.hh"
#include "XslFFReweighting/XSLBall04_etalnu.hh"
#include "XslFFReweighting/XSLBall05.hh"
#include "XslFFReweighting/XSLFnal04.hh"
#include "XslFFReweighting/XSLHpqcd04.hh"

#include <cstdlib>
#include <cmath>
using std::cout;
using std::endl;
///////////////

XSLMCTruthAnalyzer*
XSLReader::FillNtupleWithMCTruthEventVariables( AbsEvent* anEvent, HepTuple* ntuple)
{

  XSLMCTruthAnalyzer* EvtTruth = new XSLMCTruthAnalyzer(*_MCTruthList);
  int typeB1=EvtTruth->typeB1();
  int typeB2=EvtTruth->typeB2();
  ntuple->column("Tru_EvtType",EvtTruth->eventType());
  ntuple->column("Tru_typeB1",typeB1 );
  ntuple->column("Tru_typeB2",typeB2 );

  //if we have the bad Xu decay mode, this is still signal!
  //(However, we want to preserve the information of the - sign in the ntuples)
  typeB1=abs(typeB1); 
  typeB2=abs(typeB2);


  //if sigModeInt==0, we look for any signal mode
  //if sigModeInt!=0, we look for this particular decay mode only
  int mode_e=_sigModeInt;
  int mode_mu=mode_e*10+mode_e;

  int NbSigDec=0;
  HepLorentzVector p4SigNu(0,0,0,0);
  bool B1IsSig=false,B2IsSig=false;
  if(typeB1==mode_e||typeB1==mode_mu||(_sigModeInt==0&&typeB1<=77)) 
    {
      NbSigDec+=1;
      p4SigNu=EvtTruth->Nu1_Lab().p4();
      B1IsSig=true;
    }
  if(typeB2==mode_e||typeB2==mode_mu||(_sigModeInt==0&&typeB2<=77)) 
    {
      NbSigDec+=1;
      p4SigNu=EvtTruth->Nu2_Lab().p4(); //This reset the value given if typeB1==sig, but it's OK (note that in signal coll the sig B is the 2nd one)
      B2IsSig=true;
    }

  ///////////////////////////////////////////
  // Xtra Missing stuff
  ///////////////////////////////////////////
  HepLorentzVector Neutrons = EvtTruth->p4Neutrons();
  HepLorentzVector KLongs = EvtTruth->p4KLongs();
  HepLorentzVector XtraNus = EvtTruth->p4Neutrinos();
  int nbNus=EvtTruth->nbNeutrinos();
  
  if(NbSigDec>0) { nbNus-=1; XtraNus = XtraNus-p4SigNu; }

  //KL0
  ntuple->column("Tru_nKL",EvtTruth->nbKLongs() );
  ntuple->column("Tru_eKL",(float)KLongs.e() );
  ntuple->column("Tru_pKL",(float)KLongs.rho() );
  ntuple->column("Tru_thKL",(float)KLongs.theta() );
  ntuple->column("Tru_phiKL",(float)KLongs.phi() );

  //neutrons
  ntuple->column("Tru_nNeutron",EvtTruth->nbNeutrons() );
  ntuple->column("Tru_eNeutron",(float)Neutrons.e() );
  ntuple->column("Tru_pNeutron",(float)Neutrons.rho() );
  ntuple->column("Tru_thNeutron",(float)Neutrons.theta() );
  ntuple->column("Tru_phiNeutron",(float)Neutrons.phi() );

  //Extra Neutrinos
  ntuple->column("Tru_nXnu", nbNus );
  ntuple->column("Tru_eXnu",(float)XtraNus.e() );
  ntuple->column("Tru_pXnu",(float)XtraNus.rho() );
  ntuple->column("Tru_thXnu",(float)XtraNus.theta() );
  ntuple->column("Tru_phiXnu",(float)XtraNus.phi() );

  //Extra tracks not in the MC tree
  //These are filled in FillNtupleWithBumpListBlock and  FillNtupleWithTrkListBlock, which has to be executed first!
  ntuple->column("Tru_nBadRecoed",_nBadRecoed);
  ntuple->column("Tru_eBadRecoed",(float)_BadRecoed.e());
  ntuple->column("Tru_pBadRecoed",(float)_BadRecoed.rho());
  ntuple->column("Tru_thBadRecoed",(float)_BadRecoed.theta());
  ntuple->column("Tru_phiBadRecoed",(float)_BadRecoed.phi());

  //Missing trk/neutrals (Klongs not included so far... (considered missing))
  int nMissingRecoed=0;
  HepLorentzVector MissingRecoed=HepLorentzVector(0,0,0,0);

  BtaCandidate *TruthCand,*TruthCand2;
  HepAListIterator<BtaCandidate> iterMC(*_MCTruthList);
  HepAListIterator<BtaCandidate> iterMC2(*_MatchedTruthList);
  while ( 0 != (TruthCand = iterMC()) )         
    {
      bool notFound=true;
      
      int lund = abs( (int)TruthCand->pdtEntry()->lundId() );
      int UID = (int)TruthCand->uid();
      int MotherLund = -666;
      BtaCandidate* moth=0;
      if(TruthCand->theMother()!=0) { moth=TruthCand->theMother(); MotherLund = abs( (int)TruthCand->theMother()->pdtEntry()->lundId() ); }

      if(TruthCand->nDaughters()==0 && !(lund==12||lund==14||lund==16||lund==130)
	 && !(MotherLund==PdtLund::pi_plus || MotherLund==PdtLund::K_plus || MotherLund==PdtLund::mu_minus || MotherLund==PdtLund::K_L0 ) )
	//Stable particles (proton, gamma, electron,...other?), but not neutrinos or Klongs! 
	{
	  iterMC2.rewind();
	  while ( 0 != (TruthCand2 = iterMC2()) && notFound ) 
	    { 
	      if((int)TruthCand2->uid()==UID) { notFound=false;  } 
	    }
	  
	  if(notFound)
	    { 
	      nMissingRecoed+=1; 
	      MissingRecoed+=TruthCand->p4(); 
	    }
	}

      if( notFound!=false && (lund==PdtLund::pi_plus || lund==PdtLund::K_plus || lund==PdtLund::mu_minus )
	       && !(MotherLund==PdtLund::pi_plus || MotherLund==PdtLund::K_plus || MotherLund==PdtLund::mu_minus || MotherLund==PdtLund::K_L0 ) )
	//these are not necessarly real stable particles from EvtGen's point of view but they are from the reconstruction one!
	{
	  bool notFound=true;
	  iterMC2.rewind();
	  while ( 0 != (TruthCand2 = iterMC2()) && notFound ) 
	    { 
	  if((int)TruthCand2->uid()==UID) { notFound=false; } 
	    }
	  
	  if(notFound)
	    { 
	      nMissingRecoed+=1; 
	      MissingRecoed+=TruthCand->p4();  
	    }
	}
      
    }// while ( 0 != (TruthCand = iterMC()) )  
  ntuple->column("Tru_nMissedRecoed",nMissingRecoed);
  ntuple->column("Tru_eMissedRecoed",(float)MissingRecoed.e());
  ntuple->column("Tru_pMissedRecoed",(float)MissingRecoed.rho());
  ntuple->column("Tru_thMissedRecoed",(float)MissingRecoed.theta());
  ntuple->column("Tru_phiMissedRecoed",(float)MissingRecoed.phi());
 

  ///////////////////////////////
  // Signal kin info
  //////////////////////////////

  //BtaPrintTree printTree;
  //cout << printTree.print(*sigB_Lab) << endl;

  //Truth block
  HTValOrderedVector<float> pSigBLab,mLep,pLepLab,thLepLab,pNuLab,pXuLab,mX;
  HTValOrderedVector<float> q2Ups,thLUps,thVUps,chiUps,pSigBUps,pLepUps,thLepUps,pNuUps,pXuUps;
  HTValOrderedVector<float> pSigB_B,thSigB_B,phiSigB_B,pLepB,thLepB,phiLepB;
  HTValOrderedVector<float> pNuB,thNuB,phiNuB,pXuB,thXuB,phiXuB,q2B,thLB,chiB,thVB;
  HTValOrderedVector<float> q2YAve,thLYAve,thVYAve,chiYAve,q2PDG,pVDaugB,thVDaugB,phiVDaugB;
  HTValOrderedVector<float> FLATQ2ToLCSR,FLATQ2ToISGW2;
  HTValOrderedVector<float> FLATQ2ToBall04,FLATQ2ToHPQCD04,FLATQ2ToFNAL04,FtoAlpha52,FtoAlpha61,FtoAlpha70;
  HTValOrderedVector<float> ISGW2ToLCSR,ISGW2ToBall04,ISGW2ToHPQCD04,ISGW2ToFNAL04,ItoAlpha52,ItoAlpha61,ItoAlpha70;
  HTValOrderedVector<int> lepLund,modInt,indB,indXu,indLep,indNu;

  int pos=1;
  while(pos<3 && (B1IsSig||B2IsSig))
    {
      //LAB FRAME
      BtaCandidate sigBLab;
      BtaCandidate LepLab;
      BtaCandidate NuLab;
      BtaCandidate XuLab;
      BtaCandidate VDaugLab;
      int type=0;
      string mode;
      bool doIt=false;

      if(B2IsSig&&pos==1) //Always B2 first for easier usage of signal collections
	{
	  sigBLab=EvtTruth->B2();
	  LepLab=EvtTruth->Lep2_Lab();
	  NuLab=EvtTruth->Nu2_Lab();
	  XuLab=EvtTruth->Xu2_Lab();
	  VDaugLab=EvtTruth->VDaug2_Lab();
	  type=typeB2;
	  mode=EvtTruth->modeB2();
	  doIt=true;
	}
      else if(B1IsSig&&pos==2)
	{
	  sigBLab=EvtTruth->B1();
	  LepLab=EvtTruth->Lep1_Lab();
	  NuLab=EvtTruth->Nu1_Lab();
	  XuLab=EvtTruth->Xu1_Lab();
	  VDaugLab=EvtTruth->VDaug1_Lab();
	  type=typeB1;
	  mode=EvtTruth->modeB1();
	  doIt=true;
	}

      if(doIt)
	{
	  if(mode=="pilnu") modInt.append(1);
	  else if(mode=="pi0lnu") modInt.append(2);
	  else if(mode=="eta2lnu") modInt.append(3);
	  else if(mode=="eta3lnu") modInt.append(4);
	  else if(mode=="etaplnuE2PP") modInt.append(5);
	  else if(mode=="etaplnuE3PP") modInt.append(6);
	  else if(mode=="etaplnuRG") modInt.append(7);
	  else if(mode=="rhoClnu") modInt.append(8);
	  else if(mode=="rho0lnu") modInt.append(9);
	  else if(mode=="omegalnu") modInt.append(10);
	  else if(mode=="pi0lnu_o") modInt.append(-2);
	  else if(mode=="etalnu_o") modInt.append(-3);
	  else if(mode=="etaplnu_o") modInt.append(-5);
	  else if(mode=="omegalnu_o") modInt.append(-10);
	  else modInt.append(-666);

	  indB.append( CandIndexInList( &sigBLab, _MCTruthList) );
	  indXu.append( CandIndexInList( &XuLab, _MCTruthList) );
	  indLep.append( CandIndexInList( &LepLab, _MCTruthList) );
	  indNu.append( CandIndexInList( &NuLab, _MCTruthList) );

	  //This trick is to avoid the effect of PHOTOS
	  //HepLorentzVector LepLab2 = sigBLab.p4()-XuLab.p4()-NuLab.p4();
	  
	  /*
	  pSigBLab.append( (float)sigBLab.p() );
	  //mLep.append( (float)LepLab.mass() );
	  pLepLab.append( (float)LepLab.p() );
	  thLepLab.append( (float)LepLab.p3().theta() );
	  pNuLab.append( (float)NuLab.p() );
	  pXuLab.append( (float)XuLab.p() );
	  mX.append( (float)XuLab.mass() );
	  */

	  
	  //BOOST
	  //Here's how we want to procede for the boost:
	  //First, we boost in the rotated Ups(4S) frame and coordinates using BtaBooster's rotateAndBoost
	  //Then, we boost tot the B frame while KEEPING the Ups(4S) coordiantes using standard boost
	  
	  //HepLorentzVector Ups = gblEnv->getPep()->pepBeams()->total4Momentum();  //reco style
	  //HepLorentzVector Ups(-0.1102553491300846,0.0,5.8772159072547439,12.101750056479895); //new reco style 
	                                                                            //(this way there is no problem for cosBY with OffPeak)
	  HepLorentzVector Ups = EvtTruth->Ups4S().p4();  //pure MC truth
	  XSLRotateAndBoost rAndB = XSLRotateAndBoost();
	  
	  ////////////////////////
	  //boosting to the Ups(4s) frame
	  HepLorentzVector sigBUps = rAndB.BoostToFrame( sigBLab.p4(), Ups );
	  HepLorentzVector LepUps = rAndB.BoostToFrame( LepLab.p4(), Ups );
	  HepLorentzVector NuUps = rAndB.BoostToFrame( NuLab.p4(), Ups );
	  HepLorentzVector XuUps = rAndB.BoostToFrame( XuLab.p4(), Ups );
	  HepLorentzVector VDaugUps = rAndB.BoostToFrame( VDaugLab.p4(), Ups );
	  
	  ////////////////////////////////////////
	  //Kin in the Ups4S frame:
	  HepLorentzVector UpsAsBLab = 0.5*Ups;
	  
	  XSLKin* UpsFrame;
	  if(!(type==1||type==11)) UpsFrame = new XSLKin(UpsAsBLab,LepLab.p4(),XuLab.p4(),VDaugLab.p4());
	  else  UpsFrame = new XSLKin(UpsAsBLab,LepLab.p4(),XuLab.p4());
	  
	  q2Ups.append( (float)UpsFrame->q2() );
	  thLUps.append( (float)UpsFrame->theta_l() );
	  if(_sigModeInt==1){
	    thVUps.append( (float)UpsFrame->theta_v() );
	    chiUps.append( (float)UpsFrame->chi() );
	  }
	  pSigBUps.append( (float)sigBUps.vect().mag() );
	  pLepUps.append( (float)LepUps.vect().mag() );
	  thLepUps.append( (float)LepUps.theta() );
	  pNuUps.append( (float)NuUps.vect().mag() );
	  pXuUps.append( (float)XuUps.vect().mag() );
	  
	  
	  //////////////////////////////
	  //B Frame
	  
	  //Boosting to B frame...
	  HepLorentzVector sigB_B = rAndB.BoostToFrame( XuUps, sigBUps );
	  HepLorentzVector LepB = rAndB.BoostToFrame( LepUps, sigBUps );
	  HepLorentzVector NuB = rAndB.BoostToFrame( NuUps, sigBUps );
	  HepLorentzVector XuB = rAndB.BoostToFrame( XuUps, sigBUps );
	  HepLorentzVector VDaugB = rAndB.BoostToFrame( VDaugUps, sigBUps );
	  
	  
	  /////////////////////////////////////////////////////
	  //Kin in the B frame:
	  XSLKin* BFrame;
	  if(!(type==1||type==11)) BFrame = new XSLKin(sigBLab.p4(),LepLab.p4(),XuLab.p4(),VDaugLab.p4());
	  else  BFrame = new XSLKin(sigBLab.p4(),LepLab.p4(),XuLab.p4());
	  
	  pSigB_B.append( (float)sigB_B.vect().mag() );
	  thSigB_B.append( (float)sigB_B.theta() );
	  phiSigB_B.append( (float)sigB_B.phi() );
	  pLepB.append( (float)LepB.vect().mag() );
	  thLepB.append( (float)LepB.theta() );
	  phiLepB.append( (float)LepB.phi() );
	  pNuB.append( (float)NuB.vect().mag() );
	  thNuB.append( (float)NuB.theta() );
	  phiNuB.append( (float)NuB.phi() );
	  pXuB.append( (float)XuB.vect().mag() );
	  thXuB.append( (float)XuB.theta() );
	  phiXuB.append( (float)XuB.phi() );
	  
	  float qq=(float)BFrame->q2();
	  q2B.append( qq );
	  thLB.append( (float)BFrame->theta_l() );
	  chiB.append( (float)BFrame->chi() );
	  thVB.append( (float)BFrame->theta_v() );	  
	  
	  
	  ////////////////////////////////////////
	  //Kin in the YAverage "Frame"
	  
	  //New object!
	  YAveBooster YAve(LepLab.p4(),XuLab.p4(),(int)sigBLab.charge());
	  YAve.ComputeKin();
	  q2YAve.append( (float)YAve.q2() );
	  thLYAve.append( (float)YAve.thL() );
	  if(_sigModeInt==1){
	    thVYAve.append( (float)YAve.thV() );
	    chiYAve.append( (float)YAve.chi() );  
	  }
	  
	  //q2 computed as (B-Xu)^2 in the B frame with PDG mean mass
	  double mB = sigBLab.p4().m();
	  double mXu = XuLab.pdtEntry()->mass();
	  double pXuB = XuB.vect().mag();
	  double eXuBframe = sqrt( mXu*mXu + pXuB*pXuB );
	  double Q2PDG =  mB*mB + mXu*mXu - 2*mB*eXuBframe;
	  q2PDG.append( (float)Q2PDG );
	  int LepLund = abs( (int)LepLab.pdtEntry()->lundId() );
	  lepLund.append( LepLund );
	  
	  
	  //(in case of pilnu, XuDaug=Xu)
	  pVDaugB.append( (float)VDaugB.vect().mag() );      
	  thVDaugB.append( (float)VDaugB.theta() );      
	  phiVDaugB.append( (float)VDaugB.phi() );      
	  
	  float w1=0,w2=0,w3=0,w4=0,w5=0,w61=0,w70=0,w52=0;
	  float w11=0,w33=0,w44=0,w55=0,ww61=0,ww70=0,ww52=0;
	  if(qq>0){
	    if((type>=1&&type<=4)||(type>=11&&type<=44))
	      {
		XSLEvtFFWeight* Ball01 = new XSLBall01(BFrame,mode); 
		w1=Ball01->FromFLATQ2ToThisModel();
		w11=Ball01->FromISGW2ToThisModel();
		delete Ball01;
		
		XSLEvtFFWeight* ISGW2 = new XSLPseudoScalarISGW2(BFrame,mode); 
		if(qq>0) w2=ISGW2->FromFLATQ2ToThisModel();
		delete ISGW2;

		if(type==1||type==2||type==11||type==22) 
		  { // B->pi/pi0nu		    
		    XSLEvtFFWeight* Ball04 = new XSLBall04_pilnu(BFrame,mode); 
		    w3=Ball04->FromFLATQ2ToThisModel();
		    w33=Ball04->FromISGW2ToThisModel();
		    delete Ball04;

		    XSLEvtFFWeight* HPQCD04 = new XSLHpqcd04(BFrame,mode); 
		    w4=HPQCD04->FromFLATQ2ToThisModel();
		    w44=HPQCD04->FromISGW2ToThisModel();
		    delete HPQCD04;

		    XSLEvtFFWeight* FNAL04 = new XSLFnal04(BFrame,mode); 
		    w5=FNAL04->FromFLATQ2ToThisModel();
		    w55=FNAL04->FromISGW2ToThisModel();
		    delete FNAL04;
		}
		else if(type==3||type==4||type==33||type==44) 
		  { // B->eta(')lnu
		    XSLEvtFFWeight* Ball04 = new XSLBall04_etalnu(BFrame,mode); 
		    w3=Ball04->FromFLATQ2ToThisModel();
		    w33=Ball04->FromISGW2ToThisModel();
		    delete Ball04;		    
		  }
	      }
	    else if((type>=5&&type<=7)||(type>=55&&type<=77))
	      {      
		XSLEvtFFWeight* Ball98 = new XSLBall98(BFrame,mode); 
		w1=Ball98->FromFLATQ2ToThisModel();
		w11=Ball98->FromISGW2ToThisModel();
		delete Ball98;	      

		XSLEvtFFWeight* ISGW2 = new XSLVectorISGW2(BFrame,mode); 
		w2=ISGW2->FromFLATQ2ToThisModel();
		delete ISGW2;	       		
		
		XSLEvtFFWeight* Ball05 = new XSLBall05(BFrame,mode); 
		w3=Ball05->FromFLATQ2ToThisModel();
		w33=Ball05->FromISGW2ToThisModel();
		delete Ball05;
	      }  
	    else cout<<"PROBLEM in Truth Block!!"<<endl;
	  }//if(qq>0)
	  FLATQ2ToLCSR.append(w1);      
	  FLATQ2ToISGW2.append(w2);      
	  FLATQ2ToBall04.append(w3);
	  FLATQ2ToHPQCD04.append(w4);
	  FLATQ2ToFNAL04.append(w5);
	  FtoAlpha52.append(w52);
	  FtoAlpha61.append(w61);
	  FtoAlpha70.append(w70);

	  ISGW2ToLCSR.append(w11);      
	  ISGW2ToBall04.append(w33);
	  ISGW2ToHPQCD04.append(w44);
	  ISGW2ToFNAL04.append(w55);
	  ItoAlpha52.append(ww52);
	  ItoAlpha61.append(ww61);
	  ItoAlpha70.append(ww70);
	  
	  //Clearing the memory...
	  delete BFrame;
	  delete UpsFrame;
	  //delete YYFrame;
	}//if(doIt)

      //increase the position in the loop!
      pos+=1;
    }

  ntuple->column("TruSig_Nb",NbSigDec,0 ,"SigDec",HTRange<int>(0,3));
 
  //Kin in lab frame -> les indices font la job...
  /*
  ntuple->column("TruSig_pSigBLab",pSigBLab,"TruSig_Nb",0 , "SigDec");
  //ntuple->column("TruSig_mLep",mLep,"TruSig_Nb",0 , "SigDec");
  ntuple->column("TruSig_pLepLab",pLepLab,"TruSig_Nb",0 , "SigDec");
  ntuple->column("TruSig_thLepLab",thLepLab,"TruSig_Nb",0 , "SigDec");
  ntuple->column("TruSig_pNuLab",pNuLab,"TruSig_Nb",0 , "SigDec");
  ntuple->column("TruSig_pXuLab",pXuLab,"TruSig_Nb",0 , "SigDec");
  ntuple->column("TruSig_mXu",mX,"TruSig_Nb",0 , "SigDec");
  */

  ntuple->column("TruSig_indB",indB,"TruSig_Nb",0 , "SigDec");
  ntuple->column("TruSig_indXu",indXu,"TruSig_Nb",0 , "SigDec");
  ntuple->column("TruSig_indLep",indLep,"TruSig_Nb",0 , "SigDec");
  ntuple->column("TruSig_indNu",indNu,"TruSig_Nb",0 , "SigDec");
  ntuple->column("TruSig_mode",modInt,"TruSig_Nb",0 , "SigDec");

  //Kin in Ups frame
  ntuple->column("TruSig_q2Ups",q2Ups,"TruSig_Nb",0 , "SigDec");
  ntuple->column("TruSig_thLUps",thLUps,"TruSig_Nb",0 , "SigDec");
  ntuple->column("TruSig_pSigBUps",pSigBUps,"TruSig_Nb",0 , "SigDec");
  ntuple->column("TruSig_pLepUps",pLepUps,"TruSig_Nb",0 , "SigDec");
  ntuple->column("TruSig_thLepUps",thLepUps,"TruSig_Nb",0 , "SigDec");
  ntuple->column("TruSig_pNuUps",pNuUps,"TruSig_Nb",0 , "SigDec");
  ntuple->column("TruSig_pXuUps",pXuUps,"TruSig_Nb",0 , "SigDec");
  if(_sigModeInt==1){
    ntuple->column("TruSig_thVUps",thVUps,"TruSig_Nb",0 , "SigDec");
    ntuple->column("TruSig_chiUps",chiUps,"TruSig_Nb",0 , "SigDec");
  }

  //Kin in B Frame
  ntuple->column("TruSig_pSigB_B",pSigB_B,"TruSig_Nb",0 , "SigDec");
  ntuple->column("TruSig_thSigB_B",thSigB_B,"TruSig_Nb",0 , "SigDec");
  ntuple->column("TruSig_phiSigB_B",phiSigB_B,"TruSig_Nb",0 , "SigDec");
  ntuple->column("TruSig_pLepB",pLepB,"TruSig_Nb",0 , "SigDec");
  ntuple->column("TruSig_thLepB",thLepB,"TruSig_Nb",0 , "SigDec");
  ntuple->column("TruSig_phiLepB",phiLepB,"TruSig_Nb",0 , "SigDec");
  ntuple->column("TruSig_pNuB",pNuB,"TruSig_Nb",0 , "SigDec");
  ntuple->column("TruSig_thNuB",thNuB,"TruSig_Nb",0 , "SigDec");
  ntuple->column("TruSig_phiNuB",phiNuB,"TruSig_Nb",0 , "SigDec");
  ntuple->column("TruSig_pXuB",pXuB,"TruSig_Nb",0 , "SigDec");
  ntuple->column("TruSig_thXuB",thXuB,"TruSig_Nb",0 , "SigDec");
  ntuple->column("TruSig_phiXuB",phiXuB,"TruSig_Nb",0 , "SigDec");
  ntuple->column("TruSig_q2B",q2B,"TruSig_Nb",0 , "SigDec");
  ntuple->column("TruSig_thLB",thLB,"TruSig_Nb",0 , "SigDec");
  ntuple->column("TruSig_chiB",chiB,"TruSig_Nb",0 , "SigDec");
  ntuple->column("TruSig_thVB",thVB,"TruSig_Nb",0 , "SigDec");


  ////////////////////////////////////////
  //Kin in the YAverage "Frame"
  ntuple->column("TruSig_q2YAve",q2YAve,"TruSig_Nb",0 , "SigDec");
  ntuple->column("TruSig_thLYAve",thLYAve,"TruSig_Nb",0 , "SigDec");
  if(_sigModeInt==1){
    ntuple->column("TruSig_thVYAve",thVYAve,"TruSig_Nb",0 , "SigDec");
    ntuple->column("TruSig_chiYAve",chiYAve,"TruSig_Nb",0 , "SigDec");  
  }

  //Special stuff...
  ntuple->column("TruSig_q2PDG",q2PDG,"TruSig_Nb",0 , "SigDec");
  ntuple->column("TruSig_lepLund",lepLund,"TruSig_Nb",0 , "SigDec");
  ntuple->column("TruSig_pVDaugB",pVDaugB,"TruSig_Nb",0 , "SigDec");      
  ntuple->column("TruSig_thVDaugB",thVDaugB,"TruSig_Nb",0 , "SigDec");      
  ntuple->column("TruSig_phiVDaugB",phiVDaugB,"TruSig_Nb",0 , "SigDec");      


  //FF reweighting
  ntuple->column("TruSig_FLATQ2ToLCSR",FLATQ2ToLCSR,"TruSig_Nb",0 , "SigDec");      
  ntuple->column("TruSig_FLATQ2ToISGW2",FLATQ2ToISGW2,"TruSig_Nb",0 , "SigDec");      
  ntuple->column("TruSig_FLATQ2ToBall04",FLATQ2ToBall04,"TruSig_Nb",0 , "SigDec");      
  ntuple->column("TruSig_FLATQ2ToHPQCD04",FLATQ2ToHPQCD04,"TruSig_Nb",0 , "SigDec");      
  ntuple->column("TruSig_FLATQ2ToFNAL04",FLATQ2ToFNAL04,"TruSig_Nb",0 , "SigDec");      
  ntuple->column("TruSig_FLATQ2ToAlpha52",FtoAlpha52,"TruSig_Nb",0 , "SigDec");      
  ntuple->column("TruSig_FLATQ2ToAlpha61",FtoAlpha61,"TruSig_Nb",0 , "SigDec");      
  ntuple->column("TruSig_FLATQ2ToAlpha70",FtoAlpha70,"TruSig_Nb",0 , "SigDec");      

  ntuple->column("TruSig_ISGW2ToLCSR",ISGW2ToLCSR,"TruSig_Nb",0 , "SigDec");      
  ntuple->column("TruSig_ISGW2ToBall04",ISGW2ToBall04,"TruSig_Nb",0 , "SigDec");      
  ntuple->column("TruSig_ISGW2ToHPQCD04",ISGW2ToHPQCD04,"TruSig_Nb",0 , "SigDec");      
  ntuple->column("TruSig_ISGW2ToFNAL04",ISGW2ToFNAL04,"TruSig_Nb",0 , "SigDec");      
  ntuple->column("TruSig_ISGW2ToAlpha52",ItoAlpha52,"TruSig_Nb",0 , "SigDec");      
  ntuple->column("TruSig_ISGW2ToAlpha61",ItoAlpha61,"TruSig_Nb",0 , "SigDec");      
  ntuple->column("TruSig_ISGW2ToAlpha70",ItoAlpha70,"TruSig_Nb",0 , "SigDec");      

  return EvtTruth;
}

void
XSLReader::FillNtupleWithMCTruthBlockInfo(AbsEvent* anEvent, HepTuple* ntuple)
{
  HTValOrderedVector<float> pLab,thLab,phiLab,mass;
  HTValOrderedVector<int> lund,motherIndex;//,nDaug;

  BtaCandidate* cand;
  HepAListIterator<BtaCandidate> iter_list(*_MCTruthList);    
  while ( 0 != ( cand = iter_list()) ) 
    {
      pLab.append((float)cand->p() );
      thLab.append((float)cand->p3().theta() );
      phiLab.append((float)cand->p3().phi() );
      mass.append((float)cand->mass() );
      lund.append((int)cand->pdtEntry()->lundId());
      //nDaug.append((int)cand->nDaughters());
      int ind=-666;
      if(cand->theMother()!=0) ind = CandIndexInList(cand->theMother(), _MCTruthList);
      motherIndex.append(ind);
    }
      
  ntuple->column("TruBlk_Nb",_MCTruthList->length(),0 ,"Bla",HTRange<int>(0,100));
  ntuple->column("TruBlk_pLab",pLab,"TruBlk_Nb",0 , "Bla");
  ntuple->column("TruBlk_thLab",thLab,"TruBlk_Nb",0 , "Bla");
  ntuple->column("TruBlk_phiLab",phiLab,"TruBlk_Nb",0 , "Bla");
  ntuple->column("TruBlk_m",mass,"TruBlk_Nb",0 , "Bla");
  ntuple->column("TruBlk_lund",lund,"TruBlk_Nb",0 , "Bla");
  //ntuple->column("TruBlk_nDaug",nDaug,"TruBlk_Nb",0 , "Bla");
  ntuple->column("TruBlk_motherInd",motherIndex,"TruBlk_Nb",0 , "Bla");

  return;
}



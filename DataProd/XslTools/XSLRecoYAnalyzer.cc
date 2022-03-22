//*****************************************************************
//   
//   Creation: David Cote, Universite de Montreal, 02/16/04
// 
#include "BaBar/BaBar.hh"

#include "XslTools/XSLRecoYAnalyzer.hh"

#include "AbsEnv/AbsEnv.hh"
#include "PepEnv/PepEnv.hh"
#include "PepEnvData/PepBeams.hh"
#include "BetaCoreTools/BtaBooster.hh"
#include "PDT/Pdt.hh"
#include "CLHEP/Alist/ConstAList.h"
#include "CLHEP/Alist/ConstAIterator.h"

#include "ProbTools/probab.hh"

#include "VtxFitter/VtxFitterOper.hh" 
#include "Beta/BtaAbsVertex.hh"
#include "FastVtx/BtaOpFastVtx.hh"
#include "FastVtx/BtaOpFastVtxV0.hh"
#include "CompositionUtils/BtaV0Finder.hh"

#include "Beta/EventInfo.hh"
#include "BetaCoreTools/BtaOpMakeTree.hh"
using std::cout;
using std::endl;
using std::string;

XSLRecoYAnalyzer::XSLRecoYAnalyzer()
{}

XSLRecoYAnalyzer::~XSLRecoYAnalyzer()
{}


void
XSLRecoYAnalyzer::Init(BtaCandidate* YLab)
{
  //////////////
  //Initializing
  _FamilyFiguredOut=false;
  HepLorentzVector in = HepLorentzVector(0,0,0,0);

  _VDLab = BtaCandidate(in);
  _VDUps = BtaCandidate(in);
  _LepLab = BtaCandidate(in);
  _LepUps = BtaCandidate(in);
  _XuLab = BtaCandidate(in);
  _XuUps = BtaCandidate(in);
  _YUps = BtaCandidate(in);
  _YLab = BtaCandidate(in);

  _fille1Lab = BtaCandidate(in);
  _fille2Lab = BtaCandidate(in);
  _filleX0Lab = BtaCandidate(in);
  _filleX0Ups = BtaCandidate(in);

  _filleX0LabForFit = BtaCandidate(in);

  _pfille1Lab = BtaCandidate(in);
  _pfille2Lab = BtaCandidate(in);
  _pfillePi0Lab = BtaCandidate(in);
  _pfillePi0Ups = BtaCandidate(in);

  _ppfilleGam1Lab = BtaCandidate(in);
  _ppfilleGam2Lab = BtaCandidate(in);

  //////////
  //Building stuff

  _YLab = BtaCandidate(*YLab);
  _LepLab = GetLep(_YLab);
  _XuLab = GetXu(_YLab);

  _Ups = gblEnv->getPep()->pepBeams()->total4Momentum();

  BtaBooster* CMS = new BtaBooster( _Ups );
  _YUps = CMS->boostTo( _YLab );
  delete CMS;

  _LepUps = GetLep(_YUps);
  _XuUps = GetXu(_YUps);

  //B
  _mB=0;  
  if(_XuLab.charge()==0) _mB=Pdt::lookup(PdtLund::B_plus)->mass();
  else _mB=Pdt::lookup(PdtLund::B0)->mass();  

  //E Beam
  _eBeamCM = sqrt(_Ups.mag2())/2.; 

  //correction factor for OffPeak
  _f=1.0;
  if( _eBeamCM<_mB ) {
    // 5.2895 is the mean run1-4 CM beam energy
    // see: http://babar-hn.slac.stanford.edu:5090/HyperNews/get/semi_lept_decays/177/2/1.html
    _f=5.2895/_eBeamCM;
  }


  //mode
  FigureOutXuFamily();

  int lund = abs( (int)_XuLab.pdtEntry()->lundId() );
  int nD = _XuLab.nDaughters();
  if(lund==PdtLund::pi_plus) _mode="pilnu";
  else if(lund==PdtLund::pi0 ) _mode="pi0lnu";
  else if(lund==PdtLund::eta && nD==2) _mode="eta2lnu";
  else if(lund==PdtLund::eta && nD==3) _mode="eta3lnu";
  else if(lund==PdtLund::eta_prime && nD==2 ) _mode="etaplnuRG"; 
  else if(lund==PdtLund::eta_prime && nD==3 )
    {
      int nDEta = _filleX0Lab.nDaughters();
      if(nDEta==2) _mode = "etaplnuE2PP"; 
      else if(nDEta==3) _mode = "etaplnuE3PP"; 
    }
  else if(lund==PdtLund::rho_plus) _mode="rhoClnu";
  else if(lund==PdtLund::rho0) _mode="rho0lnu";
  else if(lund==PdtLund::omega && nD==3 ) _mode="omegalnu";

  //The random pi0/eta daughter...
  /*
  if(_mode=="pi0lnu"||_mode=="eta2lnu")
    {
      unsigned int XuUid = (unsigned int) _XuLab.uid();
      unsigned int seed = XuUid*lowerID;
      srand( seed );
      float dice = rand()%2; //dice=0 or 1
      _DaughterToPick = ( dice<0.5) ? HigherEnergy:LowerEnergy;
    }
  else _DaughterToPick = NotApplicable;
  */

  return;
}


double 
XSLRecoYAnalyzer::openAngleL0()
{
  //This means Xu vs Lep opening angle 
  double ang = _XuUps.p3().angle( _LepUps.p3() );
  return ang;
}

double 
XSLRecoYAnalyzer::HelicityAngleL1()
{
  //Note: this is theta_V!! You might not want to use it...

  double hel=-666;
  if(_XuLab.nDaughters()<1) return hel;  

  //B->Xu l nu, Xu -> Daug X 
  //Starting from B frame, hel = angle between Daug in Xu frame and Xu in B frame

  HepLorentzVector Moth=_XuUps.p4();
  HepLorentzVector Daug = VDUps().p4();
  Daug.boost(-Moth.boostVector());
  
  hel = Daug.vect().angle(Moth.vect());
  return hel;
}

double 
XSLRecoYAnalyzer::HelicityAngleL2()
{
  double hel=-666;
  if(_XuLab.nDaughters()<1) return hel;  


  //Xu -> CompMoth X, CompMoth -> Daug Y 
  //Starting from Xu frame, hel = angle between Daug in CompMoth frame and CompMoth in Xu frame
  HepLorentzVector Moth(0), Daug(0);
  int nD= _filleX0Lab.nDaughters();
  if(!(nD==2||nD==3)) return hel;

  if( _mode=="omegalnu"||_mode=="eta3lnu"||_mode=="etaplnuE2PP"||_mode=="etaplnuE3PP"||_mode=="etaplnuRG"||_mode=="rhoClnu" ) 
    {
      Moth = _filleX0Lab.p4();
      if(nD==3) Daug = _pfillePi0Lab.p4();
      else if(nD==2) 
	{
	  bool pf1Higher=(_pfille1Lab.p3().phi()>_pfille2Lab.p3().phi()) ? true:false;	  
	  if(pf1Higher) Daug = _pfille1Lab.p4();
	  else Daug = _pfille2Lab.p4();
	}
      else { cout<<"Problem in XSLRecoYAnalyzer::HelicityAngleUpsL2()!!!"<<endl; return hel; }
    }
  else return hel;
  
  Moth.boost(-_XuLab.p4().boostVector());
  Daug.boost(-_XuLab.p4().boostVector());
  Daug.boost(-Moth.boostVector());
  
  hel = Daug.vect().angle(Moth.vect());
  return hel;
}


double 
XSLRecoYAnalyzer::HelicityAngleL3()
{
  double hel=-666;
  int nD= _pfillePi0Lab.nDaughters();
  if(_mode!="etaplnuE3PP" || nD!=2 ) return hel;

  //eta -> pi0 X, pi0 -> Gam Y 
  //Starting from eta frame, hel = angle between Gam in pi0 frame and pi0 in eta frame
  if(_pfillePi0Lab.nDaughters()==0) return hel; 
  HepLorentzVector Moth = _pfillePi0Lab.p4();
  HepLorentzVector Daug(0);
  	  
  bool ppf1Higher=(_ppfilleGam1Lab.p3().phi()>_ppfilleGam2Lab.p3().phi()) ? true:false;
  if(ppf1Higher) Daug = _ppfilleGam1Lab.p4();
  else Daug = _ppfilleGam2Lab.p4();
  
  Moth.boost(-_filleX0Lab.p4().boostVector());
  Daug.boost(-_filleX0Lab.p4().boostVector());
  Daug.boost(-Moth.boostVector());
  
  hel = Daug.vect().angle(Moth.vect());
  return hel;
}

double
XSLRecoYAnalyzer::mXuLep2()
{
  HepLorentzVector XuLep = _XuLab.p4()+_LepLab.p4();
  return XuLep.m2();
}


void 
XSLRecoYAnalyzer::DalitzL1( double &m12, double &m23 )
{
  //For Dalitz plot of eta/eta'/omega -> 3 bodys 
  // ->> m12 and m23 are mass squared! 
  m12=-666;
  m23=-666;

  if(_mode=="eta3lnu"||_mode=="etaplnuE2PP"||_mode=="etaplnuE3PP"||_mode=="omegalnu") 
    {
      HepLorentzVector a = _fille1Lab.p4()+_fille2Lab.p4();
      m12=a.m2();
      HepLorentzVector b = _filleX0Lab.p4()+_fille2Lab.p4();
      m23=b.m2();
    }

  return;
}

double
XSLRecoYAnalyzer::cosBY()
{
  double eBeam = _f*_eBeamCM;
  double mY = _f*_YUps.mass();
  double pY = _f*_YUps.p();
  double eY=sqrt(mY*mY+pY*pY);

  double pB=pBCM();
  //the _f is for the offres case. 
  //(This is what's done in the skim, in XSLRecoYAnalyzer::cosBY, XSLRecoYAnalyzer::pBCM and YAveBooster::Init)

  double cosBY = (2*eY*eBeam-_mB*_mB-mY*mY)/(2*pB*pY);
  return cosBY;
}

double 
XSLRecoYAnalyzer::pBCM()
{
  //the _f is for the offres case. 
  //(This is what's done in the skim, in XSLRecoYAnalyzer::cosBY, XSLRecoYAnalyzer::pBCM and YAveBooster::Init)
  double pBCM=sqrt( _f*_f*_eBeamCM*_eBeamCM-_mB*_mB );
  return pBCM;
}

BtaCandidate
XSLRecoYAnalyzer::GetXu( BtaCandidate Y )
{
  
  BtaCandidate hadron;
  BtaCandidate* fille;
  HepAListIterator<BtaCandidate>iter_fille = Y.daughterIterator();
  int d=0;
  while ( 0 != (fille =iter_fille()) && d==0 )   
    {  
      //The first daughter is always the hadron by construction
      //WARNING! this is the case only because the Y are created like this: BtaCandidate* Y = comb.create(*hadron,*lepton);
      //The hadron would however be the second daughter if Y would be created like: BtaCandidate* Y = comb.create(*lepton,*hadron);
      if(d==0 ) hadron = BtaCandidate(*fille); 
      d+=1;
    }

  return hadron;
}

BtaCandidate
XSLRecoYAnalyzer::GetLep( BtaCandidate Y )
{
  
  BtaCandidate lepton;
  BtaCandidate* fille;
  HepAListIterator<BtaCandidate>iter_fille = Y.daughterIterator();
  int d=0;
  while ( 0 != (fille =iter_fille()) && d<=1 )   
    {  
      //The first daughter is always the hadron by construction
      //WARNING! this is the case only because the Y are created like this: BtaCandidate* Y = comb.create(*hadron,*lepton);
      //The hadron would however be the second daughter if Y would be created like: BtaCandidate* Y = comb.create(*lepton,*hadron);
      if(d==1) lepton = BtaCandidate(*fille); 
      d+=1;
    }

  return lepton;
}


void
XSLRecoYAnalyzer::GetVD()
{
  //pilnu or merged pi0/eta...
  if(_XuLab.nDaughters()<1) { _VDLab = BtaCandidate(_XuLab); return; }

  //From Leif Wilden (/u/ec/wilden/vub/RELEASE/XslUser/Vub.cc):
  //rho0 -> take pi+, rho+/- -> take pi0, 
  //omega -> 3pi: boost the pi+ and the pi- in the omega rest frame, then take the Xproduct
  //omega -> anything else: take the first pion, whatever it is 
  
  //Mainly for testing of the FF reweighting:
  //eta, eta',pi0 -> two body, take any daughter
  //eta, eta' -> three body, take the pi+ pi- Xproduct
   

  //For pi0 et eta2, we want randomly fille1 ou fille2 as VDaug... ('cause there's a bias in the reco... :-(  )
  //New suggestion (Dave Brown), we take the daughter with higher phi...

  if(_mode=="pi0lnu"||_mode=="eta2lnu") 
    {
      bool f1Higher=(_fille1Lab.p3().phi()>_fille2Lab.p3().phi()) ? true:false;
      if(f1Higher) _VDLab = BtaCandidate(_fille1Lab);
      else _VDLab = BtaCandidate(_fille2Lab);
    }
  else if(_mode=="etaplnuRG") _VDLab = BtaCandidate(_fille1Lab);
  else if(_mode=="rhoClnu") _VDLab = BtaCandidate(_filleX0Lab);
  else if(_mode=="rho0lnu") _VDLab = BtaCandidate(_fille1Lab);
  else if(_mode=="eta3lnu"||_mode=="etaplnuE2PP"||_mode=="etaplnuE3PP"||_mode=="omegalnu")
    {
      HepLorentzVector piPlusHad(_fille1Lab.p4());
      HepLorentzVector piMinusHad(_fille2Lab.p4());
      piPlusHad.boost(-_XuLab.p4().boostVector());
      piMinusHad.boost(-_XuLab.p4().boostVector());
      HepLorentzVector p4VDaugXuFrame = HepLorentzVector(piPlusHad.vect().cross(piMinusHad.vect()),0.);
      
      //boost back in the lab frame for uniformity with other decay modes...
      HepLorentzVector p4VDaugLabFrame(p4VDaugXuFrame);
      p4VDaugLabFrame.boost(_XuLab.p4().boostVector());
      
      //Create the BtaCandidate (with no type, but we don't care)
      _VDLab = BtaCandidate(p4VDaugLabFrame);
    }
  else 
    {       
      HepAListIterator<BtaCandidate> dauIter(_XuLab.daughterIterator()); 
      _VDLab = BtaCandidate(*dauIter()); 
      cout<<"Unexpected case in XSLRecoYAnalyzer::GetVD()"<<endl; 
    }


  return;
}


BtaCandidate
XSLRecoYAnalyzer::VDLab()
{
  if(_VDLab.p()==0) GetVD();
  return _VDLab;
}

BtaCandidate
XSLRecoYAnalyzer::VDUps()
{
  if(_VDLab.p()==0) GetVD();
  if(_VDUps.p()==0)
    {
      BtaBooster* CMS = new BtaBooster(_Ups);
      _VDUps = CMS->boostTo(_VDLab);
      delete CMS;
    }

  return _VDUps;
}

void
XSLRecoYAnalyzer::FigureOutXuFamily()
{
  if(_FamilyFiguredOut) return;  //Only need to do this once! ;-)
  _FamilyFiguredOut=true;

  BtaCandidate* Cand;

  //Filles
  //Note: no daughter is set in the case of a merged pi0!
  HepAListIterator<BtaCandidate> Xu_Daug = _XuLab.daughterIterator();
  while ( 0 != (Cand = Xu_Daug()) )   
    {
      int lund = Cand->pdtEntry()->lundId();

      if(lund==PdtLund::pi_plus) _fille1Lab = BtaCandidate(*Cand); 
      else if(lund==PdtLund::pi_minus) _fille2Lab = BtaCandidate(*Cand); 
      else if(lund==PdtLund::pi0) _filleX0Lab = BtaCandidate(*Cand); 
      else if(lund==PdtLund::gamma && _fille1Lab.p()==0 ) _fille1Lab = BtaCandidate(*Cand); 
      else if(lund==PdtLund::gamma ) _fille2Lab = BtaCandidate(*Cand); 
      else if(lund==PdtLund::eta ) _filleX0Lab = BtaCandidate(*Cand); 
      else if(lund==PdtLund::rho0 ) _filleX0Lab = BtaCandidate(*Cand); 
    }

  int lundXu = _XuLab.pdtEntry()->lundId();
  if(lundXu==PdtLund::pi0) return; //No more family than the (one) two gamma(s) in this case


  //Petites filles
  BtaBooster* CMS = new BtaBooster(_Ups);
  if(_filleX0Lab.p()!=0)  _filleX0Ups = CMS->boostTo(_filleX0Lab);

  HepAListIterator<BtaCandidate> petiteFille = _filleX0Lab.daughterIterator();
  while ( 0 != (Cand = petiteFille()) )   
    {
      int lund = Cand->pdtEntry()->lundId();

      if(lund==PdtLund::pi_plus) _pfille1Lab = BtaCandidate(*Cand); 
      else if(lund==PdtLund::pi_minus) _pfille2Lab = BtaCandidate(*Cand); 
      else if(lund==PdtLund::pi0) _pfillePi0Lab = BtaCandidate(*Cand); 
      else if(lund==PdtLund::gamma && _pfille1Lab.p()==0 ) _pfille1Lab = BtaCandidate(*Cand); 
      else if(lund==PdtLund::gamma ) _pfille2Lab = BtaCandidate(*Cand); 
    }

  
  //Toutes petites filles 
  if(_pfillePi0Lab.p()==0) { delete CMS; return; }//Le seul cas ou cela est possible est: eta'->eta pi+ pi-, eta->pi+ pi- pi0, pi0->gamma gamma
  else _pfillePi0Ups = CMS->boostTo(_pfillePi0Lab);

  HepAListIterator<BtaCandidate> toutePetiteFille = _pfillePi0Lab.daughterIterator();
  while ( 0 != (Cand = toutePetiteFille()) )   
    {
      int lund = Cand->pdtEntry()->lundId();

      if(lund==PdtLund::gamma && _ppfilleGam1Lab.p()==0 ) _ppfilleGam1Lab = BtaCandidate(*Cand); 
      else if(lund==PdtLund::gamma ) _ppfilleGam2Lab = BtaCandidate(*Cand); 
    }

  delete CMS;
  return;
}


BtaCandidate
XSLRecoYAnalyzer::fittedB_Tree(EventInfo* evtInfo, Hep3Vector p3NuLab, BtaCandidate &fittedNu, BtaCandidate &fittedXu,  BtaCandidate &fittedLep,
			       BtaCandidate &f1, BtaCandidate &f2,  BtaCandidate &fX0)
{
  //The beam-energy constraint should only be used for decays of X from Ups(4S)->X Xbar, for example B0

  HepLorentzVector in = HepLorentzVector(0,0,0,0);
  BtaCandidate fittedB = BtaCandidate(in);

  BtaCandidate* newXu = XuToFit(evtInfo);
  if(newXu!=0)
    {
      //Lepton
      BtaCandidate newLep = BtaCandidate(_LepLab);
      newLep.invalidateFit();

      //Neutrino
      int lund = _LepLab.pdtEntry()->lundId();
      static const PdtEntry* Pdtnue = Pdt::lookup(PdtLund::nu_e);
      static const PdtEntry* PdtAntinue = Pdt::lookup(PdtLund::anti_nu_e);
      static const PdtEntry* Pdtnumu = Pdt::lookup(PdtLund::nu_mu);
      static const PdtEntry* PdtAntinumu = Pdt::lookup(PdtLund::anti_nu_mu);

      BtaCandidate newNu;
      if(lund==PdtLund::e_plus) newNu= BtaCandidate(p3NuLab,Pdtnue);
      else if(lund==PdtLund::e_minus) newNu= BtaCandidate(p3NuLab,PdtAntinue);
      else if(lund==PdtLund::mu_plus) newNu= BtaCandidate(p3NuLab,Pdtnumu);
      else if(lund==PdtLund::mu_minus) newNu= BtaCandidate(p3NuLab,PdtAntinumu);
      setMassConstraint(newNu);

      BtaOpMakeTree comb;
      BtaCandidate* toFit = comb.create(*newXu,newLep,newNu);      
      static const PdtEntry* PdtB0 = Pdt::lookup(PdtLund::B0);
      static const PdtEntry* PdtAntiB0 = Pdt::lookup(PdtLund::anti_B0);
      static const PdtEntry* PdtBplus = Pdt::lookup(PdtLund::B_plus);
      static const PdtEntry* PdtBminus = Pdt::lookup(PdtLund::B_minus);
      if(newLep.charge()>0 && newXu->charge()!=0 ) toFit->setType(PdtB0);
      else if(newLep.charge()<0 && newXu->charge()!=0 ) toFit->setType(PdtAntiB0);
      else if(newLep.charge()>0 && newXu->charge()==0 ) toFit->setType(PdtBplus);
      else if(newLep.charge()<0 && newXu->charge()==0 ) toFit->setType(PdtBminus);
      else cout<<"Problem in XSLRecoYAnalyzer::fittedB !!"<<endl;
      setMassConstraint(*toFit);
      setEnergyConstraint(*toFit,evtInfo) ;
      
      //if(evtInfo->primaryVertexCandidate()!=0) setBeamConstraint(*toFit,evtInfo->primaryVertexCandidate());
      setBeamConstraint(*toFit,evtInfo,true);  // <-- should try not giving it for >=2 trk modes
      //toFit->removeConstraint(BtaConstraint::Mass) ;
      //toFit->removeConstraint(BtaConstraint::Beam) ;
      //toFit->removeConstraint(BtaConstraint::Energy) ;
      toFit->invalidateFit(); //not needed?

      VtxFitterOper fitter(*toFit,VtxFitterOper::TreeFitter); //FastVtx,TreeFitter,Cascade,etc.
      fitter.fitAll();
      fittedB = fitter.getFitted(*toFit);

      fittedNu = fitter.getFitted(newNu);
      fittedXu = fitter.getFitted(*newXu);
      fittedLep = fitter.getFitted(newLep);  

      if(_mode=="eta2lnu"||_mode=="pi0lnu"||_mode=="rho0lnu") //Xu->Gam1/pi+ Gam2/pi- 
	{
	  f1 = fitter.getFitted(_fille1Lab);
	  f2 = fitter.getFitted(_fille2Lab);
	}
      else if(_mode=="eta3lnu"||_mode=="omegalnu"||_mode=="etaplnuE2PP"||_mode=="etaplnuE3PP") //Xu->pi+ pi- pi0/eta
	{
	  f1 = fitter.getFitted(_fille1Lab);
	  f2 = fitter.getFitted( _fille2Lab);
	  fX0 =fitter.getFitted(_filleX0LabForFit);    	  
	}
      else if(_mode=="etaplnuRG")
        {
	  f1 =  fitter.getFitted(_fille1Lab);
	  fX0 =  fitter.getFitted(_filleX0LabForFit);
	}
      else if(_mode=="rhoClnu")
	{
	  double ch = newXu->charge();
	  if(ch>0) f1 = fitter.getFitted(_fille1Lab);
	  else f2 =  fitter.getFitted(_fille2Lab);
	  fX0 =  fitter.getFitted(_filleX0LabForFit);   
	}

      delete newXu;
      delete toFit;
    }


  return fittedB; 
}

BtaCandidate
XSLRecoYAnalyzer::fittedY_Tree(EventInfo* evtInfo, BtaCandidate &fittedXu,  BtaCandidate &fittedLep,
			       BtaCandidate &f1, BtaCandidate &f2,  BtaCandidate &fX0)
{
  //The beam-energy constraint should only be used for decays of X from Ups(4S)->X Xbar, for example B0


  HepLorentzVector in = HepLorentzVector(0,0,0,0);
  BtaCandidate fittedY = BtaCandidate(in);

  BtaCandidate* newXu = XuToFit(evtInfo);
  if(newXu!=0)
    {
      //Lepton
      BtaCandidate newLep = BtaCandidate(_LepLab);
      newLep.invalidateFit();

      BtaOpMakeTree comb;
      BtaCandidate* toFit = comb.create(*newXu,newLep);      

      static const PdtEntry* PdtDstar0 = Pdt::lookup(PdtLund::D_star0);
      static const PdtEntry* PdtAntiDstar0 = Pdt::lookup(PdtLund::anti_D_star0);
      static const PdtEntry* PdtDstarplus = Pdt::lookup(PdtLund::D_star_plus);
      static const PdtEntry* PdtDstarminus = Pdt::lookup(PdtLund::D_star_minus);

      //Trick to account the lack of a neutrino (we pick D* only because it's a resonnance)
      if(newLep.charge()>0 && newXu->charge()!=0 ) toFit->setType(PdtDstar0);
      else if(newLep.charge()<0 && newXu->charge()!=0 ) toFit->setType(PdtAntiDstar0);
      else if(newLep.charge()>0 && newXu->charge()==0 ) toFit->setType(PdtDstarplus);
      else if(newLep.charge()<0 && newXu->charge()==0 ) toFit->setType(PdtDstarminus);
      else cout<<"Problem in XSLRecoYAnalyzer::fittedY_Tree !!"<<endl;

      if(_mode=="eta2lnu"||_mode=="pi0lnu") setBeamConstraint(*toFit,evtInfo,true);  // only for < 2 trk modes
      //toFit->removeConstraint(BtaConstraint::Beam);
      toFit->invalidateFit(); //not needed?

      VtxFitterOper fitter(*toFit,VtxFitterOper::TreeFitter); //FastVtx,TreeFitter,Cascade,etc.
      fitter.fitAll();
      fittedY = fitter.getFitted(*toFit);

      fittedXu = fitter.getFitted(*newXu);
      fittedLep = fitter.getFitted(newLep);  
     
      /*
      cout<<"originalXu.p4() "<<newXu->p3().x()<<" "<<newXu->p3().y()<<" "<<newXu->p3().z()<<" "<<newXu->energy()<<endl;
      cout<<"fittedXu.p4() "<<fittedXu.p3().x()<<" "<<fittedXu.p3().y()<<fittedXu.p3().z()<<fittedXu.energy()<<endl;
      cout<<"originalLep.p4() "<<newLep.p3().x()<<" "<<newLep.p3().y()<<" "<<newLep.p3().z()<<" "<<newLep.energy()<<endl;
      cout<<"fittedLep.p4() "<<fittedLep.p3().x()<<" "<<fittedLep.p3().y()<<" "<<fittedLep.p3().z()<<" "<<fittedLep.energy()<<endl;
      */

      if(_mode=="eta2lnu"||_mode=="pi0lnu"||_mode=="rho0lnu") //Xu->Gam1/pi+ Gam2/pi- 
	{
	  f1 = fitter.getFitted(_fille1Lab);
	  f2 = fitter.getFitted(_fille2Lab);
	}
      else if(_mode=="eta3lnu"||_mode=="omegalnu"||_mode=="etaplnuE2PP"||_mode=="etaplnuE3PP") //Xu->pi+ pi- pi0/eta
	{
	  f1 = fitter.getFitted(_fille1Lab);
	  f2 = fitter.getFitted( _fille2Lab);
	  fX0 =fitter.getFitted(_filleX0LabForFit);    	  
	}
      else if(_mode=="etaplnuRG")
        {
	  f1 =  fitter.getFitted(_fille1Lab);
	  fX0 =  fitter.getFitted(_filleX0LabForFit);
	}
      else if(_mode=="rhoClnu")
	{
	  double ch = newXu->charge();
	  if(ch>0) f1 = fitter.getFitted(_fille1Lab);
	  else f2 =  fitter.getFitted(_fille2Lab);
	  fX0 =  fitter.getFitted(_filleX0LabForFit);   
	}

      delete newXu;
      delete toFit;
    }


  return fittedY; 
}

double
XSLRecoYAnalyzer::ProbChi2Y_Fast()
{
  //The beam-energy constraint should only be used for decays of X from Ups(4S)->X Xbar, for example B0
  double probChi2=-666;
 
  BtaOpMakeTree comb;
  BtaCandidate* toFit=0;
  BtaCandidate* eta=0;

  int nD = _LepLab.nDaughters();
  if(nD!=0) return probChi2;

  if(_mode=="pilnu") toFit = comb.create(_XuLab,_LepLab);
  else if(_mode=="eta3lnu") toFit = comb.create(_fille1Lab,_fille2Lab,_LepLab);
  else if(_mode=="etaplnuRG") toFit = comb.create(_pfille1Lab,_pfille2Lab,_LepLab);
  else if(_mode=="etaplnuE2PP") toFit = comb.create(_fille1Lab,_fille2Lab,_LepLab);
  else if(_mode=="etaplnuE3PP") 
    {
      static const PdtEntry* PdtEta = Pdt::lookup(PdtLund::eta);
      eta = comb.create(_pfille1Lab,_pfille2Lab,_pfillePi0Lab); //I don't think the pi0 is used by the fitter...
      eta->setType(PdtEta);
      toFit = comb.create(*eta,_fille1Lab,_fille2Lab,_LepLab);
    }
  else if(_mode=="rhoClnu"&&_XuLab.charge()>0) toFit = comb.create(_fille1Lab,_LepLab); 
  else if(_mode=="rhoClnu"&&_XuLab.charge()<0) toFit = comb.create(_fille2Lab,_LepLab);
  else if(_mode=="rho0lnu") toFit = comb.create(_fille1Lab,_fille2Lab,_LepLab);
  else if(_mode=="omegalnu") toFit = comb.create(_fille1Lab,_fille2Lab,_LepLab);
  
  else return probChi2;
  
  static const PdtEntry* PdtB0 = Pdt::lookup(PdtLund::B0);
  static const PdtEntry* PdtAntiB0 = Pdt::lookup(PdtLund::anti_B0);
  static const PdtEntry* PdtBplus = Pdt::lookup(PdtLund::B_plus);
  static const PdtEntry* PdtBminus = Pdt::lookup(PdtLund::B_minus);
  if(_LepLab.charge()>0 && _XuLab.charge()!=0 ) toFit->setType(PdtB0);
  else if(_LepLab.charge()<0 && _XuLab.charge()!=0 ) toFit->setType(PdtAntiB0);
  else if(_LepLab.charge()>0 && _XuLab.charge()==0 ) toFit->setType(PdtBplus);
  else if(_LepLab.charge()<0 && _XuLab.charge()==0 ) toFit->setType(PdtBminus);
  else{
    cout<<"****** WHAT THE FUCK!!!???  ******"<<endl;
    cout<<"Lep Charge: "<<_LepLab.charge()<<endl;
    cout<<"Xu Charge: "<<_XuLab.charge()<<endl;
  }
  
  toFit->invalidateFit(); //invalidate all fits through the whole decay tree. This is MANDATORY.

  VtxFitterOper fitter(*toFit,VtxFitterOper::FastVtx); //FastVtx,TreeFitter,etc.
  fitter.fitAll();
  BtaCandidate fittedYFast = fitter.getFitted(*toFit);
  
  if(fittedYFast.decayVtx()!=0)
    {
      int nDofYFast=fittedYFast.decayVtx()->nDof();
      double chi2YFast=fittedYFast.decayVtx()->chiSquared();
      probChi2=probab(nDofYFast,chi2YFast);	 
    }
  
  if(eta!=0) delete eta;  
  if(toFit!=0) delete toFit;

  return probChi2; 
}

BtaCandidate
XSLRecoYAnalyzer::fittedXu(EventInfo* evtInfo,  BtaCandidate &f1, BtaCandidate &f2,  BtaCandidate &fX0)
{ 


  HepLorentzVector in = HepLorentzVector(0,0,0,0);
  BtaCandidate fittedXu = BtaCandidate(in);

  BtaCandidate* toFit = XuToFit(evtInfo);
  if(toFit!=0)
    {
      VtxFitterOper::algType algType  = VtxFitterOper::TreeFitter;
      //if(_mode=="pi0lnu"||_mode=="eta2lnu"||_mode=="rhoClnu") algType = VtxFitterOper::Cascade;
      //else algType = VtxFitterOper::TreeFitter;

      //if(evtInfo->primaryVertexCandidate()!=0) setBeamConstraint(*toFit,evtInfo->primaryVertexCandidate());
      if(_mode=="eta2lnu" || _mode =="pi0lnu" || _mode =="rhoClnu") setBeamConstraint(*toFit,evtInfo,true);  // only for < 2 trk modes

      VtxFitterOper fitter(*toFit,algType); //FastVtx,TreeFitter,Cascade,Add4,etc.
      fitter.fitAll();
      fittedXu = fitter.getFitted(*toFit);
      
      if(_mode=="eta2lnu"||_mode=="pi0lnu"||_mode=="rho0lnu") //Xu->Gam1/pi+ Gam2/pi- 
	{
	  f1 = fitter.getFitted(_fille1Lab);
	  f2 = fitter.getFitted(_fille2Lab);
	}
      else if(_mode=="eta3lnu"||_mode=="omegalnu"||_mode=="etaplnuE2PP"||_mode=="etaplnuE3PP") //Xu->pi+ pi- pi0/eta
	{
	  f1 = fitter.getFitted(_fille1Lab);
	  f2 = fitter.getFitted( _fille2Lab);
	  fX0 =fitter.getFitted(_filleX0LabForFit);    	  
	}
      else if(_mode=="etaplnuRG")
        {
	  f1 =  fitter.getFitted(_fille1Lab);
	  fX0 =  fitter.getFitted(_filleX0LabForFit);
	}
      else if(_mode=="rhoClnu")
	{
	  double ch = toFit->charge();
	  if(ch>0) f1 = fitter.getFitted(_fille1Lab);
	  else f2 =  fitter.getFitted(_fille2Lab);
	  fX0 =  fitter.getFitted(_filleX0LabForFit);   
	}

      delete toFit;
    }


  return fittedXu; 
}

BtaCandidate*
XSLRecoYAnalyzer::PiToFit(EventInfo* evtInfo)
{


  //The beam-energy constraint should only be used for decays of X from Ups(4S)->X Xbar, for example B0
  BtaCandidate* toFit = new BtaCandidate(_XuLab);
  toFit->invalidateFit();


  return toFit; //to be deleted by the caller (if not 0)
}

BtaCandidate*
XSLRecoYAnalyzer::Pi0ToFit(EventInfo* evtInfo)
{


  BtaCandidate* toFit = new BtaCandidate(_XuLab);
  //setMomentumConstraint(*toFit);
  setMassConstraint(*toFit);
  toFit->invalidateFit();


  return toFit; //to be deleted by the caller (if not 0)
}

BtaCandidate*
XSLRecoYAnalyzer::Eta2ToFit(EventInfo* evtInfo)
{


  BtaCandidate* toFit = new BtaCandidate(_XuLab);
  //setMomentumConstraint(*toFit);
  setMassConstraint(*toFit);
  toFit->invalidateFit();


  return toFit; //to be deleted by the caller (if not 0)
}

BtaCandidate*
XSLRecoYAnalyzer::Eta3ToFit(EventInfo* evtInfo)
{


  _filleX0LabForFit = BtaCandidate(_filleX0Lab);
  setMassConstraint(_filleX0LabForFit);
  _filleX0LabForFit.invalidateFit();

  BtaOpMakeTree comb;
  BtaCandidate* toFit = comb.create( _filleX0LabForFit,_fille1Lab,_fille2Lab);
  static const PdtEntry* PdtEta = Pdt::lookup(PdtLund::eta);
  toFit->setType(PdtEta);
  setMassConstraint(*toFit);
  toFit->invalidateFit(); //not needed?


  return toFit; //to be deleted by the caller (if not 0)
}

BtaCandidate*
XSLRecoYAnalyzer::EtapRGToFit(EventInfo* evtInfo)
{


  _filleX0LabForFit = BtaCandidate(_filleX0Lab);
  _filleX0LabForFit.invalidateFit();

  BtaOpMakeTree comb;
  BtaCandidate* toFit = comb.create(_filleX0LabForFit,_fille1Lab);
  static const PdtEntry* PdtEtap = Pdt::lookup(PdtLund::eta_prime);
  toFit->setType(PdtEtap);
  setMassConstraint(*toFit);
  toFit->invalidateFit(); //not needed?


  return toFit; //to be deleted by the caller (if not 0)
}

BtaCandidate*
XSLRecoYAnalyzer::EtapE3PPToFit(EventInfo* evtInfo)
{


  BtaCandidate pi0 = BtaCandidate(_pfillePi0Lab);
  setMassConstraint(pi0);
  pi0.invalidateFit();

  BtaOpMakeTree comb;
  BtaCandidate* eta = comb.create(pi0,_pfille1Lab,_pfille2Lab);
  static const PdtEntry* PdtEta = Pdt::lookup(PdtLund::eta);
  eta->setType(PdtEta);
  setMassConstraint(*eta);
  eta->invalidateFit(); //not needed?
  _filleX0LabForFit = BtaCandidate(*eta);
  delete eta;

  BtaCandidate* toFit = comb.create(_filleX0LabForFit,_fille1Lab,_fille2Lab);
  static const PdtEntry* PdtEtap = Pdt::lookup(PdtLund::eta_prime);
  toFit->setType(PdtEtap);
  setMassConstraint(*toFit);
  toFit->invalidateFit(); //not needed?


  return toFit; //to be deleted by the caller (if not 0)
}

BtaCandidate*
XSLRecoYAnalyzer::EtapE2PPToFit(EventInfo* evtInfo)
{


  _filleX0LabForFit = BtaCandidate(_filleX0Lab);
  setMassConstraint(_filleX0LabForFit);
  _filleX0LabForFit.invalidateFit();

  BtaOpMakeTree comb;
  BtaCandidate* toFit = comb.create(_filleX0LabForFit,_fille1Lab,_fille2Lab);
  static const PdtEntry* PdtEtap = Pdt::lookup(PdtLund::eta_prime);
  toFit->setType(PdtEtap);

  setMassConstraint(*toFit);
  toFit->invalidateFit();  //not needed?


  return toFit; //to be deleted by the caller (if not 0)
}

BtaCandidate*
XSLRecoYAnalyzer::RhoCToFit(EventInfo* evtInfo)
{


  _filleX0LabForFit = BtaCandidate(_filleX0Lab);
  //setMomentumConstraint(_filleX0LabForFit);
  setMassConstraint(_filleX0LabForFit);
  _filleX0LabForFit.invalidateFit();
  
  BtaOpMakeTree comb;
  BtaCandidate* toFit;
  if(_XuLab.charge()>0)
    { 
      toFit = comb.create(_filleX0LabForFit,_fille1Lab);
      static const PdtEntry* PdtRhop = Pdt::lookup(PdtLund::rho_plus);
      toFit->setType(PdtRhop);
    }
  else
    {
      toFit = comb.create(_filleX0LabForFit,_fille2Lab);
      static const PdtEntry* PdtRhom = Pdt::lookup(PdtLund::rho_minus);
      toFit->setType(PdtRhom);      
    }
  toFit->invalidateFit();


  return toFit; //to be deleted by the caller (if not 0)
}

BtaCandidate*
XSLRecoYAnalyzer::Rho0ToFit(EventInfo* evtInfo)
{


  BtaCandidate* toFit = new BtaCandidate(_XuLab);
  toFit->invalidateFit();


  return toFit; //to be deleted by the caller (if not 0)
}

BtaCandidate*
XSLRecoYAnalyzer::OmegaToFit(EventInfo* evtInfo)
{


  _filleX0LabForFit = BtaCandidate(_filleX0Lab);
  setMassConstraint(_filleX0LabForFit);
  _filleX0LabForFit.invalidateFit();

  BtaOpMakeTree comb;
  BtaCandidate* toFit = comb.create(_filleX0LabForFit,_fille1Lab,_fille2Lab);
  static const PdtEntry* PdtOme = Pdt::lookup(PdtLund::omega);
  toFit->setType(PdtOme);
  setMassConstraint(*toFit);
  toFit->invalidateFit();  


  return toFit; //to be deleted by the caller (if not 0)
}


BtaCandidate*
XSLRecoYAnalyzer::XuToFit(EventInfo* evtInfo)
{


  BtaCandidate* toFit=0;

  int nD = _XuLab.nDaughters();

  if(_mode=="pilnu") toFit = PiToFit(evtInfo);
  else if(_mode=="pi0lnu" && nD==2) toFit = Pi0ToFit(evtInfo);
  else if(_mode=="eta2lnu"&& nD==2) toFit = Eta2ToFit(evtInfo);
  else if(_mode=="eta3lnu") toFit = Eta3ToFit(evtInfo);
  else if(_mode=="etaplnuRG") toFit = EtapRGToFit(evtInfo); 
  else if(_mode=="etaplnuE2PP") toFit = EtapE2PPToFit(evtInfo); 
  else if(_mode=="etaplnuE3PP") toFit = EtapE3PPToFit(evtInfo); 
  else if(_mode=="rhoClnu") toFit = RhoCToFit(evtInfo); 
  else if(_mode=="rho0lnu") toFit = Rho0ToFit(evtInfo);
  else if(_mode=="omegalnu") toFit = OmegaToFit(evtInfo);

  //note: toFit==0 for merged pi0's (otherwise  fittedXu = fitter.getFitted(*toFit); apparently core dump...)
  return toFit; //to be deleted by the caller (if not 0)
}



BtaCandidate
XSLRecoYAnalyzer::fittedEtaL1(EventInfo* evtInfo )
{  


  int lund = _XuLab.pdtEntry()->lundId();
  int nD = _XuLab.nDaughters();
  if(!(lund==PdtLund::eta_prime&&nD==3)) return 0;


  BtaCandidate pi0 = BtaCandidate(_pfillePi0Lab);
  setOriginConstraint(pi0,evtInfo); //This is for the kinematic fit Add4 or Cascade
  setMassConstraint(pi0);
  pi0.invalidateFit();

  BtaOpMakeTree comb;
  BtaCandidate* eta = comb.create(pi0,_pfille1Lab,_pfille2Lab);
  eta->setType(Pdt::lookup(PdtLund::eta));
  setMassConstraint(*eta);
  eta->invalidateFit(); //not needed?
  
  VtxFitterOper fitter(*eta,VtxFitterOper::TreeFitter); //FastVtx,TreeFitter,etc.
  fitter.fitAll();
  BtaCandidate fittedEta = fitter.getFitted(*eta);

    
  delete eta;
  return fittedEta;
}


double
XSLRecoYAnalyzer::ProbChi2TreeToCutOn(EventInfo* evtInfo)
{

  //This fit is like fitted YTree, but without the Xu mass constraint that we're interested to use for our fit...
  //For pi0/eta2, it's probably not a good idea to use this...
  //For pi/rhoC/rho0, we reproduce the results of fittedYTree
  //For eta3/etap/omega, we define a variation of fittedYTree to cut on.

  double prob=-666;

  BtaCandidate* newXu = XuToFit(evtInfo);
  if(newXu!=0)
    {
      //the Key: removing Xu's mass contraint! :-)
      newXu->removeConstraint(BtaConstraint::Mass);

      //Lepton
      BtaCandidate newLep = BtaCandidate(_LepLab);
      newLep.invalidateFit();

      BtaOpMakeTree comb;
      BtaCandidate* toFit = comb.create(*newXu,newLep);      

      static const PdtEntry* PdtDstar0 = Pdt::lookup(PdtLund::D_star0);
      static const PdtEntry* PdtAntiDstar0 = Pdt::lookup(PdtLund::anti_D_star0);
      static const PdtEntry* PdtDstarplus = Pdt::lookup(PdtLund::D_star_plus);
      static const PdtEntry* PdtDstarminus = Pdt::lookup(PdtLund::D_star_minus);

      //Trick to account the lack of a neutrino (we pick D* only because it's a resonnance)
      if(newLep.charge()>0 && newXu->charge()!=0 ) toFit->setType(PdtDstar0);
      else if(newLep.charge()<0 && newXu->charge()!=0 ) toFit->setType(PdtAntiDstar0);
      else if(newLep.charge()>0 && newXu->charge()==0 ) toFit->setType(PdtDstarplus);
      else if(newLep.charge()<0 && newXu->charge()==0 ) toFit->setType(PdtDstarminus);
      else cout<<"Problem in XSLRecoYAnalyzer::fittedY_Tree !!"<<endl;

      if(_mode=="eta2lnu"||_mode=="pi0lnu") setBeamConstraint(*toFit,evtInfo,true);  // only for < 2 trk modes
      //toFit->removeConstraint(BtaConstraint::Beam);
      toFit->invalidateFit(); //not needed?

      VtxFitterOper fitter(*toFit,VtxFitterOper::TreeFitter); //FastVtx,TreeFitter,Cascade,etc.
      fitter.fitAll();
      BtaCandidate fittedY = fitter.getFitted(*toFit);
      
      int nDof=fittedY.decayVtx()->nDof();
      double chi2=fittedY.decayVtx()->chiSquared();
      prob=probab(nDof,chi2);	    

      delete newXu;
      delete toFit;
    }

  
  return prob;

}


BtaCandidate 
XSLRecoYAnalyzer::VDLabFromFittedDaughters(BtaCandidate f1, BtaCandidate f2,  BtaCandidate fX0)
{
  //This function is to be used to get VDLab of re-fitted Xu and daughters

  HepLorentzVector in = HepLorentzVector(0,0,0,0);
  BtaCandidate  VDLabBis = BtaCandidate(in) ;   
  //pilnu or merged pi0/eta...
  if(_XuLab.nDaughters()<1) {  VDLabBis = BtaCandidate(_XuLab); return VDLabBis ; }

  //From Leif Wilden (/u/ec/wilden/vub/RELEASE/XslUser/Vub.cc):
  //rho0 -> take pi+, rho+/- -> take pi0, 
  //omega -> 3pi: boost the pi+ and the pi- in the omega rest frame, then take the Xproduct
  //omega -> anything else: take the first pion, whatever it is 
  

  if(_mode=="pi0lnu"||_mode=="eta2lnu") 
    {
      bool f1Higher=(f1.p3().phi()>f2.p3().phi()) ? true:false;
      if(f1Higher) VDLabBis = BtaCandidate(f1);
      else VDLabBis = BtaCandidate(f2);
    }
  else if(_mode=="etaplnuRG") VDLabBis = BtaCandidate(f1);
  else if(_mode=="rhoClnu") VDLabBis = BtaCandidate(fX0);
  else if(_mode=="rho0lnu")  VDLabBis = BtaCandidate(f1);
  else if(_mode=="eta3lnu"||_mode=="etaplnuE2PP"||_mode=="etaplnuE3PP"||_mode=="omegalnu")
    {
      HepLorentzVector piPlusHad(f1.p4());
      HepLorentzVector piMinusHad(f2.p4());
      piPlusHad.boost(-_XuLab.p4().boostVector());
      piMinusHad.boost(-_XuLab.p4().boostVector());
      HepLorentzVector p4VDaugXuFrame = HepLorentzVector(piPlusHad.vect().cross(piMinusHad.vect()),0.);
      
      //boost back in the lab frame for uniformity with other decay modes...
      HepLorentzVector p4VDaugLabFrame(p4VDaugXuFrame);
      p4VDaugLabFrame.boost(_XuLab.p4().boostVector());
      
      //Create the BtaCandidate (with no type, but we don't care)
      VDLabBis = BtaCandidate(p4VDaugLabFrame);
    }
  else 
    {       
      HepAListIterator<BtaCandidate> dauIter(_XuLab.daughterIterator()); 
      VDLabBis = BtaCandidate(*dauIter()); 
      cout<<"Unexpected case in XSLRecoYAnalyzer::GetVDLabBis()"<<endl; 
    }    

  return VDLabBis;
}

double 
XSLRecoYAnalyzer::deltaE( HepLorentzVector PmissUps )
{
  double deltaE = _XuUps.energy()+_LepUps.energy()+PmissUps.rho()-_eBeamCM;
  return deltaE;
}

double 
XSLRecoYAnalyzer::deltaE2( HepLorentzVector PmissUps )
{
  double deltaE2 = _XuUps.energy()+_LepUps.energy()+PmissUps.e()-_eBeamCM;
  return deltaE2;
}

double 
XSLRecoYAnalyzer::delThetaNu( HepLorentzVector PmissUps )
{
  double csBY = cosBY();
  if(fabs(csBY)>1) return -666;
  double sinBY = sqrt(1-csBY*csBY);
  
  double pNu = _eBeamCM - _YUps.energy();

  double thYMiss = _YUps.p3().angle(PmissUps.vect());
  double cosBMiss = cos(thYMiss)*csBY+sin(thYMiss)*sinBY;

  double num = PmissUps.rho()*pBCM()*cosBMiss - PmissUps.vect()*_YUps.p3();
  double denom = PmissUps.rho()*pNu;
  if(denom==0) return -665;

  double delThNu = num/denom;
  return delThNu;
}

double 
XSLRecoYAnalyzer::mES( HepLorentzVector PmissUps )
{
  double eBeam2 = _eBeamCM*_eBeamCM;
  
  Hep3Vector p = PmissUps.vect()+_XuUps.p3()+_LepUps.p3();

  double arg = eBeam2 - p.mag2();
  double mES = sqrt(fabs(arg));
  
  if(arg>=0) return mES;
  else return -mES;
}
 





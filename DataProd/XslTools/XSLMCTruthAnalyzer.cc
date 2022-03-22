//*****************************************************************
//   
//   Creation: David Cote, Universite de Montreal, 12/11/03
//
#include "BaBar/BaBar.hh"



#include "XslTools/XSLMCTruthAnalyzer.hh"

#include <math.h>
#include "ErrLogger/ErrLog.hh"
#include "AbsEnv/AbsEnv.hh"
#include "GenEnv/GenEnv.hh"
#include "PDT/Pdt.hh"

#include "PepEnv/PepEnv.hh"
#include "PepEnvData/PepBeams.hh"
#include "BetaCoreTools/BtaOpMakeTree.hh"
#include "BetaCoreTools/BtaPrintTree.hh"
using std::cout;
using std::endl;
using std::string;

XSLMCTruthAnalyzer::XSLMCTruthAnalyzer() 
{
  _eventType=-666;
  _typeB1=-666;
  _typeB2=-666;
  _modeB1="unknown";
  _modeB2="unknown";
}

XSLMCTruthAnalyzer::XSLMCTruthAnalyzer(HepAList<BtaCandidate> inputMcList)
{ 
  _mcList=inputMcList;
  AnalyzeEvent();
  AnalyzeMissingMomentum();
}


XSLMCTruthAnalyzer::~XSLMCTruthAnalyzer()
{}


void
XSLMCTruthAnalyzer::AnalyzeEvent()
{
  //1-599: ee -> Ups(4S)
  //600: ccbar
  //700: uds
  //800: tautau
  //900: mumu
  //1000: Bhabha
  //1100: gamma gamma
  //1200: Radiative Bhabhas
  //1300: other(??)
  _eventType=-666;

  _typeB1=-666;
  _typeB2=-666;
  _modeB1="unknown";
  _modeB2="unknown";

  //generic multi-purpose BtaCandidate*:
  BtaCandidate* TruthCand;

  HepAListIterator<BtaCandidate> iterMC(_mcList);

  //So, we first get the Ups(4S)...
  bool UpsFound=false;
  while ( 0 != (TruthCand = iterMC()) && UpsFound==false )         
    {
      if(TruthCand->pdtEntry()->lundId()==PdtLund::Upsilon_4S) 
	{ _Ups4S = BtaCandidate(*TruthCand ); UpsFound=true; }
    }

  if(UpsFound)
    {
      _eventType=1;

      //Then we get the two B's...
      int p=1;
      HepAListIterator<BtaCandidate> UpsDaughter = _Ups4S.daughterIterator();
      while ( 0 != (TruthCand = UpsDaughter()) )   
	{
	  bool isB=( TruthCand->pdtEntry()->lundId()==PdtLund::B0
		     || TruthCand->pdtEntry()->lundId()==PdtLund::anti_B0
		     || TruthCand->pdtEntry()->lundId()==PdtLund::B_plus 
		     || TruthCand->pdtEntry()->lundId()==PdtLund::B_minus );
	    
	  if(p==1 && isB) { _B1 = BtaCandidate(*TruthCand); }
	  if(p==2 && isB) { _B2 = BtaCandidate(*TruthCand); }
	  if(isB) p+=1;
	}

      //BType also sets _Xu,_lep,_nu if B->pi/pi0/eta/etap/rhoC/rho0/omega l nu
      //Note: the order B1 and then B2 is important for the XuFamily to be set for B2 by default in case of two sig B's
      _typeB1 = BType(&_B1,_Xu1_Lab,_Lep1_Lab,_Nu1_Lab);
      if(_typeB1<=77)
	{
	  _modeB1 = FigureOutXuFamily(_Xu1_Lab, _typeB1);  //This sets Xu1 family
	  _VDaug1_Lab = GetXuDaughter(_Xu1_Lab, _modeB1);  //GetXuDaughter needs FigureOutXuFamily first
	}

      _typeB2 = BType(&_B2,_Xu2_Lab,_Lep2_Lab,_Nu2_Lab);
      if(_typeB2<=77)
	{
	  _modeB2 = FigureOutXuFamily(_Xu2_Lab, _typeB2);   //This could reset the family from Xu1 to Xu2
	  _VDaug2_Lab = GetXuDaughter(_Xu2_Lab, _modeB2);  //GetXuDaughter needs FigureOutXuFamily first
	}

    }//if(UpsFound)
  else
    {      
      //1-599: ee -> Ups(4S)
      //600: ccbar
      //700: uds
      //800: tautau
      //900: mumu
      //1000: Bhabha
      //1100: gamma gamma
      //1200: Radiative Bhabhas
      //1300: other(??)
      
      bool cont=true;
      bool gammaFound=false;
      bool eFound=false;
      BtaCandidate* top=0;
      iterMC.rewind();
      while ( 0 != (TruthCand = iterMC()) && cont )         
	{
	  //Get a BtaCandidate without mother...
	  bool again=true;
	  top=TruthCand;
	  while(again) 
	    { 
	      if(top->theMother()!=0) 
		{
		  if(top->theMother()->pdtEntry()->lundId()==PdtLund::vpho) { again=false; }
		  else { top=top->theMother(); again=true; } 
		}
	      else{ again=false; }
	    } 
	  int topLund = abs( (int)top->pdtEntry()->lundId() );
	  
	  //BtaPrintTree printTree;
	  //cout << printTree.print(*top) << endl;		

	  if(topLund==PdtLund::c) { _eventType=600; cont=false; }
	  else if(topLund==PdtLund::d || topLund==PdtLund::u || topLund==PdtLund::s) 
	    { _eventType=700; cont=false; }
	  else if(topLund==PdtLund::tau_minus) { _eventType=800; cont=false; }
	  else if(topLund==PdtLund::mu_minus) { _eventType=900; cont=false; }
	  else if(topLund==PdtLund::e_minus) { eFound=true; cont=true; }
	  else if(topLund==PdtLund::gamma) { gammaFound=true; cont=true; }
	  else cont=true; 
	}

      if(cont==true)
	{
	  if(eFound&&!gammaFound) _eventType=1000;
	  else if(!eFound&&gammaFound) _eventType=1100;
	  else if(eFound&&gammaFound) _eventType=1200;
	  else  _eventType=1300; 
	}
    }//else if(!UpsFound)

  return;
}


int
XSLMCTruthAnalyzer::BType( BtaCandidate* B, BtaCandidate &XuLab, BtaCandidate &LepLab, BtaCandidate &NuLab)
{
  int type=-666;
  // -77 a 77: B-> signal mode  (negative means isBadXuDecayMode)  
  //88: B-> other Xu l nu
  //100-299: B-> Xc l nu
  //300: B-> Xu tau nu_tau
  //350: B-> Xc tau nu
  //400: B-> hadronic
  //500: B-> leptonic
  //550: B-> Other(baryons,??)

  bool LepFound=false,TauFound=false,XudFound=false,XscFound=false,nonResXu=false,JETSET=false,muFound=false;


  BtaCandidate *tau=0;
  BtaCandidate* TruthCand;
  HepAListIterator<BtaCandidate> BDaughter = B->daughterIterator();
  while ( 0 != (TruthCand = BDaughter()) )   
    {
      int lund = abs( (int)TruthCand->pdtEntry()->lundId() );
      int quark = lund%1000;

      if(lund<5) JETSET=true;
      else if(lund==11) LepFound=true; 
      else if(lund==13) { LepFound=true; muFound=true; }
      else if(lund==15) { TauFound=true; tau = new BtaCandidate(*TruthCand); }
      else if(lund==41 || lund==42){ nonResXu=true; }

      //exceptions to the "quark" last three digits system...
      else if(lund==130) XscFound=true;  //K0L
      else if(lund==331) XudFound=true;  //eta_prime
      else if(lund==4101||lund==4103||lund==2101||lund==2103||lund==4301||lund==4303||lund==2203) JETSET=true; 
      //JETSET "exotic" particles for B decays:
      // cd_0, cd_1, ud_0, ud_1, cs_0, cs_1, uu_1 
      // 4101, 4103, 2101, 2103, 4301, 4303, 2203

      //The last three digits has to be between 100-299
      else if(quark>=100&&quark<=299) XudFound=true; 
      //The last three digits has to be between 300-499
      else if(quark>=300&&quark<=499) XscFound=true; 
    }

  if(LepFound&&!TauFound&&XudFound&&!XscFound) type=88;  //B->Xulnu
  else if(LepFound&&!TauFound&&!XudFound&&XscFound) type=100; //B->Xc l nu
  else if(LepFound&&!TauFound&&XudFound&&XscFound){ //B->Xc Xu l nu
    if(muFound) type=201; 
    else type=101; 
  }
  else if(!LepFound&&TauFound&&XudFound&&!XscFound) type=300; //B->Xu tau nu
  else if(!LepFound&&TauFound&&!XudFound&&XscFound) type=350; //B->Xc tau nu
  else if(!LepFound&&!TauFound&&(XudFound||XscFound)) type=400; //B->hadrons
  else if(LepFound&&!TauFound&&!XudFound&&!XscFound){
    if(nonResXu){ type=500; } //B->nonResXu l nu
    else{ type=510; } //B->l nu (gamma) //not supposed to be generated in generic SP4-8
  }
  else if(JETSET){ type=520; }
  else type=550; //B->other

  //AnalyzeBtoXulnu returns 1-77 type. It can also return -type for BadXuDecayMode.
  //Otherwise, type stays =88
  if(type==88) type=AnalyzeBtoXulnu(B, XuLab, LepLab, NuLab, type); 


  // B->D(*) e nu_e = 1xx
  // B->D(*) mu nu_mu = 2xx
  //AnalyzeBtoXclnu returns (for the mu case): 
  //                    200  :  B->otherXc l nu, 
  //                    201  :  B->Xc Xu l nu, 
  //                    210  :  B->D0 l nu, D0->K pi
  //                    211  :  B->D0 l nu, D0->K pi pi0
  //                    212  :  B->D0 l nu, D0->K 3pi
  //                    213  :  B->D0 l nu, D0->e X
  //                    214  :  B->D0 l nu, D0->mu X
  //                    215  :  B->D0 l nu, D0->other
  //                    220  :  B->D0* l nu, D0*->D0 pi0, D0->K pi
  //                    221  :  B->D0* l nu, D0*->D0 pi0, D0->K pi pi0
  //                    222  :  B->D0* l nu, D0*->D0 pi0, D0->K 3pi
  //                    223  :  B->D0* l nu, D0*->D0 pi0, D0->e X
  //                    224  :  B->D0* l nu, D0*->D0 pi0, D0->mu X
  //                    225  :  B->D0* l nu, D0*->D0 pi0, D0->other
  //                    230  :  B->D0* l nu, D0*->D0 gamma, D0->K pi
  //                    231  :  B->D0* l nu, D0*->D0 gamma, D0->K pi pi0
  //                    232  :  B->D0* l nu, D0*->D0 gamma, D0->K 3pi
  //                    233  :  B->D0* l nu, D0*->D0 gamma, D0->e X
  //                    234  :  B->D0* l nu, D0*->D0 gamma, D0->mu X
  //                    235  :  B->D0* l nu, D0*->D0 gamma, D0->other
  //                    240  :  B->D*+ l nu, D*+->D0 pi+, D0->K pi
  //                    241  :  B->D*+ l nu, D*+->D0 pi+, D0->K pi pi0
  //                    242  :  B->D*+ l nu, D*+->D0 pi+, D0->K 3pi
  //                    243  :  B->D*+ l nu, D*+->D0 pi+, D0->e X
  //                    244  :  B->D*+ l nu, D*+->D0 pi+, D0->mu X
  //                    245  :  B->D*+ l nu, D*+->D0 pi+, D0->other
  //                    256  :  B->D*+ l nu, D*+->D+ pi0, D+->K- pi+ pi+
  //                    257  :  B->D*+ l nu, D*+->D+ pi0, D+->e X
  //                    258  :  B->D*+ l nu, D*+->D+ pi0, D+->mu X
  //                    259  :  B->D*+ l nu, D*+->D+ pi0, D+->other
  //                    266  :  B->D*+ l nu, D*+->D+ gamma, D+->K- pi+ pi+
  //                    267  :  B->D*+ l nu, D*+->D+ gamma, D+->e X
  //                    268  :  B->D*+ l nu, D*+->D+ gamma, D+->mu X
  //                    269  :  B->D*+ l nu, D*+->D+ gamma, D+->other
  //                    276  :  B->D+ l nu, D+->K- pi+ pi+
  //                    277  :  B->D+ l nu, D+->e X
  //                    278  :  B->D+ l nu, D+->mu X
  //                    279  :  B->D+ l nu, D+->other
  //                    280  :  B->D**+ l nu (other)
  //                    281  :  B->D1+ l nu
  //                    282  :  B->D2*+ l nu
  //                    283  :  B->D_0*+ l nu
  //                    284  :  B->D1'+ l nu
  //                    290  :  B->D**0 l nu (other)
  //                    291  :  B->D10 l nu
  //                    292  :  B->D2*0 l nu
  //                    293  :  B->D_0*0 l nu
  //                    294  :  B->D1'0 l nu
  else if(type==100) type=AnalyzeBtoXclnu(B, type);


  //AnalyzeTau Returns 0,1,2 
  //                    300  :  B->Xu Tau, Tau->other
  //                    301  :  B->Xu Tau, Tau->X e
  //                    302  :  B->Xu Tau, Tau->X mu
  //                    350  :  B->Xc Tau, Tau->other
  //                    351  :  B->Xc Tau, Tau->X e
  //                    352  :  B->Xc Tau, Tau->X mu
  else if(type==300||type==350) type+=AnalyzeTau(tau);


  //                    400   : B->Xu X 
  //                    410   : B->Xu X, X->Xu Y 
  //                    420   : B->Xu X, X->e Y 
  //                    425   : B->Xu X, X->Xu e Y 
  //                    430   : B->Xu X, X->mu Y 
  //                    435   : B->Xu X, X->Xu mu Y 
  //                    450   : B->X, X->Xu Y  
  //                    460   : B->X, X->e Y 
  //                    465   : B->X, X->Xu e Y 
  //                    470   : B->X, X->mu Y 
  //                    475   : B->X, X->Xu mu Y 
  //                    490   : B->other hadrons 
  else if(type==400) type+=AnalyzeBtoHadrons(B);

  //                    500   : B->Xu0/Xu_plus l nu
  //                    510   : B->l nu (gamma) 
  //                    520   : B->JETSET 

  if(tau!=0) delete tau;
  return type;
}

int
XSLMCTruthAnalyzer::AnalyzeBtoHadrons( BtaCandidate* B )
{
  bool BtoXu=false,XtoXu=false,XtoElectron=false,XtoMuon=false;
  BtaCandidate *TruthCand,*TruthCand2;

  HepAListIterator<BtaCandidate> BDaughter = B->daughterIterator();
  while ( 0 != (TruthCand = BDaughter()) )   
    {
      int lund = abs( (int)TruthCand->pdtEntry()->lundId() );
      if(lund==PdtLund::pi_plus || lund==PdtLund::pi0
	 || lund==PdtLund::eta || lund==PdtLund::eta_prime 
	 || lund==PdtLund::rho_plus || lund==PdtLund::rho0 
	 || lund==PdtLund::omega ) BtoXu=true;  //Impossible to decay to lepton in this case 
      else
	{
	  HepAListIterator<BtaCandidate> XDaughter = TruthCand->daughterIterator();
	  while ( 0 != (TruthCand2 = XDaughter()) )   
	    {
	      int lund2 = abs( (int)TruthCand2->pdtEntry()->lundId() );
	      if(lund2==PdtLund::pi_plus || lund2==PdtLund::pi0
		 || lund2==PdtLund::eta || lund2==PdtLund::eta_prime 
		 || lund2==PdtLund::rho_plus || lund2==PdtLund::rho0 
		 || lund2==PdtLund::omega ) XtoXu=true;  //Impossible to decay to lepton in this case 
	      else if(lund2==PdtLund::e_minus) XtoElectron=true;
	      else if(lund2==PdtLund::mu_minus) XtoMuon=true;
	    }// while ( 0 != (TruthCand2 = XDaughter()) )   
	}//else
    }// while ( 0 != (TruthCand = BDaughter()) )



  //         00   : B->Xu X 
  //         10   : B->Xu X, X->Xu Y 
  //         20   : B->Xu X, X->e Y 
  //         25   : B->Xu X, X->Xu e Y 
  //         30   : B->Xu X, X->mu Y 
  //         35   : B->Xu X, X->Xu mu Y 
  //         50   : B->X, X->Xu Y  
  //         60   : B->X, X->e Y 
  //         65   : B->X, X->Xu e Y 
  //         70   : B->X, X->mu Y 
  //         75   : B->X, X->Xu mu Y 
  //         90   : B->other hadrons 
  int t=90;

  if(BtoXu && !XtoXu && !XtoElectron && !XtoMuon) t=0;
  else if(BtoXu && XtoXu && !XtoElectron && !XtoMuon) t=10;
  else if(BtoXu && !XtoXu && XtoElectron && !XtoMuon) t=20;
  else if(BtoXu && XtoXu && XtoElectron && !XtoMuon) t=25;
  else if(BtoXu && !XtoXu && !XtoElectron && XtoMuon) t=30;
  else if(BtoXu && XtoXu && !XtoElectron && XtoMuon) t=35;

  else if(!BtoXu && XtoXu && !XtoElectron && !XtoMuon) t=50;
  else if(!BtoXu && !XtoXu && XtoElectron && !XtoMuon) t=60;
  else if(!BtoXu && XtoXu && XtoElectron && !XtoMuon) t=65;
  else if(!BtoXu && !XtoXu && !XtoElectron && XtoMuon) t=70;
  else if(!BtoXu && XtoXu && !XtoElectron && XtoMuon) t=75;
  else t=90;

  return t;
}


int
XSLMCTruthAnalyzer::AnalyzeTau( BtaCandidate* tau )
{
  int t=0;
  BtaCandidate* TruthCand;

  HepAListIterator<BtaCandidate> TauDaughter = tau->daughterIterator();
  while ( 0 != (TruthCand = TauDaughter()) && t==0 )   
    {
      int lund = abs( (int)TruthCand->pdtEntry()->lundId() );
      if(lund==11) t=1; 
      else if(lund==13) t=2; 
    }

  return t;
}


int
XSLMCTruthAnalyzer::AnalyzeBtoXclnu( BtaCandidate* B, int type )
{
  if(type!=100){ cout<<"Something's WRONG in XSLMCTruthAnalyzer::AnalyzeBtoXclnu !!!!!!!!"<<endl; exit(1); }

  BtaCandidate* TruthCand;
  BtaCandidate* D=0;
  bool DFound=false,LepFound=false,cont=true;
  int lundLep = 0;

  HepAListIterator<BtaCandidate> BDaughter = B->daughterIterator();
  while ( 0 != (TruthCand = BDaughter()) && cont )   
    {
      int lund = abs( (int)TruthCand->pdtEntry()->lundId() );
      int quark = lund%1000;
      
      if(lund==11||lund==13) { lundLep=lund; LepFound=true; }
      else if(quark>=400&&quark<=499) 
	{ D = new BtaCandidate(*TruthCand); DFound=true; }
      
      if(DFound&&LepFound) cont=false;
    }

  if(lundLep==13) type=200;
  else if(lundLep!=11){ cout<<"WRONG lepton in XSLMCTruthAnalyzer::AnalyzeBtoXclnu !!!!!!!!"<<endl; exit(1); }

  if(DFound) type+=AnalyzeD( D );

  if(D!=0) delete D;
  return type;
}

int
XSLMCTruthAnalyzer::AnalyzeD( BtaCandidate* D )
{
  //t goes from 0 to 99
  int t=0;

  BtaCandidate *D0=0,*Dplus=0;
  BtaCandidate* TruthCand;
  int lundD = abs( (int)D->pdtEntry()->lundId() );


  //Longue sequence of if{} else if{} pour determiner le chiffre des dizaines
  if(lundD==PdtLund::D0) { t=10; D0 = new BtaCandidate( *D ); }
  else if(lundD==PdtLund::D_plus) { t=70; Dplus = new BtaCandidate( *D ); }
  else if(lundD==PdtLund::D_star0)
    {
      HepAListIterator<BtaCandidate> D_Daughter = D->daughterIterator();
      while ( 0 != (TruthCand = D_Daughter()) )   
	{
	  int lund = abs( (int)TruthCand->pdtEntry()->lundId() );
	  if(lund==PdtLund::D0) D0 = new BtaCandidate( *TruthCand ); 
	  else if(lund==PdtLund::pi0) t=20;
	  else if(lund==PdtLund::gamma && D->nDaughters()==2) t=30;
	}
    }// else if(lundD==PdtLund::D_star0)
  else if(lundD==PdtLund::D_star_plus )
    {
      bool cont=true;
      HepAListIterator<BtaCandidate> D_Daughter = D->daughterIterator();
      while ( 0 != (TruthCand = D_Daughter()) && cont )   
	{
	  int lund = abs( (int)TruthCand->pdtEntry()->lundId() );
	  if(lund==PdtLund::D0) 
	    { 
	      D0 = new BtaCandidate( *TruthCand ); 
	      t=40; 
	      cont=false; 
	    }
	  else if(lund==PdtLund::D_plus) Dplus = new BtaCandidate( *TruthCand ); 
	  else if(lund==PdtLund::pi0) t=50;
	  else if(lund==PdtLund::gamma && D->nDaughters()==2) t=60;
	}
    }// else if(lundD==PdtLund::D_star_plus)
  else if(abs((int)D->charge())==1)
    {
      //D**+/-, we ignore the D decay in this case
      if(lundD==PdtLund::D_1_plus) t=81;
      else if(lundD==PdtLund::D_2_star_plus) t=82;
      else if(lundD==PdtLund::D_0_star_plus) t=83;
      else if(lundD==PdtLund::D_prime_1_plus) t=84;
      else t=80;
    }// else if(abs(D->charge())==1)
  else if(D->charge()==0)  
    {
      //D**0, we ignore the D decay in this case
      if(lundD==PdtLund::D_10) t=91;
      else if(lundD==PdtLund::D_2_star0) t=92;
      else if(lundD==PdtLund::D_0_star0) t=93;
      else if(lundD==PdtLund::D_prime_10) t=94;
      else t=90;
    }// else if(D->charge()==0)

  //Final touch (pour les unites): the D0/D+ decay modes...
  if(D0!=0) { t+= AnalyzeD0(D0); delete D0; }
  if(Dplus!=0) { t+= AnalyzeDplus(Dplus); delete Dplus; }
  return t;
}

int
XSLMCTruthAnalyzer::AnalyzeD0( BtaCandidate* D0 )
{
  //t goes from 0 to 5
  //                    210  :  B->D0 l nu, D0->K pi
  //                    211  :  B->D0 l nu, D0->K pi pi0
  //                    212  :  B->D0 l nu, D0->K 3pi
  //                    213  :  B->D0 l nu, D0->e X
  //                    214  :  B->D0 l nu, D0->mu X
  //                    215  :  B->D0 l nu, D0->other
  int t=0;
  int nKp=0,npip=0,npi0=0,ne=0,nmu=0,no=0;

  BtaCandidate* TruthCand;
  HepAListIterator<BtaCandidate> D0_Daughter = D0->daughterIterator();
  while ( 0 != (TruthCand = D0_Daughter()) )   
    {
      int lund = abs( (int)TruthCand->pdtEntry()->lundId() );
      if(lund==PdtLund::K_plus) nKp+=1; 
      else if(lund==PdtLund::pi_plus) npip+=1; 
      else if(lund==PdtLund::pi0) npi0+=1; 
      else if(lund==PdtLund::e_minus) ne+=1; 
      else if(lund==PdtLund::mu_minus) nmu+=1; 
      else if(lund!=PdtLund::gamma) no+=1;       
      //we don't end the while even if there's an unknown particle at this point in case there's a lepton later on
    }
  

  //We assume that charge conservation is OK!
  if(ne>0) t=3;
  else if(nmu>0) t=4;
  else if(nKp==1&&npip==1&&npi0==0&&ne==0&&nmu==0&&no==0) t=0;
  else if(nKp==1&&npip==1&&npi0==1&&ne==0&&nmu==0&&no==0) t=1;
  else if(nKp==1&&npip==3&&npi0==0&&ne==0&&nmu==0&&no==0) t=2;
  else t=5;
    
  return t;
}

int
XSLMCTruthAnalyzer::AnalyzeDplus( BtaCandidate* Dplus )
{
  //t goes from 6 to 9
  //                    276  :  B->D+ l nu, D+->K- pi+ pi+
  //                    277  :  B->D+ l nu, D+->e X
  //                    278  :  B->D+ l nu, D+->mu X
  //                    279  :  B->D+ l nu, D+->other
  int t=0;
  int nKp=0,npip=0,ne=0,nmu=0,no=0;
  double Kp=-1,pip=1;
  if(Dplus->pdtEntry()->lundId()==PdtLund::D_minus) { Kp=1; pip=-1; }

  BtaCandidate* TruthCand;
  HepAListIterator<BtaCandidate> Dp_Daughter = Dplus->daughterIterator();
  while ( 0 != (TruthCand = Dp_Daughter()) )   
    {
      int lund = abs( (int)TruthCand->pdtEntry()->lundId() );
      if(lund==Kp) nKp+=1; 
      else if(lund==pip) npip+=1; 
      else if(lund==PdtLund::e_minus) ne+=1; 
      else if(lund==PdtLund::mu_minus) nmu+=1; 
      else if(lund!=PdtLund::gamma) no+=1;       
      //we don't end the while even if there's an unknown particle at this point in case there's a lepton later on
    }
  

  //We assume that charge conservation is OK!
  if(ne>0) t=7;
  else if(nmu>0) t=8;
  else if(nKp==1&&npip==2&&ne==0&&nmu==0&&no==0) t=6;
  else t=9;
  
  return t;
}

int
XSLMCTruthAnalyzer::AnalyzeBtoXulnu( BtaCandidate* B, BtaCandidate &XuLab, BtaCandidate &LepLab, BtaCandidate &NuLab, int type )
{
  BtaCandidate* TruthCand;
  BtaCandidate* Xu=0;
  BtaCandidate* lep=0;
  BtaCandidate* nu=0;
  
  int nXud=0,nLep=0,nNu=0;
  bool maybeSignal=true;
  HepAListIterator<BtaCandidate> BDaughter = B->daughterIterator();
  while ( 0 != (TruthCand = BDaughter()) && maybeSignal )   
    {
      int lund = abs( (int)TruthCand->pdtEntry()->lundId() );

      //Warning! The Xu, lep, nu are set even though this B might not be signal
      //One has to ask typeB==1-77 before taking Xu, lep, nu!
      if(lund==11||lund==13) { lep = new BtaCandidate(*TruthCand); nLep+=1; }
      else if(lund==12||lund==14) { nu = new BtaCandidate(*TruthCand); nNu+=1; }
      else if(lund==111||lund==211||lund==221||lund==331||lund==113||lund==213||lund==223) 
	{ Xu = new BtaCandidate(*TruthCand); nXud+=1; }
      else if(lund!=22) maybeSignal=false;

      if(nLep>1||nNu>1||nXud>1) maybeSignal=false;
    }

  if( maybeSignal && abs( (int)lep->pdtEntry()->lundId()+ (int)nu->pdtEntry()->lundId() )==1 )
    {
      int XudLund= abs( (int)Xu->pdtEntry()->lundId() );
      int lepLund= abs( (int)lep->pdtEntry()->lundId() );
      double chrg = abs((int)Xu->charge()+(int)lep->charge()); //a little bit of redundancy doesn't hurt... ;-)

      if(XudLund==PdtLund::pi_plus && lepLund==PdtLund::e_minus && chrg==0) type=1;
      else if(XudLund==PdtLund::pi_plus && lepLund==PdtLund::mu_minus && chrg==0) type=11;
      else if(XudLund==PdtLund::pi0 && lepLund==PdtLund::e_minus && chrg==1) type=2;
      else if(XudLund==PdtLund::pi0 && lepLund==PdtLund::mu_minus && chrg==1) type=22;
      else if(XudLund==PdtLund::eta && lepLund==PdtLund::e_minus && chrg==1) type=3;
      else if(XudLund==PdtLund::eta && lepLund==PdtLund::mu_minus && chrg==1) type=33;
      else if(XudLund==PdtLund::eta_prime && lepLund==PdtLund::e_minus && chrg==1) type=4;
      else if(XudLund==PdtLund::eta_prime && lepLund==PdtLund::mu_minus && chrg==1) type=44;
      else if(XudLund==PdtLund::rho_plus && lepLund==PdtLund::e_minus && chrg==0) type=5;
      else if(XudLund==PdtLund::rho_plus && lepLund==PdtLund::mu_minus && chrg==0) type=55;
      else if(XudLund==PdtLund::rho0 && lepLund==PdtLund::e_minus && chrg==1) type=6;
      else if(XudLund==PdtLund::rho0 && lepLund==PdtLund::mu_minus && chrg==1) type=66;
      else if(XudLund==PdtLund::omega && lepLund==PdtLund::e_minus && chrg==1) type=7;
      else if(XudLund==PdtLund::omega && lepLund==PdtLund::mu_minus && chrg==1) type=77;
      else cout<<"PROBLEM IN MCTRUTH ANALYZER::AnalyzeSigB !!!!!!!!!!!!!!!!!"<<endl;
    }

  if(type==2||type==22||type==3||type==33||type==4||type==44||type==7||type==77) 
    {
      bool BadXuDecayMode=false;
      int XuLund = (int)Xu->pdtEntry()->lundId();
      BadXuDecayMode=IsBadXuDecayMode(Xu,XuLund);
      if(BadXuDecayMode) type= -type;
    }
  
  if(type<=77) 
    {
      XuLab = BtaCandidate(*Xu);
      LepLab = BtaCandidate(*lep);
      NuLab = BtaCandidate(*nu);
    }

  if(Xu!=0) delete Xu;
  if(lep!=0) delete lep;
  if(nu!=0)  delete nu;

  return type;
}


bool
XSLMCTruthAnalyzer::IsBadXuDecayMode( BtaCandidate* Xu, int XuLund )
{
  //A "Bad decay mode" is a decay mode that we don't try to reconstruct with real data
  bool IsBadDecayMode=false;

  //First taking care of the pi0 -> e+ e- Gamma case...
  int nD = Xu->nDaughters();
  if(XuLund==111 && nD!=2) IsBadDecayMode=true;

  //Now the other cases...
  BtaCandidate* Cand;
  BtaCandidate Eta;
  HepAListIterator<BtaCandidate> XuDaughter = Xu->daughterIterator();

  //Look for omega->pi+ pi- pi0 (gamma(s))
  //look for eta->gamma gamma or eta->pi+ pi-pi0 (gamma(s))
  int nPip=0,nPim=0,nPi0=0,nGam=0,nEta=0,nRho0=0;
  while ( 0 != (Cand = XuDaughter()) && IsBadDecayMode==false )   
    {
      PdtLund::LundType lundD=Cand->pdtEntry()->lundId();
      if(lundD==PdtLund::pi_plus) nPip+=1;
      else if(lundD==PdtLund::pi_minus) nPim+=1;
      else if(lundD==PdtLund::pi0) nPi0+=1;
      else if(lundD==PdtLund::eta) { nEta+=1; Eta = BtaCandidate(*Cand); }
      else if(lundD==PdtLund::rho0) nRho0+=1;
      else if(lundD==PdtLund::gamma) nGam+=1;
      else IsBadDecayMode=true; //once it's set to true, it can never be put back to false again
    }
  
  //Note: there's no ISR/FSR in Xu decays! (so we must ask nGam==0)
  bool etapipi=(nGam==0&&nPip==1&&nPim==1&&nPi0==0&&nEta==1&&nRho0==0);
  bool pipipi0=(nGam==0&&nPip==1&&nPim==1&&nPi0==1&&nEta==0&&nRho0==0);
  bool GammaGamma=(nGam==2&&nPip==0&&nPim==0&&nPi0==0&&nEta==0&&nRho0==0);
  bool rho0Gamma=(nGam==1&&nPip==0&&nPim==0&&nPi0==0&&nEta==0&&nRho0==1);

  //omega->pi+ pi- pi0 only!
  if(XuLund==223 && !pipipi0) IsBadDecayMode=true;

  //eta->gamma gamma and eta->pi+ pi- pi0 both accepted!
  if(XuLund==221 && !(pipipi0||GammaGamma) ) IsBadDecayMode=true;

  //eta' -> eta pi+ pi- and eta'->rho0 gamma both accepted
  if(XuLund==331 && !(etapipi||rho0Gamma) ) IsBadDecayMode=true;

  //Final test: eta'->eta pi+ pi-, the eta has to decay to Gam Gam or pi0 pi+ pi- ...
  if(XuLund==331 && etapipi ) 
    {
      HepAListIterator<BtaCandidate> EtaIter = Eta.daughterIterator();
      BtaCandidate* EtaD;
      int nPipEta=0,nPimEta=0,nPi0Eta=0,nGamEta=0;
      while ( 0 != (EtaD = EtaIter()) && IsBadDecayMode==false )   
	{
	  PdtLund::LundType lunEtaD=EtaD->pdtEntry()->lundId();
	  if(lunEtaD==PdtLund::pi_plus) nPipEta+=1;
	  else if(lunEtaD==PdtLund::pi_minus) nPimEta+=1;
	  else if(lunEtaD==PdtLund::pi0) nPi0Eta+=1;
	  else if(lunEtaD==PdtLund::gamma) nGamEta+=1;
	  else IsBadDecayMode=true; 
	}
      bool GG=(nGamEta==2&&nPimEta==0&&nPipEta==0&&nPi0Eta==0);
      bool PPP=(nGamEta==0&&nPimEta==1&&nPipEta==1&&nPi0Eta==1);
      if(!(GG||PPP)) IsBadDecayMode=true;
    }

  return IsBadDecayMode;
}

BtaCandidate
XSLMCTruthAnalyzer::GetXuDaughter(BtaCandidate XuLab, string mode)
{
  HepLorentzVector in = HepLorentzVector(0,0,0,0);
  BtaCandidate VDLab(in);

  //no merged pi0/eta in MCTruth...
  if(mode=="pilnu") { VDLab = BtaCandidate(XuLab); return VDLab; }

  //From Leif Wilden (/u/ec/wilden/vub/RELEASE/XslUser/Vub.cc):
  //rho0 -> take pi+, rho+/- -> take pi0, 
  //omega -> 3pi: boost the pi+ and the pi- in the omega rest frame, then take the Xproduct
  //omega -> anything else: take the first pion, whatever it is 
   

  //New suggestion (Dave Brown), we take the daughter with higher phi for pi0/eta...
  if(mode=="pi0lnu"||mode=="eta2lnu") 
    {
      bool f1Higher=(_fille1Lab.p3().phi()>_fille2Lab.p3().phi()) ? true:false;
      if(f1Higher) VDLab = BtaCandidate(_fille1Lab);
      else VDLab = BtaCandidate(_fille2Lab);
    }
  else if(mode=="etaplnuRG") VDLab = BtaCandidate(_fille1Lab);
  else if(mode=="rhoClnu") VDLab = BtaCandidate(_filleX0Lab);
  else if(mode=="rho0lnu") VDLab = BtaCandidate(_fille1Lab);
  else if(mode=="eta3lnu"||mode=="etaplnuE2PP"||mode=="etaplnuE3PP"||mode=="omegalnu")
    {
      HepLorentzVector piPlusHad(_fille1Lab.p4());
      HepLorentzVector piMinusHad(_fille2Lab.p4());
      piPlusHad.boost(-XuLab.p4().boostVector());
      piMinusHad.boost(-XuLab.p4().boostVector());
      HepLorentzVector p4VDaugXuFrame = HepLorentzVector(piPlusHad.vect().cross(piMinusHad.vect()),0.);
      
      //boost back in the lab frame for uniformity with other decay modes...
      HepLorentzVector p4VDaugLabFrame(p4VDaugXuFrame);
      p4VDaugLabFrame.boost(XuLab.p4().boostVector());
      
      //Create the BtaCandidate (with no type, but we don't care)
      VDLab = BtaCandidate(p4VDaugLabFrame);
    }
  else 
    {       
      //Possible for bad xu decay mode
      HepAListIterator<BtaCandidate> dauIter(XuLab.daughterIterator()); 
      VDLab = BtaCandidate(*dauIter()); 

      //Cette commande batarde veut dire: "si mode ne contient pas "_o"..."
      if( mode.find("_o",0)==string::npos ) cout<<"Unexpected case in XSLMCTruthAnalyzer::GetVD()."<<endl; 
    }
    
  //Voila! :-)
  return VDLab;
}

string
XSLMCTruthAnalyzer::FigureOutXuFamily(BtaCandidate Xu, int type)
{
  //Init...
  string mode="unknown";
  HepLorentzVector in = HepLorentzVector(0,0,0,0);
  _fille1Lab = BtaCandidate(in);
  _fille2Lab = BtaCandidate(in);
  _filleX0Lab = BtaCandidate(in);

  _pfille1Lab = BtaCandidate(in);
  _pfille2Lab = BtaCandidate(in);
  _pfillePi0Lab = BtaCandidate(in);

  _ppfilleGam1Lab = BtaCandidate(in);
  _ppfilleGam2Lab = BtaCandidate(in);

      
  if(type==1||type==11) { mode="pilnu"; return mode; }
  int nD=Xu.nDaughters();
  int nDD=0;

  BtaCandidate* Cand;
  //Filles
  //There is no merged pi0 in the MCTruth
  HepAListIterator<BtaCandidate> Xu_Daug = Xu.daughterIterator();
  while ( 0 != (Cand = Xu_Daug()) )   
    {
      PdtLund::LundType lund = Cand->pdtEntry()->lundId();

      if(lund==PdtLund::pi_plus) _fille1Lab = BtaCandidate(*Cand); 
      else if(lund==PdtLund::pi_minus) _fille2Lab = BtaCandidate(*Cand); 
      else if(lund==PdtLund::pi0) _filleX0Lab = BtaCandidate(*Cand); 
      else if(lund==PdtLund::gamma && _fille1Lab.p()==0 ) _fille1Lab = BtaCandidate(*Cand); 
      else if(lund==PdtLund::gamma ) _fille2Lab = BtaCandidate(*Cand); 
      else if(lund==PdtLund::eta ) { _filleX0Lab = BtaCandidate(*Cand); nDD=_filleX0Lab.nDaughters(); }
      else if(lund==PdtLund::rho0 ) _filleX0Lab = BtaCandidate(*Cand); 
      else if(type>0) cout<<"Problem!! Probably a bug in XSLMCTruthAnalyzer::FigureOutXuFamily..."<<endl;
    }

  if(type==2||type==22) { mode="pi0lnu"; return mode; } //No more family than the (one) two gamma(s) in this case
  if((type==3||type==33) &&nD==2) { mode="eta2lnu"; return mode; }


  //Petites filles
  HepAListIterator<BtaCandidate> petiteFille = _filleX0Lab.daughterIterator();
  while ( 0 != (Cand = petiteFille()) )   
    {
      PdtLund::LundType lund = Cand->pdtEntry()->lundId();

      if(lund==PdtLund::pi_plus) _pfille1Lab = BtaCandidate(*Cand); 
      else if(lund==PdtLund::pi_minus) _pfille2Lab = BtaCandidate(*Cand); 
      else if(lund==PdtLund::pi0) _pfillePi0Lab = BtaCandidate(*Cand); 
      else if(lund==PdtLund::gamma && _pfille1Lab.p()==0 ) _pfille1Lab = BtaCandidate(*Cand); 
      else if(lund==PdtLund::gamma ) _pfille2Lab = BtaCandidate(*Cand); 
    }

  if((type==3||type==33) &&nD==3) { mode="eta3lnu"; return mode; }
  if((type==4||type==44) &&nD==2) { mode="etaplnuRG"; return mode; }
  if((type==4||type==44) &&nD==3 &&nDD==2) { mode="etaplnuE2PP"; return mode; }
  if(type==5||type==55) { mode="rhoClnu"; return mode; } 
  if(type==6||type==66) { mode="rho0lnu"; return mode; } 
  if(type==7||type==77) { mode="omegalnu"; return mode; } 
  
  
  //Toutes petites filles 
  //Le seul cas ou cela est possible est: eta'->eta pi+ pi-, eta->pi+ pi- pi0, pi0->gamma gamma
  HepAListIterator<BtaCandidate> toutePetiteFille = _pfillePi0Lab.daughterIterator();
  while ( 0 != (Cand = toutePetiteFille()) )   
    {
      PdtLund::LundType lund = Cand->pdtEntry()->lundId();

      if(lund==PdtLund::gamma && _ppfilleGam1Lab.p()==0 ) _ppfilleGam1Lab = BtaCandidate(*Cand); 
      else if(lund==PdtLund::gamma ) _ppfilleGam2Lab = BtaCandidate(*Cand); 
    }

  if((type==4||type==44) &&nD==3 &&nDD==3) { mode="etaplnuE3PP"; return mode; }

  //Last case: Bad Xu decay mode...
  if(type==-2||type==-22) { mode="pi0lnu_o"; return mode; }
  if(type==-3||type==-33) { mode="etalnu_o"; return mode; }
  if(type==-4||type==-44) { mode="etaplnu_o"; return mode; }
  if(type==-7||type==-77) { mode="omegalnu_o"; return mode; }
  
  //else
  cout<<"Problems in XSLMCTruthAnalyzer::FigureOutXuFamily!!!"<<endl;

  return mode;

}


void 
XSLMCTruthAnalyzer::SetB1Family()
{
  if(_typeB2<0||_typeB2>77) return; //nothing to do in this case
  else if(_typeB1>=1&&_typeB1<=77) string tmp = FigureOutXuFamily(_Xu1_Lab, _typeB1); //this will replace the Xu family from Xu2 to Xu1
  else cout<<"SetB1Family:: No family to set!!"<<endl;
  return;
}

void 
XSLMCTruthAnalyzer::SetB2Family()
{
  if(_typeB1<0||_typeB1>77) return; //nothing to do in this case
  else if(_typeB2>=1&&_typeB2<=77) string tmp = FigureOutXuFamily(_Xu2_Lab, _typeB2); //this could replace the Xu family from Xu1 to Xu2
                                                                               //but is normally set to Xu2 by default
  else cout<<"SetB2Family:: No family to set!!"<<endl;
  return;
}

void 
XSLMCTruthAnalyzer::AnalyzeMissingMomentum()
{  
  //We assume that everything is already built
  //We consider any neutrino coming from a B, a D or a Tau (the signal nu has to substracted by the caller if wanted)
  //We consider all Klongs

  //Extra and missing trk/neutrals are done somewhere else (in the BetaMiniApp for example)

  _nNus=0;
  _nKLongs=0; 
  _nNeutrons=0; 
  _Nus=HepLorentzVector(0,0,0,0);
  _KLongs=HepLorentzVector(0,0,0,0);
  _Neutrons=HepLorentzVector(0,0,0,0);


  BtaCandidate* TruthCand;
  HepAListIterator<BtaCandidate> iterMC(_mcList);
  while ( 0 != (TruthCand = iterMC()) )         
    {
      int lund = abs( (int)TruthCand->pdtEntry()->lundId() );
      if(lund==12||lund==14||lund==16) //these are neutrinos
	{
	  int MotherLund = abs( (int)TruthCand->theMother()->pdtEntry()->lundId() );
	  int quark = MotherLund%1000;

	  if(quark>=400 || MotherLund==15) { _nNus+=1; _Nus+=TruthCand->p4(); }
	}
      
      else if(lund==PdtLund::K_L0) { _nKLongs+=1; _KLongs+=TruthCand->p4(); }
      else if(lund==PdtLund::n0) { _nNeutrons+=1; _Neutrons+=TruthCand->p4(); }
  
    }// while ( 0 != (TruthCand = iterMC()) )  
  
  return;
}










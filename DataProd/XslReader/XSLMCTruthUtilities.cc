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

//XSL code
#include "XslTools/XSLRecoYAnalyzer.hh"
#include "XslTools/XSLMCTruthAnalyzer.hh"

using std::cout;
using std::endl;

///////////////

float
XSLReader::IsSignal( XSLRecoYAnalyzer* recoY, XSLMCTruthAnalyzer* EvtTruth )
{
  float code=0;

  /*
    //A ajouter: 1 good signal XX from Decay Tree
                 2 true XX (but not the one from Decay Tree)
		 3 bad XX from direct truth match
		 4 bad XX from daughters
		 5 bad XX from both methods 
		 (or some other numbers but that's the idea...)

    code values: 
    1, good signal with electron
    2, good signal with muon
    -xxxxxxxxx1 good electron, combanitorial background
    -xxxxxxxxx2 bad electron, combanitorial background
    -xxxxxxxxx3 good muon, combanitorial background
    -xxxxxxxxx4 bad muon, combanitorial background
    -xxxxxxxxx6 good electron,regular background
    -xxxxxxxxx7 bad electron, regular background
    -xxxxxxxxx8 good muon, regular background
    -xxxxxxxxx9 bad muon, regular background
    
    -xxxxxxxx1x good Xu
    -xxxxxxxx2x bad Xu from direct Truth Match
    -xxxxxxxx3x bad Xu from Daughters
    -xxxxxxxx4x bad Xu from both methods

    -xxxxxxx1xx good fille1 (Gam1/Pip)
    -xxxxxxx2xx bad fille1 (Gam1/Pip)
    -xxxxxx1xxx good fille2 (Gam2 ou Pim)
    -xxxxxx2xxx bad fille2 (Gam2 ou Pim)
    -xxxxx1xxxx good filleX0 (Pi0/Eta/Rho0) 
    -xxxxx2xxxx bad filleX0 (Pi0/Eta/Rho0) from direct Truth Match
    -xxxxx3xxxx bad filleX0 (Pi0/Eta/Rho0) from Daughters
    -xxxxx4xxxx bad filleX0 (Pi0/Eta/Rho0) from both methods

    -xxxx1xxxxx good pfille1 (Gam1 ou Pip)
    -xxxx2xxxxx bad pfille1 (Gam1 ou Pip)
    -xxx1xxxxxx good pfille2 (Gam2 ou Pim)
    -xxx2xxxxxx bad pfille2 (Gam2 ou Pim)
    -xx1xxxxxxx good pfillePi0 
    -xx2xxxxxxx bad pfillePi0 from direct Truth Match
    -xx3xxxxxxx bad pfillePi0 from Daughters
    -xx4xxxxxxx bad pfillePi0 from both methods

    -x1xxxxxxxx good ppfille1 (Gam1 ou Pip)
    -x2xxxxxxxx bad ppfille1 (Gam1 ou Pip)
    -1xxxxxxxxx good ppfille2 (Gam2 ou Pim)
    -2xxxxxxxxx bad ppfille2 (Gam2 ou Pim)    
    
    in any case, 0 means irrelevant
  */
  
  //Some utilities...
  string recoMode=recoY->mode(); //of the B2 in case of 2 sig B
  bool mu=false,comb=false;
  int nDD=-1,nDDD=-1;
  float cons,dp,op;
  int nm;
  BtaCandidate lr,xur,f1r,f2r,f0r,pf1r,pf2r,pf0r,ppf1r,ppf2r;
  BtaCandidate lt,xut,f1t,f2t,f0t,pf1t,pf2t,pf0t,ppf1t,ppf2t;

  int lep=0,xu=0,f1=0,f2=0,fx0=0,pf1=0,pf2=0,pfx0=0,ppf1=0,ppf2=0;
  bool mergedL1=false,mergedL2=false;
  bool anyLund=false;
  bool noQualityCut=false;

  lr=recoY->LepLab();
  if(abs( (int)lr.pdtEntry()->lundId() )==PdtLund::mu_minus) mu=true; 

  xur=recoY->XuLab();
  double ch = xur.charge();
  PdtLund::LundType xurlun = xur.pdtEntry()->lundId();
  int nD=xur.nDaughters();
  bool merged = (recoMode=="pi0lnu"||recoMode=="eta2lnu")&&nD==0;

  if(recoMode==EvtTruth->modeB2()) 
    {
      comb=true;
      lt=EvtTruth->Lep2_Lab();
      xut=EvtTruth->Xu2_Lab();      
    }
  else if(recoMode==EvtTruth->modeB1()) 
    {
      comb=true;
      lt=EvtTruth->Lep1_Lab();
      xut=EvtTruth->Xu1_Lab();      
      EvtTruth->SetB1Family(); //this sets the Xu1 daughters (set to Xu2 by default)
    }

  if(comb) //sig candidate or combinatorial bkg
    {
      bool mcIsMotherComb=false;
      
      lep = ( TruthMatchTheseCand(&lr,&lt,cons,mcIsMotherComb,anyLund,noQualityCut) ) ? 1:2;
      if(mu) lep+=2; //3,4
      if(merged) xu = ( TruthMatchTheseMergedCand(&xur,&xut,mcIsMotherComb,anyLund,noQualityCut) ) ? 1:2; 
      else xu = ( TruthMatchTheseCand(&xur,&xut,cons,mcIsMotherComb,anyLund,noQualityCut) ) ? 1:2; 
      
      if(nD>0)
	{
	  //Some of these could be null but we take care of it later
	  f1r=recoY->fille1Lab();
	  f1t=EvtTruth->fille1Lab();
	  f2r=recoY->fille2Lab();
	  f2t=EvtTruth->fille2Lab();
	  f0r=recoY->filleX0Lab();
	  f0t=EvtTruth->filleX0Lab();	  

	  if(recoMode=="eta2lnu"||recoMode=="pi0lnu"||recoMode=="rho0lnu") //Xu->Gam1/pi+ Gam2/pi- 
	    {
	      nDD=0;
	      f1 = ( TruthMatchTheseCand(&f1r,&f1t,cons,mcIsMotherComb,anyLund,noQualityCut) ) ? 1:2;
	      f2 = ( TruthMatchTheseCand(&f2r,&f2t,cons,mcIsMotherComb,anyLund,noQualityCut) ) ? 1:2;
	      if(f1*f2==4)
		{
		  //Swap the truth gamma... it's true that there's no reason why Gam1Truth would be Gam1Rec and vice-versa
		  f1 = ( TruthMatchTheseCand(&f1r,&f2t,cons,mcIsMotherComb,anyLund,noQualityCut) ) ? 1:2;
		  f2 = ( TruthMatchTheseCand(&f2r,&f1t,cons,mcIsMotherComb,anyLund,noQualityCut) ) ? 1:2;		  
		}
	    }
	  else if(recoMode=="eta3lnu"||recoMode=="omegalnu"||recoMode=="etaplnuE2PP"||recoMode=="etaplnuE3PP") //Xu->pi+ pi- pi0/eta
	    {
	      nDD=recoY->filleX0Lab().nDaughters();
	      mergedL1 = (nDD==0) ? true:false;

	      f1 = ( TruthMatchTheseCand(&f1r,&f1t,cons,mcIsMotherComb,anyLund,noQualityCut) ) ? 1:2;
	      f2 = ( TruthMatchTheseCand(&f2r,&f2t,cons,mcIsMotherComb,anyLund,noQualityCut) ) ? 1:2;

	      if(mergedL1) fx0 = ( TruthMatchTheseMergedCand(&f0r,&f0t,mcIsMotherComb,anyLund,noQualityCut) ) ? 1:2;
	      else fx0 = ( TruthMatchTheseCand(&f0r,&f0t,cons,mcIsMotherComb,anyLund,noQualityCut) ) ? 1:2; 
	    }
	  else if(recoMode=="etaplnuRG") 
	    {
	      nDD=recoY->filleX0Lab().nDaughters();
	      f1 = ( TruthMatchTheseCand(&f0r,&f0t,cons,mcIsMotherComb,anyLund,noQualityCut) ) ? 1:2;
	      fx0 = ( TruthMatchTheseCand(&f0r,&f0t,cons,mcIsMotherComb,anyLund,noQualityCut) ) ? 1:2;	      
	    }
	  else if(recoMode=="rhoClnu") 
	    {
	      nDD=recoY->filleX0Lab().nDaughters();
	      mergedL1 = (nDD==0) ? true:false;

	      if(ch>0) f1 = ( TruthMatchTheseCand(&f1r,&f1t,cons,mcIsMotherComb,anyLund,noQualityCut) ) ? 1:2;
	      else f2 = ( TruthMatchTheseCand(&f2r,&f2t,cons,mcIsMotherComb,anyLund,noQualityCut) ) ? 1:2;
	      if(mergedL1) fx0=( TruthMatchTheseMergedCand(&f0r,&f0t,mcIsMotherComb,anyLund,noQualityCut) ) ? 1:2;
	      else fx0 = ( TruthMatchTheseCand(&f0r,&f0t,cons,mcIsMotherComb,anyLund,noQualityCut) ) ? 1:2;  
	    }
	  
	  if(nDD>0)
	    {	      
	      pf1r=recoY->pfille1Lab();
	      pf1t=EvtTruth->pfille1Lab();
	      pf2r=recoY->pfille2Lab();
	      pf2t=EvtTruth->pfille2Lab();

	      pf1 = ( TruthMatchTheseCand(&pf1r,&pf1t,cons,mcIsMotherComb,anyLund,noQualityCut) ) ? 1:2;
	      pf2 = ( TruthMatchTheseCand(&pf2r,&pf2t,cons,mcIsMotherComb,anyLund,noQualityCut) ) ? 1:2;

	      if(pf1t.charge()==0&&pf2t.charge()==0&&pf1*pf2==4)
		{
		  //Swap the truth gamma... it's true that there's no reason why Gam1Truth would be Gam1Rec and vice-versa
		  pf1 = ( TruthMatchTheseCand(&pf1r,&pf2t,cons,mcIsMotherComb,anyLund,noQualityCut) ) ? 1:2;
		  pf2 = ( TruthMatchTheseCand(&pf2r,&pf1t,cons,mcIsMotherComb,anyLund,noQualityCut) ) ? 1:2;
		}

	      if(recoMode=="etaplnuE3PP")
		{
		  pf0r=recoY->pfillePi0Lab();
		  pf0t=EvtTruth->pfillePi0Lab();	  

		  nDDD=pf0r.nDaughters();
		  mergedL2 = (nDDD==0) ? true:false;

		  if(mergedL2)
		    {pfx0=(TruthMatchTheseMergedCand(&pf0r,&pf0t,mcIsMotherComb,anyLund,noQualityCut) ) ? 1:2;}
		  else     
		    {
		      ppf1r=recoY->ppfilleGam1Lab();
		      ppf1t=EvtTruth->ppfilleGam1Lab();
		      ppf2r=recoY->ppfilleGam2Lab();
		      ppf2t=EvtTruth->ppfilleGam2Lab();

		      pfx0=(TruthMatchTheseCand(&pf0r,&pf0t,cons,mcIsMotherComb,anyLund,noQualityCut) ) ? 1:2;

		      ppf1=(TruthMatchTheseCand(&ppf1r,&ppf1t,cons,mcIsMotherComb,anyLund,noQualityCut) ) ? 1:2;
		      ppf2=(TruthMatchTheseCand(&ppf2r,&ppf2t,cons,mcIsMotherComb,anyLund,noQualityCut) ) ? 1:2;
		      if(ppf1*ppf2==4)
			{
			  //Swap the truth gamma... it's true that there's no reason why Gam1Truth would be Gam1Rec and vice-versa
			  ppf1=(TruthMatchTheseCand(&ppf1r,&ppf2t,cons,mcIsMotherComb,anyLund,noQualityCut) ) ? 1:2;
			  ppf2=(TruthMatchTheseCand(&ppf2r,&ppf1t,cons,mcIsMotherComb,anyLund,noQualityCut) ) ? 1:2;
			}		      
		    }
		}//if(recoMode=="etaplnuE3PP")	      
	    }//if(nDD>0)
	}//if(nD>0)
    }//if combinatorial bkg


  else //regular bkg - could not possibly be signal
    {
      BtaCandidate *lt=0,*xut=0,*f1t=0,*f2t=0,*fx0t=0,*pf1t=0,*pf2t=0,*pfx0t=0,*ppf1t=0,*ppf2t=0;
      bool mcIsMotherReg=true;

      lt = TruthMatch(&lr,cons,op,dp,nm,anyLund,noQualityCut); 
      lep = (lt!=0) ? 6:7;
      if(mu) lep+=2; //8,9
      if(merged) xut = TruthMatchMerged(&xur,op,dp,nm,anyLund,noQualityCut); 
      else xut = TruthMatch(&xur,cons,op,dp,nm,anyLund,noQualityCut); 
      xu = (xut!=0) ? 1:2;

      if(nD>0)
	{
	  //Some of these could be null but we take care of it later
	  f1r=recoY->fille1Lab();
	  f2r=recoY->fille2Lab();
	  f0r=recoY->filleX0Lab();
	  PdtLund::LundType f0rlun = PdtLund::null;

	  if(recoMode=="eta2lnu"||recoMode=="pi0lnu"||recoMode=="rho0lnu") //Xu->Gam1/pi+ Gam2/pi- 
	    {
	      nDD=0;
	      f1t = TruthMatch(&f1r,cons,op,dp,nm,anyLund,noQualityCut);
	      f1 = (f1t!=0) ? 1:2;

	      f2t = TruthMatch(&f2r,cons,op,dp,nm,anyLund,noQualityCut);
	      f2 = (f2t!=0) ? 1:2;

	      //Xcheck of xu truth matching #1, step 2 and 3 are in case xut is null
	      if(f1*f2!=1) xu+=2;
	      else if(f1t->theMother()==0 || f2t->theMother()==0) xu+=2;
	      else if(f1t->theMother()->uid() != f2t->theMother()->uid() ) xu+=2;
	      else if(!anyLund && xurlun!=f2t->theMother()->pdtEntry()->lundId() ) xu+=2;		  
	    }
	  else if(recoMode=="eta3lnu"||recoMode=="omegalnu"||recoMode=="etaplnuE2PP"||recoMode=="etaplnuE3PP") //Xu->pi+ pi- pi0/eta
	    {
	      nDD=f0r.nDaughters();
	      mergedL1 = (nDD==0) ? true:false;
	      f0rlun = f0r.pdtEntry()->lundId();

	      f1t = TruthMatch(&f1r,cons,op,dp,nm,anyLund,noQualityCut);
	      f1 = (f1t!=0) ? 1:2;

	      f2t = TruthMatch(&f2r,cons,op,dp,nm,anyLund,noQualityCut);
	      f2 = (f2t!=0) ? 1:2;

	      if(mergedL1) fx0t = TruthMatchMerged(&f0r,op,dp,nm,anyLund,noQualityCut);	
	      else fx0t = TruthMatch(&f0r,cons,op,dp,nm,anyLund,noQualityCut);	      
	      fx0 = (fx0t!=0) ? 1:2;

	      //Xcheck of xu truth matching #1, step 2 and 3 are in case xut is null
	      if(f1*f2*fx0!=1) xu+=2;
	      else if(f1t->theMother()==0 || f2t->theMother()==0 || fx0t->theMother()==0 ) xu+=2;
	      else if(f1t->theMother()->uid() != f2t->theMother()->uid() ) xu+=2;
	      else if(f1t->theMother()->uid() != fx0t->theMother()->uid() ) xu+=2;
	      else if(!anyLund && xurlun!=f2t->theMother()->pdtEntry()->lundId() ) xu+=2;		  
	    }
	  else if(recoMode=="etaplnuRG")
	    {
	      nDD=f0r.nDaughters();
	      if(nDD!=2) cout<<"Problem in etaplnuRG IsSignal"<<endl;
	      f0rlun = f0r.pdtEntry()->lundId();

	      f1t = TruthMatch(&f1r,cons,op,dp,nm,anyLund,noQualityCut);
	      f1 = (f1t!=0) ? 1:2;

	      fx0t = TruthMatch(&f0r,cons,op,dp,nm,anyLund,noQualityCut);	      
	      fx0 = (fx0t!=0) ? 1:2;

	      //Xcheck of xu truth matching #1, step 2 and 3 are in case xut is null
	      if(f1*fx0!=1) xu+=2;
	      else if(f1t->theMother()==0 || fx0t->theMother()==0 ) xu+=2;
	      else if(f1t->theMother()->uid() != fx0t->theMother()->uid() ) xu+=2;
	      else if(!anyLund && xurlun!=f1t->theMother()->pdtEntry()->lundId() ) xu+=2;		  
	    }
	  else if(recoMode=="rhoClnu")
	    {
	      nDD=f0r.nDaughters();
	      mergedL1 = (nDD==0) ? true:false;
	      f0rlun = f0r.pdtEntry()->lundId();

	      if(ch>0)
		{ f1t = TruthMatch(&f1r,cons,op,dp,nm,anyLund,noQualityCut);
		f1 = (f1t!=0) ? 1:2;}
	      else 
		{f2t = TruthMatch(&f2r,cons,op,dp,nm,anyLund,noQualityCut);
		f2 = (f2t!=0) ? 1:2;}

	      if(mergedL1) fx0t = TruthMatchMerged(&f0r,op,dp,nm,anyLund,noQualityCut);	
	      else fx0t = TruthMatch(&f0r,cons,op,dp,nm,anyLund,noQualityCut);	      
	      fx0 = (fx0t!=0) ? 1:2;
	      
	      //Xcheck of xu truth matching #1, step 2 and 3 are in case xut is null
	      if(ch>0)
		{
		  if(f1*fx0!=1) xu+=2;
		  else if(f1t->theMother()==0 || fx0t->theMother()==0 ) xu+=2;
		  else if(f1t->theMother()->uid() != fx0t->theMother()->uid() ) xu+=2;
		  else if(!anyLund && xurlun!=f1t->theMother()->pdtEntry()->lundId() ) xu+=2;		  
		}
	      else 
		{
		  if(f2*fx0!=1) xu+=2;
		  else if(f2t->theMother()==0 || fx0t->theMother()==0 ) xu+=2;
		  else if(f2t->theMother()->uid() != fx0t->theMother()->uid() ) xu+=2;
		  else if(!anyLund && xurlun!=f2t->theMother()->pdtEntry()->lundId() ) xu+=2;		  
		}
	    }//else if(recoMode=="rhoClnu")
	  
	  //final Xcheck of Xu truth-matching for all modes at the same time
	  if(xut!=0 && xu<3)
	    {
	      if(f1t!=0 && !TruthMatchTheseCand(&f1r,xut,cons,mcIsMotherReg,anyLund,noQualityCut) ) xu+=2;  
	      else if(f2t!=0 && !TruthMatchTheseCand(&f2r,xut,cons,mcIsMotherReg,anyLund,noQualityCut) ) xu+=2;  
	      else if(fx0t!=0 && !TruthMatchTheseCand(&f0r,xut,cons,mcIsMotherReg,anyLund,noQualityCut) ) xu+=2;  
	    }

	  if(nDD>0)
	    {
	      pf1r=recoY->pfille1Lab();
	      pf2r=recoY->pfille2Lab();

	      if(recoMode=="etaplnuE3PP")
		{
		  pf0r=recoY->pfillePi0Lab();
		  PdtLund::LundType pf0rlun = pf0r.pdtEntry()->lundId();
		  nDDD=pf0r.nDaughters();
		  mergedL2 = (nDDD==0) ? true:false;

		  pf1t = TruthMatch(&pf1r,cons,op,dp,nm,anyLund,noQualityCut);
		  pf1 = (pf1t!=0) ? 1:2;

		  pf2t = TruthMatch(&pf2r,cons,op,dp,nm,anyLund,noQualityCut);
		  pf2 = (pf2t!=0) ? 1:2;	      

		  if(mergedL2) pfx0t = TruthMatchMerged(&pf0r,op,dp,nm,anyLund,noQualityCut);	
		  else pfx0t = TruthMatch(&pf0r,cons,op,dp,nm,anyLund,noQualityCut);	      
		  pfx0 = (pfx0t!=0) ? 1:2;

		  //Xcheck of fX0 truth matching #1, step 2 and 3 are in case fx0t is null
		  if(pf1*pf2*pfx0!=1) fx0+=2;
		  else if(pf1t->theMother()==0 || pf2t->theMother()==0 || pfx0t->theMother()==0 ) fx0+=2;
		  else if(pf1t->theMother()->uid() != pf2t->theMother()->uid() ) fx0+=2;
		  else if(pf1t->theMother()->uid() != pfx0t->theMother()->uid() ) fx0+=2;
		  else if(!anyLund && f0rlun != pf2t->theMother()->pdtEntry()->lundId() ) fx0+=2;		  

		  if(!mergedL2)
		    {
		      ppf1r=recoY->ppfilleGam1Lab();
		      ppf2r=recoY->ppfilleGam2Lab();

		      ppf1t = TruthMatch(&ppf1r,cons,op,dp,nm,anyLund,noQualityCut);
		      ppf1 = (ppf1t!=0) ? 1:2;
		      ppf2t = TruthMatch(&ppf2r,cons,op,dp,nm,anyLund,noQualityCut);		      
		      ppf2 = (ppf2t!=0) ? 1:2;
		      
		      //Xcheck of pfX0 truth matching #1, step 2 and 3 are in case pfx0t is null
		      if(ppf1*ppf2!=1) pfx0+=2;
		      else if(ppf1t->theMother()==0 || ppf2t->theMother()==0 ) pfx0+=2;
		      else if(ppf1t->theMother()->uid() != ppf2t->theMother()->uid() ) pfx0+=2;
		      else if(!anyLund && pf0rlun != ppf2t->theMother()->pdtEntry()->lundId() ) pfx0+=2;
		      else if(pfx0t!=0 &&  pfx0<3) //final Xcheck
			{
			  if(ppf1t!=0 && !TruthMatchTheseCand(&ppf1r,pfx0t,cons,mcIsMotherReg,anyLund,noQualityCut) ) pfx0+=2;  
			  else if(ppf2t!=0 && !TruthMatchTheseCand(&ppf2r,pfx0t,cons,mcIsMotherReg,anyLund,noQualityCut) ) pfx0+=2;  
			}
		    }
		}//if(recoMode=="etaplnuE3PP")	      
	      else
		{
		  pf1t = TruthMatch(&pf1r,cons,op,dp,nm,anyLund,noQualityCut);
		  pf1 = (pf1t!=0) ? 1:2;
		  pf2t = TruthMatch(&pf2r,cons,op,dp,nm,anyLund,noQualityCut);
		  pf2 = (pf2t!=0) ? 1:2;	      

		  //Xcheck of fX0 truth matching #1, step 2 and 3 are in case fx0t is null
		  if(pf1*pf2!=1) fx0+=2;
		  else if(pf1t->theMother()==0 || pf2t->theMother()==0 ) fx0+=2;
		  else if(pf1t->theMother()->uid() != pf2t->theMother()->uid() ) fx0+=2;
		  else if(!anyLund && f0rlun != pf2t->theMother()->pdtEntry()->lundId() ) fx0+=2;		  
		}//else
	      
	      //final Xcheck of fx0 truth-matching for all modes at the same time
	      if(fx0t!=0 && fx0<3)
		{
		  if(pf1t!=0 && !TruthMatchTheseCand(&pf1r,fx0t,cons,mcIsMotherReg,anyLund,noQualityCut) ) fx0+=2;  
		  else if(pf2t!=0 && !TruthMatchTheseCand(&pf2r,fx0t,cons,mcIsMotherReg,anyLund,noQualityCut) ) fx0+=2;  
		  else if(pfx0t!=0 && !TruthMatchTheseCand(&pf0r,fx0t,cons,mcIsMotherReg,anyLund,noQualityCut) ) fx0+=2;  
		}	      

	    }//if(nDD>0)
	}//if(nD>0)

      if(lt!=0) delete lt;
      if(xut!=0) delete xut;
      if(f1t!=0) delete f1t;
      if(f2t!=0) delete f2t;
      if(fx0t!=0) delete fx0t;
      if(pf1t!=0) delete pf1t;
      if(pf2t!=0) delete pf2t;
      if(pfx0t!=0) delete pfx0t;
      if(ppf1t!=0) delete ppf1t;
      if(ppf2t!=0) delete ppf2t;
    }//else regular bkg

  //Summary of what we just learned

  code-=lep;
  code-=(xu*10);
  code-=(f1*100);
  code-=(f2*1000);
  code-=(fx0*10000);
  code-=(pf1*100000);
  code-=(pf2*1000000);
  code-=(pfx0*10000000);
  code-=(ppf1*100000000);
  code-=(ppf2*1000000000);

  //Signal?
  // la condition comb est necessaire: le Xu et le Lep peuvent etre bons, mais le B doit aussi avoir le bon decay
  if(comb&&xu<2&&f1<2&&f2<2&&fx0<2&&pf1<2&&pf2<2&&pfx0<2&&ppf1<2&&ppf2<2)
    { 
      if(lep==3)code=2; 
      else if(lep==1) code=1; 
    }
  

  //cross checks
  if(xu<2&&(lep==1||lep==3)&&!(f1<2&&f2<2&&fx0<2&&pf1<2&&pf2<2&&pfx0<2&&ppf1<2&&ppf2<2)) 
    { 
      cout<<"reco mode: "<<recoMode<<", Incompatible Truth matching #1!!    Xu&&Lep==true,   but NOT all the daughters are OK"<<endl;
      cout<<"  f1="<<f1<<"  f2="<<f2<<"  f0="<<fx0<<"  pf1="<<pf1<<"  pf2="<<pf2<<"  pf0="<<pfx0<<"  ppf1="<<ppf1<<"  ppf2="<<ppf2<<endl;
    }

  if(!comb&&xu>1&&xu!=4&&nD!=0) 
    { cout<<"reco mode: "<<recoMode<<", Different results for two _Xu_ truth-matching methods. xu="<<xu<<endl; }
  if(!comb&&fx0>1&&fx0!=4&&nDD!=0) 
    { cout<<"reco mode: "<<recoMode<<", Different results for two _fX0_ truth-matching methods. fx0="<<fx0<<endl; }
  if(!comb&&pfx0>1&&pfx0!=4) 
    { cout<<"reco mode: "<<recoMode<<", Different results for two _pfX0_ truth-matching methods. pfx0="<<pfx0<<endl; }

  return code;
}


int XSLReader::IsSignalOld( BtaCandidate* recoY )
{

  //Getting the recoDaughters...
  BtaCandidate* recoLep=0;
  BtaCandidate* recoHad=0;
  BtaCandidate* fille;
  HepAListIterator<BtaCandidate>Ydaugh = recoY->daughterIterator();
  int d=0;
  while ( 0 != (fille =Ydaugh()) )   
    {
      //The first daughter is always the hadron by construction
      //WARNING! this is the case only because the Y are created like this: BtaCandidate* Y = comb.create(*hadron,*lepton);
      //The hadron would however be the second daughter if Y would be created like: BtaCandidate* Y = comb.create(*lepton,*hadron);

      if(d==0) recoHad = new BtaCandidate( *fille );
      else if(d==1) recoLep = new BtaCandidate( *fille );
      else cout<<"PROBLEM in IsSignal!!"<<endl;
      d+=1;
    }

  int nM=0;
  float delp=-666,openAng=-666,cons=-666;
  BtaCandidate* MClep = TruthMatch(recoLep,cons,openAng,delp,nM);
  BtaCandidate* MChad = TruthMatch(recoHad,cons,openAng,delp,nM);
  //if(recoHad->pdtEntry()->lundId()==PdtLund::pi0 && recoHad->nDaughters()==0 ) TruthMatchMerged(recoHad,openAng,delp,nM);
  //else TruthMatch(recoHad,cons,openAng,delp,nM);

  //////////////
  int signal = -999;
  int nfilles=0;
  BtaCandidate tmp;
  bool isGoodNu = 0;
 
  if( MClep != 0 && MChad !=0 && MClep->theMother() !=0 && MChad->theMother() !=0)
    {
      int lundMChad = (int)MChad->pdtEntry()->lundId();
      int lundMClep = (int)MClep->pdtEntry()->lundId();
      if( 
	 ( abs( lundMChad ) == PdtLund::pi_plus
	   || abs( lundMChad ) == PdtLund::pi0 
	   || abs( lundMChad ) == PdtLund::eta 
	   || abs( lundMChad ) == PdtLund::eta_prime 
	   || abs( lundMChad ) == PdtLund::rho_plus 
	   || abs( lundMChad ) == PdtLund::rho0 
	   || abs( lundMChad ) == PdtLund::omega )
	 &&( abs( lundMClep ) == PdtLund::e_minus
	     || abs( lundMClep ) == PdtLund::mu_minus )
	 && MChad->theMother()->uid() == MClep->theMother()->uid()
	 && ( abs( (int)MChad->theMother()->pdtEntry()->lundId() ) == PdtLund::B0
	      || abs( (int)MChad->theMother()->pdtEntry()->lundId() ) == PdtLund::B_plus )
	 )
	{
	  // recherche du neutrino dans le MC
	  HepAListIterator<BtaCandidate> iter_fille = MClep->theMother()->daughterIterator();
	  BtaCandidate *candid_nu; 
	  while ( 0 != ( candid_nu = iter_fille() ) ) 
	    {
	      if (candid_nu->pdtEntry()->lundId()!= PdtLund::gamma) nfilles++;
	      if ( ((abs( (int)candid_nu->pdtEntry()->lundId() ) == PdtLund::nu_e) && 
		   (abs( lundMClep ) == PdtLund::e_minus))
		   || ((abs( (int)candid_nu->pdtEntry()->lundId() ) == PdtLund::nu_mu) &&
		   (abs( lundMClep ) == PdtLund::mu_minus))
		)
		{
		  isGoodNu = 1;
		  //if (!_FoundNu) _sigNuMC = new BtaCandidate(*candid_nu); //This is now done in BuildEventMCInfoFromTruthList()
		}
	    }
	  if(isGoodNu ==1)
	    {
	    if(nfilles==3)
	      {
		//_FoundNu = true;
		if( abs( lundMChad ) == PdtLund::pi_plus 
		    && abs( lundMClep ) == PdtLund::e_minus ) signal = 1; //B0->pi e nu 
		else if( abs( lundMChad ) == PdtLund::pi_plus 
			 && abs( lundMClep ) == PdtLund::mu_minus ) signal = 11; //B0->pi mu nu
 
		else if(abs( lundMChad ) == PdtLund::pi0 
			&& abs( lundMClep ) == PdtLund::e_minus ) signal = 2;//B->pi0 e nu
		else if(abs( lundMChad ) == PdtLund::pi0 
			&& abs( lundMClep ) == PdtLund::mu_minus ) signal = 22;//B->pi0 mu nu

		else if(abs( lundMChad ) == PdtLund::eta 
			&& abs( lundMClep ) == PdtLund::e_minus ) signal = 3;//B->eta e nu
		else if(abs( lundMChad ) == PdtLund::eta 
			&& abs( lundMClep ) == PdtLund::mu_minus ) signal = 33;//B->eta mu nu

		else if(abs( lundMChad ) == PdtLund::eta_prime 
			&& abs( lundMClep ) == PdtLund::e_minus ) signal = 4;//B->eta' e nu
		else if(abs( lundMChad ) == PdtLund::eta_prime 
			&& abs( lundMClep ) == PdtLund::mu_minus ) signal = 44;//B->eta' mu nu

		else if(abs( lundMChad ) == PdtLund::rho_plus 
			&& abs( lundMClep ) == PdtLund::e_minus ) signal = 5;//B0->rhoC e nu
		else if(abs( lundMChad ) == PdtLund::rho_plus 
			&& abs( lundMClep ) == PdtLund::mu_minus ) signal = 55;//B0->rhoC mu nu

		else if(abs( lundMChad ) == PdtLund::rho0 
			&& abs( lundMClep ) == PdtLund::e_minus ) signal = 6;//B->rho0 e nu
		else if(abs( lundMChad ) == PdtLund::rho0 
			&& abs( lundMClep ) == PdtLund::mu_minus ) signal = 66;//B->rho0 mu nu

		else if(abs( lundMChad ) == PdtLund::omega 
			&& abs( lundMClep ) == PdtLund::e_minus ) signal = 7;//B->omega e nu
		else if(abs( lundMChad ) == PdtLund::omega 
			&& abs( lundMClep ) == PdtLund::mu_minus ) signal = 77;//B->omega mu nu

	      } //if(nfilles==3)
	    else signal = -1; // We have pion l nu + something else
	    } // if(isGoodNu ==1)
	  else signal = -2; //we have pion + l coming from same B but nu is not OK
	}//abs( MChad->pdtEntry etc...) 
      else signal = 0; // we don't have pion & lepton coming from the same B
    }//if(MClep!=0 etc...)
  else signal = -3; // we don't know if we have signal MC matching haven't worked
 
  if(MClep!=0) delete MClep;
  if(MChad!=0) delete MChad;
  if(recoLep!=0) delete recoLep;
  if(recoHad!=0) delete recoHad;

  return signal;

}

bool
XSLReader::IndicesOfTruthMatched( BtaCandidate* rec, int ind[4], float cons[4], int &nMatch)
{
  //Return the indices in _MCTruthList of up to three truth matched candidates
  //0,1,2 are for the three standard matches, [3] is for the 1st quality match only.
  //The returned boolean indicates if one of the matches met quality criteria

  nMatch=0;
  ind[0]=-66; ind[1]=-66; ind[2]=-66; ind[3]=-66;
  cons[0]=-0.66; cons[1]=-0.66; cons[2]=-0.66; cons[3]=-0.66; 
  bool QualityMatchFound=false;

  //First taking care of the merged pi0 case.
  PdtLund::LundType lundRec = rec->pdtEntry()->lundId();
  if( (lundRec==PdtLund::pi0 || lundRec==PdtLund::eta) && rec->nDaughters()==0 ) 
    { cout<<"TruthMatch doesn't deal with merged pi0's or eta's. Please use TruthMatchMerged instead."<<endl; return false; }

  BtaCandidate* recForTruthMatch=0;
  if((lundRec==PdtLund::e_minus || lundRec==PdtLund::e_plus) && rec->nDaughters()!=0)  //We have problems with BremRecovered electrons! 
    { 
      BtaCandidate* cand;
      int g=0;
      HepAListIterator<BtaCandidate>iter_fille = rec->daughterIterator();
      while ( 0 != (cand =iter_fille()) && g>=0 )   
	{ if(cand->charge()!=0)
	  { 
	    recForTruthMatch = new BtaCandidate( *cand ); 
	    g=-1; 
	  }
	}
    }
  else recForTruthMatch = new BtaCandidate( *rec );

  //////////////////////////////
  //Matching the rec cand...
  double pRec=rec->p();
  double recCharge = rec->charge();
  BtaCandidate *Tmp=0;
  
  float consTmp=-0.66;
  int i=0;
  while( i<3 )
    {
      Tmp=0;
      consTmp=-0.66;
      
      if(recForTruthMatch->nDaughters()!=0) 
	{
	  //We can't use which!=0 for truth matching composite candidates 
	  Tmp = truthMapWithCon->mcFromRecoWithCon(recForTruthMatch,consTmp);
	  i=3;
	}
      else Tmp = truthMapWithCon->mcFromRecoWithCon(recForTruthMatch,consTmp, i);
      
      if( Tmp!=0 ) 
	{	  
	  nMatch+=1; //succesful truth-match!
	  int zeInd=CandIndexInList(Tmp,_MCTruthList);
	  ind[i]=zeInd;
	  cons[i]=consTmp;

	  //Quality match?
	  double mcCharge= Tmp->charge();
	  double pMC= Tmp->p();
	  if(mcCharge==recCharge && fabs(pMC-pRec)<_dpMax && !QualityMatchFound ){ 
	    QualityMatchFound=true;
	    ind[3]=zeInd;    //0,1,2 are for the three standard matches, [3] is for the 1st quality match only.
	    cons[3]=consTmp;
	  }//quality match?
	}//if(Tmp!=0)
      else{ i=3; } //If the ith match didn't work, the i+1th won't match either!
      i+=1;
    }
  delete recForTruthMatch;
  return QualityMatchFound;
}



BtaCandidate*
XSLReader::TruthMatch( BtaCandidate* rec, float &Consistency, float &OpeningAngle, float &delP, int &nMatch, bool anyLund, bool noQualityCut )
{
  BtaCandidate *MC = 0;

  //Merged pi0's are not taken care by this method!
  //In case of several matches, only the first good one MC cand is returned but the nMatch is counted.
  nMatch=0;

  //First taking care of the merged pi0 case.
  PdtLund::LundType lundRec = rec->pdtEntry()->lundId();
  if( (lundRec==PdtLund::pi0 || lundRec==PdtLund::eta) && rec->nDaughters()==0 ) 
    { cout<<"TruthMatch doesn't deal with merged pi0's or eta's. Please use TruthMatchMerged instead."<<endl; return MC; }

  BtaCandidate* recForTruthMatch=NULL;
  if((lundRec==PdtLund::e_minus || lundRec==PdtLund::e_plus) && rec->nDaughters()!=0)  //We have problems with BremRecovered electrons! 
    { 
      BtaCandidate* cand;
      int g=0;
      HepAListIterator<BtaCandidate>iter_fille = rec->daughterIterator();
      while ( 0 != (cand =iter_fille()) && g>=0 )   
	{ if(cand->charge()!=0)
	  { 
	    recForTruthMatch = new BtaCandidate( *cand ); 
	    g=-1; 
	  }
	}
    }
  else recForTruthMatch = new BtaCandidate( *rec );

  //////////////////////////////
  //Matching the rec cand...
  int i =0;
  bool Found=false;
  double pRec=rec->p();
  double recCharge = rec->charge();
  Hep3Vector p3Rec = rec->p3();
  BtaCandidate *Tmp=0;
  float cons;

  while( i<3 )
    {
      Tmp=0;
      
      if(recForTruthMatch->nDaughters()!=0) 
	{
	  //We can't use which!=0 for truth matching composite candidates 
	  Tmp = truthMapWithCon->mcFromRecoWithCon(recForTruthMatch,cons);
	  i=3;
	}
      else Tmp = truthMapWithCon->mcFromRecoWithCon(recForTruthMatch,cons, i);
      
      if( Tmp!=0 ) 
	{	  
	  double ch= Tmp->charge();
	  double pMC= Tmp->p();
	  PdtLund::LundType lundMC = Tmp->pdtEntry()->lundId();
	  Hep3Vector p3MC = Tmp->p3();

	  if(_MyVerbose.value())
	    {
	      cout<<"i: "<<i<<endl;
	      cout<<"LundId: "<<(int)lundMC;
	      cout<<"pMC: "<<pMC<<endl;
	      cout<<"rec->p(): "<<pRec<<endl;
	    }
	
	  if( ((lundMC==lundRec)||anyLund) && ((ch==recCharge && fabs(pMC-pRec)<_dpMax ) || noQualityCut) ) 
	    { 
	      if(Found==false)
		{
		  Found=true; 
		  OpeningAngle = p3MC.angle(p3Rec);
		  delP = pMC - pRec;      
		  Consistency = cons;
		  MC = new BtaCandidate(*Tmp); 
		}
	      nMatch+=1;
	    }
	}
      i+=1;
    }
  delete recForTruthMatch;
  
  //deleted by the caller
  return MC;
}



bool
XSLReader::TruthMatchTheseCand( BtaCandidate* rec, BtaCandidate* mc, float &Consistency, bool mcIsMother, bool anyLund, bool noQualityCut )
{
  bool Found=false;
  //Merged pi0's are not taken care by this method!
  PdtLund::LundType lundRec = rec->pdtEntry()->lundId();
  if( (lundRec==PdtLund::pi0 || lundRec==PdtLund::eta) && rec->nDaughters()==0 ) 
    { cout<<"TruthMatchTheseCand doesn't deal with merged pi0's or eta's. Please use TruthMatchTheseMergedCand instead."<<endl; return Found; }

  BtaCandidate* recForTruthMatch=0;
  if((lundRec==PdtLund::e_minus || lundRec==PdtLund::e_plus) && rec->nDaughters()!=0)  //We have problems with BremRecovered electrons! 
    { 
      BtaCandidate* cand;
      int g=0;
      HepAListIterator<BtaCandidate>iter_fille = rec->daughterIterator();
      while ( 0 != (cand =iter_fille()) && g>=0 )   
	{ if(cand->charge()!=0)
	  { 
	    recForTruthMatch = new BtaCandidate( *cand ); 
	    g=-1; 
	  }
	}
    }
  else recForTruthMatch = new BtaCandidate( *rec );

  //////////////////////////////
  //Matching the rec cand...
  BtaCandidate *Tmp;
  int i =0;
  double ch=-999;
  double pMC=-999;
  double pRec=rec->p();
  float cons;

  while( i<3 && Found==false )
    {
      Tmp=0;

      if(recForTruthMatch->nDaughters()!=0) 
	{
	  //We can't use which!=0 for truth matching composite candidates 
	  Tmp = truthMapWithCon->mcFromRecoWithCon(recForTruthMatch,cons);
	  i=3;
	}
      else Tmp = truthMapWithCon->mcFromRecoWithCon( recForTruthMatch,cons, i);
      
      if( Tmp!=0 ) 
	{	  
	  ch= Tmp->charge();
	  pMC= Tmp->p();
	  PdtLund::LundType lundMC = Tmp->pdtEntry()->lundId();

	  if(_MyVerbose.value())
	    {
	      cout<<"i: "<<i<<endl;
	      cout<<"LundId: "<<(int)lundMC<<endl;
	      cout<<"pMC: "<<pMC<<endl;
	      cout<<"rec->p(): "<<pRec<<endl;
	    }
	
	  if( ((lundMC==lundRec)||anyLund) && ((ch==rec->charge() && fabs(pMC-pRec)<_dpMax ) || noQualityCut) ) 
	    { 
	      if(!mcIsMother&&(int)Tmp->uid()==(int)mc->uid())
		{
		  Found=true; 
		  Consistency = cons;
		}
	      else if(mcIsMother&&(int)Tmp->theMother()!=0) 
		{ if((int)Tmp->theMother()->uid()==(int)mc->uid())
		  {
		    Found=true; 
		    Consistency = cons;
		  }}
	    }
	}
      i+=1;
    }
  delete recForTruthMatch;
  
  return Found;
}


BtaCandidate*
XSLReader::TruthMatchMerged( BtaCandidate* rec, float &OpeningAngle, float &delP, int &nFound, bool anyLund, bool noQualityCut )
{
  //Warning: Option truth match to this lund implicit in merged pi0's truth-matching methods!
  BtaCandidate* MC=0;
  nFound=0;
  if(rec->nDaughters()!=0) 
    {
      cout<<"TruthMatchMerged is intended to be used for merged eta/pi0 only but you're trying to use if for a composite!!"<<endl;
      return MC;
    }
  int lundRec = (int)rec->pdtEntry()->lundId();
  if(lundRec!=PdtLund::pi0 && lundRec!=PdtLund::eta)
    {
      cout<<"TruthMatchMerged is intended to be used for merged eta/pi0 only!!"<<endl;
      cout<<"You are trying to use it for a ptcle with lund: "<<lundRec<<"..."<<endl;
      return MC;      
    }

  //////////////////////////////
  //First looking at the photons matched with the merged and their mother
  BtaCandidate *Tmp=0,*mother1=0,*mother2=0,*mother3=0;
  bool firstMotherNotUsed=true,secondMotherNotUsed=true,thirdMotherNotUsed=true;
  float Consistency;

  int i =0;
  int nGood=0;
  while( i<3 )
    {
      //Nothe that mcFromReco always returns a gamma for merged pi0's
      Tmp = truthMapWithCon->mcFromRecoWithCon( rec, Consistency, i);
      
      if( Tmp!=0 ){ if(Tmp->theMother()!=0)
	{
	  if((int)Tmp->theMother()->pdtEntry()->lundId()==lundRec || anyLund) 
	    {
	      if(firstMotherNotUsed){ mother1 = new BtaCandidate(*Tmp->theMother()); firstMotherNotUsed=false; nGood+=1; }
	      else if(secondMotherNotUsed){ mother2 = new BtaCandidate(*Tmp->theMother()); secondMotherNotUsed=false; nGood+=1; }
	      else if(thirdMotherNotUsed){ mother3 = new BtaCandidate(*Tmp->theMother()); thirdMotherNotUsed=false; nGood+=1; }
	    }
	}//if(Tmp->theMother()!=0)
      }//if(Tmp!=0)
      i+=1;
    }//while(i<3)
  if(nGood<2) 
    {
      if(mother1!=0) delete mother1;
      if(mother2!=0) delete mother2;
      if(mother3!=0) delete mother3;
      return MC;
    }

  float pRec=rec->p();
  bool notFound=true;
  if(mother1!=0 && mother2!=0)
    {
      if( (int)mother1->uid()==(int)mother2->uid() 
	  && ( (int)mother1->pdtEntry()->lundId()==lundRec || anyLund) 
	  && ( ( mother1->charge()==0 && fabs(mother1->p()-pRec)<_dpMax ) || noQualityCut) )
	{
	  nFound+=1;
	  if(notFound)
	    {
	      notFound=false; 
	      OpeningAngle = mother1->p3().angle(rec->p3());
	      delP = mother1->p() - pRec;      
	      MC = new BtaCandidate(*mother1); 	      
	    }
	}
    }
  if(mother1!=0 && mother3!=0)
    {
      if( (int)mother1->uid()==(int)mother3->uid() 
	  && ( (int)mother1->pdtEntry()->lundId()==lundRec || anyLund) 
	  && ( ( mother1->charge()==0 && fabs(mother1->p()-pRec)<_dpMax ) || noQualityCut) )
	{
	  nFound+=1;
	  if(notFound)
	    {
	      notFound=false; 
	      OpeningAngle = mother1->p3().angle(rec->p3());
	      delP = mother1->p() - pRec;      
	      MC = new BtaCandidate(*mother1); 	      
	    }
	}
    }
  if(mother2!=0 && mother3!=0)
    {
      if( (int)mother2->uid()==(int)mother3->uid() 
	  && ( (int)mother2->pdtEntry()->lundId()==lundRec || anyLund) 
	  && ( ( mother2->charge()==0 && fabs(mother2->p()-pRec)<_dpMax ) || noQualityCut) )
	{
	  nFound+=1;
	  if(notFound)
	    {
	      notFound=false; 
	      OpeningAngle = mother2->p3().angle(rec->p3());
	      delP = mother2->p() - pRec;      
	      MC = new BtaCandidate(*mother1); 	      
	    }
	}
    }

  if(mother1!=0) delete mother1;
  if(mother2!=0) delete mother2;
  if(mother3!=0) delete mother3;

  //deleted by the caller (if !=0)
  return MC;
}


bool
XSLReader::TruthMatchTheseMergedCand( BtaCandidate* rec, BtaCandidate* mc, bool mcIsMother, bool anyLund, bool noQualityCut )
{
  //Warning: Option truth match to this lund implicit in merged pi0's truth-matching methods!
  bool Found=false;
  if(rec->nDaughters()!=0) 
    {
      cout<<"TruthMatchMerged is intended to be used for merged eta/pi0 only but you're trying to use if for a composite!!"<<endl;
      return Found;
    }
  int recLund = (int)rec->pdtEntry()->lundId();
  if(recLund!=PdtLund::pi0 && recLund!=PdtLund::eta)
    {
      cout<<"TruthMatchMerged is intended to be used for merged eta/pi0 only!!"<<endl;
      cout<<"You are trying to use it for a ptcle with lund: "<<recLund<<"..."<<endl;
      return Found;      
    }
  
  //////////////////////////////
  //First looking at the photons matched with the merged and their mother
  BtaCandidate *Tmp=0,*mother1=0,*mother2=0,*mother3=0;
  bool firstMotherNotUsed=true,secondMotherNotUsed=true,thirdMotherNotUsed=true;
  float Consistency;
  
  int i =0;
  int nGood=0;
  while( i<3 )
    {
      //Nothe that mcFromReco always returns a gamma for merged pi0's
      Tmp = truthMapWithCon->mcFromRecoWithCon( rec, Consistency, i);
      
      if( Tmp!=0 ) { if(Tmp->theMother()!=0)
	{	  
	  if(Tmp->theMother()->pdtEntry()->lundId()==recLund || anyLund) 
	    {
	      if(firstMotherNotUsed){ mother1 = new BtaCandidate(*Tmp->theMother()); firstMotherNotUsed=false; nGood+=1; }
	      else if(secondMotherNotUsed){ mother2 = new BtaCandidate(*Tmp->theMother()); secondMotherNotUsed=false; nGood+=1; }
	      else if(thirdMotherNotUsed){ mother3 = new BtaCandidate(*Tmp->theMother()); thirdMotherNotUsed=false; nGood+=1; }
	    }
	}//if(Tmp->theMother()!=0)
      }//if(Tmp!=0)
      i+=1;
    }//while(i<3)
  if(nGood<2) 
    {
      if(mother1!=0) delete mother1;
      if(mother2!=0) delete mother2;
      if(mother3!=0) delete mother3;
      return Found;
    }

  float pRec=rec->p();
  if(mother1!=0 && mother2!=0)
    {
      if( (int)mother1->uid()==(int)mother2->uid() && ( (mother1->charge()==0 && fabs(mother1->p()-pRec)<_dpMax) || noQualityCut ) )
	{
	  if( (!mcIsMother&&(int)mother1->uid()==(int)mc->uid()) || (mcIsMother&&(int)mother1->theMother()->uid()==(int)mc->uid()))
	    {
	      Found=true; 
	    }
	}
    }
  if(mother1!=0 && mother3!=0)
    {
      if( (int)mother1->uid()==(int)mother3->uid() && ( (mother1->charge()==0 && fabs(mother1->p()-pRec)<_dpMax) || noQualityCut ) )
	{
	  if( (!mcIsMother&&(int)mother1->uid()==(int)mc->uid()) || (mcIsMother&&(int)mother1->theMother()->uid()==(int)mc->uid()))
	    {
	      Found=true; 
	    }
	}
    }
  if(mother2!=0 && mother3!=0)
    {
      if( (int)mother2->uid()==(int)mother3->uid() && ( (mother2->charge()==0 && fabs(mother2->p()-pRec)<_dpMax) || noQualityCut ) )
	{
	  if( (!mcIsMother&&(int)mother2->uid()==(int)mc->uid()) || (mcIsMother&&(int)mother2->theMother()->uid()==(int)mc->uid()))
	    {
	      Found=true; 
	    }
	}
    }
  
  if(mother1!=0) delete mother1;
  if(mother2!=0) delete mother2;
  if(mother3!=0) delete mother3;
  return Found;
}


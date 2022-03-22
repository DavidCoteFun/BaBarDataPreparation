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

#include "AbsEvent/AbsEvent.hh"
#include "XslTools/XSLMCTruthAnalyzer.hh"
#include "HepTuple/Tuple.h"

using std::cout;
using std::endl;

AppResult
XSLReader::event( AbsEvent* anEvent )
{

  _nbEvtProcessed+=1;
  bool cont=true;

  cont = InitializeEvent(anEvent);
     
  //Special quick and ugly hack for short terms tests...
  if(_readingMode.value()=="BypassUsrEventData") 
    {
      XSLMCTruthAnalyzer* EvtTruth;
      if( _lookMCTruth.value()==true )
	{ 
	  EvtTruth = FillNtupleWithMCTruthEventVariables( anEvent, _MainNtuple ); 
	  FillNtupleWithMCTruthBlockInfo( anEvent, _MainNtuple ); 
	  }
      else EvtTruth = new XSLMCTruthAnalyzer();

      _MainNtuple->dumpData();
      delete EvtTruth;

      return AppResult::OK;
    }

  if(cont) cont = GetEventInfoFromReco(anEvent);    
  if(cont) FillMainNtuple( anEvent );

  if(_goToTestZone.value()) 
    { 
      //TestZoneOne(anEvent); 
      TestZoneTwo(anEvent); 
    }

  return AppResult::OK;
}


int 
XSLReader::GetTrkIndex( BtaCandidate candid )
{  
  BtaCandidate target;
  PdtLund::LundType lundRec = candid.pdtEntry()->lundId();
  if( (lundRec==PdtLund::e_minus||lundRec==PdtLund::e_plus) && candid.nDaughters()!=0 )  //We have problems with BremRecovered electrons! 
    { 
      BtaCandidate* c;
      int g=0;
      HepAListIterator<BtaCandidate>iter_fille = candid.daughterIterator();
      while ( 0 != (c=iter_fille()) && g>=0 )   
	{ if(c->charge()!=0)
	  { 
	    target = BtaCandidate( *c ); 
	    g=-1; 
	  }
	}
    }
  else target = candid;


  int TrkNb=0;  
  BtaCandidate* trk;
  HepAListIterator<BtaCandidate> iter_trk(*BasicTrkList);  
  while ( 0 != ( trk = iter_trk()) ) 
    {
      if(trk->uid()==target.uid()) return TrkNb; 
      TrkNb+=1;
    }
  
  return -66;
}

int 
XSLReader::GetBremPhotonIndex( BtaCandidate lepton )
{  
  if(lepton.nDaughters()==0 ){ return -66; } //no Brem recovery for this guy!
  PdtLund::LundType lundRec = lepton.pdtEntry()->lundId();
  if( lundRec!=PdtLund::e_minus && lundRec!=PdtLund::e_plus ){ return -66; } 

  BtaCandidate photon;
  BtaCandidate *c;
  float pMax=-1.0;
  HepAListIterator<BtaCandidate>iter_fille = lepton.daughterIterator();
  while ( 0 != (c=iter_fille()) ) { 
    if(c->charge()==0.0 && c->p()>pMax)  //the candidate with hihest momentum is considered
      { 	
	photon = BtaCandidate( *c ); 
	pMax=c->p();
      }
  }
  
  if(pMax<0.0){ cout<<"We didn't find a neutral daughter of the electron!! This is unexpected!!"<<endl; return -66; }

  BtaCandidate* bump;
  int BumpNb=0;  
  HepAListIterator<BtaCandidate> iter_bump(*BasicBumpList);  
  while ( 0 != ( bump = iter_bump()) ) 
    {
      if(bump->uid()==photon.uid()) return BumpNb; 
      BumpNb+=1;
    }
  
  return -66;
}



int 
XSLReader::GetBumpIndex( BtaCandidate candid )
{  
  BtaCandidate* bump;

  int BumpNb=0;  
  HepAListIterator<BtaCandidate> iter_bump(*BasicBumpList);  
  while ( 0 != ( bump = iter_bump()) ) 
    {
      if(bump->uid()==candid.uid()) return BumpNb; 
      BumpNb+=1;
    }
  
  return -66;
}

bool
XSLReader::CandInList( BtaCandidate* cand, HepAList<BtaCandidate>* list )
{
  bool answer=false;
  
  BtaCandidate* c;
  HepAListIterator<BtaCandidate> iter(*list);  
  while ( 0 != ( c = iter()) && answer==false ) 
    { 
      if( c->uid()==cand->uid() )
	{answer=true;} 
    }

  return answer;
}

bool
XSLReader::CandIsDaughterInList( BtaCandidate* cand, HepAList<BtaCandidate>* list )
{
  bool answer=false;
  
  BtaCandidate *c,*cd;
  HepAListIterator<BtaCandidate> iter(*list);  
  while ( 0 != ( c = iter()) && answer==false ) 
    { 
      HepAListIterator<BtaCandidate>iter_fille = c->daughterIterator();
      while ( 0 != (cd=iter_fille()) && answer==false )   
	{
	  if( cd->uid()==cand->uid() )
	    {answer=true;} 
	}
    }

  return answer;
}

int
XSLReader::CandIndexInList( BtaCandidate* cand, HepAList<BtaCandidate>* list )
{
  bool answer=false;
  int index=0;

  BtaCandidate* c;
  HepAListIterator<BtaCandidate> iter(*list);  
  while ( 0 != ( c = iter()) && answer==false ) 
    { 
      if( c->uid()==cand->uid() )
	{answer=true;} 
      else index+=1;
    }

  if(answer) return index;
  else return -66;
}





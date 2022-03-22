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
//#include "XslReader/KLselector.hh"

#include "XslTools/XSLRotateAndBoost.hh"

#include "AbsEnv/AbsEnv.hh"
#include "BtaEnv/BtaEnv.hh"
#include "GenEnv/GenEnv.hh"
#include "AbsEvent/AbsEvent.hh"
#include "AbsEvent/getTmpAList.hh"
#include "ErrLogger/ErrLog.hh"
#include "UsrData/UsrEventBlock.hh"
#include "HepTuple/HTValOrderedVector.h"
#include "BetaCoreTools/BtaBooster.hh"
#include "Beta/EventInfo.hh"
#include "AbsEventTag/AbsEventTag.hh"
#include "BetaCoreTools/BtaThrust.hh"
#include "PDT/Pdt.hh"

#include "HepTuple/HTValOrderedVector.h"
#include "HepTuple/Tuple.h"


//For Timestamp
#include "AbsEvent/AbsEventID.hh"
#include "OdfCommon/odfTime.hh"
#include "EidData/EidEventTriplet.hh"
#include "EidData/EidCondKeyTriplet.hh"
#include "BdbTime/BdbTime.hh"

//Tracks
#include "TrkBase/TrkRecoTrk.hh"
#include "TrkBase/TrkFit.hh"
#include "TrkBase/TrkExchangePar.hh"
#include "TrkBase/TrkAbsFit.hh"
#include "TrkBase/TrkDifTraj.hh"
#include "TrkBase/TrkPoca.hh"
#include "BetaMicroAdapter/BtaMicroAdapter.hh"
#include "BetaMicroAdapter/BtaTrkQual.hh"
#include "BetaRecoAdapter/BtaAbsRecoObject.hh"
#include "BetaEvent/BtaParam.hh"

#include "EffTable/EffTable.hh"  //Efficiency tables from the DB! :-)
#include <map>
#include "BetaPid/PidWeight.hh"

//EMC bumps
#include "BetaMicroAdapter/BtaCalQual.hh"
#include "BetaMicroAdapter/BtaPidQual.hh"
#include "EmcData/EmcDigi.hh"
//#include "AbsCalo/AbsRecoCalo.hh"
//#include "EmcData/EmcCluster.hh"
//#include "EmcData/EmcXClMoments.hh"
//#include "EmcData/EmcClusterEnergySums.hh"
//#include "EmcData/EmcClusterMoments.hh"

//Ifr
#include "BetaMicroAdapter/BtaIfrQual.hh"
//#include "IfrPidData/IfrPidInfo.hh"

using std::cout;
using std::endl;

///////////////


bool 
XSLReader::GetEventInfoFromReco( AbsEvent* anEvent )
{
  ///////////////////////
  //  EVENT VARIABLES  //
  ///////////////////////

  AbsEventID* eventID = Ifd<AbsEventID>::get( anEvent , "AbsEventID" );
  if (eventID != 0) {
    _runNb =(int)eventID->run();
    odfTime timeStamp = eventID->eventIdTriplet().timeStamp();
    _platform = (int)eventID->eventIdTriplet().platform();  //the real type is EidPlatform_t...
    _partition = (int)eventID->eventIdTriplet().partitionMask(); //the real type is EidPlatform_t...
    _upperID = (int)timeStamp.binary().upper; //the real type is d_ULong...
    _lowerID = (int)timeStamp.binary().lower; //the real type is d_ULong...
    
    _modulus = timeStamp.modulus();

    string timeString;
    if ( _lookMCTruth.value()==true ) 
      {
	HepAList< AbsEventID >* eidList( 0 );
	getTmpAList( anEvent , eidList , IfdStrKey("AbsEventIDList") );      
	BdbTime bkgTime = ((*eidList)[ 1 ])->condKeyTriplet().key();
	timeString = bkgTime.asString("%Y%m%d", BdbTime::Local);
      }
    else
      {
	BdbTime fgTime = eventID->condKeyTriplet().key();
	timeString = fgTime.asString("%Y%m%d", BdbTime::Local);
      }
    
    
    _date = atoi(timeString.c_str());
   
  } 


  //float eMaxTight=0,muMaxTight=0,e2ndTight=0,mu2ndTight=0; // Check later if useful -->no! ;-)

  AbsEventTag* tag = Ifd<AbsEventTag>::get( anEvent );
  if(tag==0) 
    { 
      cout<<"NOTHING IN AbsEvent (sph, r2, etc)"<<endl;
    }
  else
    {
      if(tag->getFloat(_sphericityall,"sphericityAll"));
      if(tag->getFloat(_r2,"R2"));
      if(tag->getFloat(_r2all,"R2All"));
      //bool ismultihad;
      //if(tag->getBool(ismultihad,"BGFMultiHadron"))	{
      //	if(ismultihad) _ismultihadint=1;
      //else _ismultihadint=0;
      //}
      //bool BCMultiHadron;
      //if(tag->getBool(BCMultiHadron,"isBCMultiHadron"))	{
      //	if(BCMultiHadron) _isBCMH=1;
      //else _isBCMH=0;
      //}
      //bool L3EMC=tag->getBool(L3EMC,"L3OutEmc");
      //_L3OutEmc = (L3EMC) ? 1:0;
      //bool L3DCH=tag->getBool(L3DCH,"L3OutDch");
      //_L3OutDch = (L3DCH) ? 1:0;

      //if(tag->getFloat(eMaxTight,"elecTight1cm")); // Check later if useful
      //if(tag->getFloat(e2ndTight,"elecTight2cm")); // Check later if useful
      //if(tag->getFloat(muMaxTight,"muonTight1cm")); // Check later if useful
      //if(tag->getFloat(mu2ndTight,"muonTight2cm")); // Check later if usefu
    }
 
  /*
    //This is not useful for practical purposes... We need to it by hand with proper tweaking/weighting anyways.
    //A negative value of those values means a negative charge!// Check later if useful
  if (fabs(eMaxTight) > fabs(muMaxTight) )
    {
      _firstLepMax=eMaxTight;
      _secondLepMax = ( fabs(e2ndTight) > fabs(muMaxTight) ) ? e2ndTight : muMaxTight;
    }
  else
    {
      _firstLepMax=muMaxTight;
      _secondLepMax = ( fabs(mu2ndTight) > fabs(eMaxTight) ) ? mu2ndTight : eMaxTight;
    }// Check later if useful
  */


  HepAList< EventInfo >* infoList = 
    Ifd<HepAList< EventInfo > >::get(anEvent, _eventInfoList.value());
  if (infoList == 0){ ErrMsg(fatal) << "Could not locate event info list of name ";}
  _eventInfo = infoList->first();
  
  _UpsLab = gblEnv->getPep()->pepBeams()->total4Momentum();
  _s = _UpsLab.mag2(); 
  _EbeamCM = sqrt(_s)/2.; 

  _beamSpot= _eventInfo->beamSpot();


  //Thrust
  BtaThrust TH( *BasicTrkList, *_eventInfo);		    
  _Thrust_Axis = TH.thrust_axis();
  _Thrust_Mag = TH.thrust();


  //////////////////////////////////////
  //  SECTION EVENT INFO FROM TRACKS  //
  //////////////////////////////////////
  _PtotLab=HepLorentzVector(0,0,0,0);

  BtaCandidate* candid;
  HepAListIterator<BtaCandidate> iter_trk(*BasicTrkList);
  while ( 0 != ( candid = iter_trk()) ) {  _PtotLab += candid->p4();   }
 
  // 2. Photons
  HepAListIterator<BtaCandidate> iter_neut(*BasicBumpList);
  while ( 0 != ( candid = iter_neut() ) )  {  _PtotLab += candid->p4();   }
 

  /////////////
  //  Pmiss  //
  /////////////

  _PmissLab = _UpsLab - _PtotLab;

  XSLRotateAndBoost rAndB = XSLRotateAndBoost();
  _PmissUps = rAndB.BoostToFrame( _PmissLab, _UpsLab );
  _PtotUps = rAndB.BoostToFrame( _PtotLab, _UpsLab );


  ///////////
  //  EMC  //
  ///////////
  HepAList<EmcDigi>* digiList = Ifd<HepAList<EmcDigi> >::get(anEvent, _EmcDigiList.value());
  if (digiList != NULL ) {
    EmcDigi *aDigi;
    HepAListIterator<EmcDigi> iterEMC(*digiList);
    
    _eTotCal=0;
    while( 0 != ( aDigi = iterEMC() ) )
      {
	_eTotCal += aDigi->energy();
      }
  }
  else cout<<"couldn't find EmcDigiList"<<endl;


  return true;
}


void
XSLReader::FillNtupleWithEventVariables( HepTuple* ntuple )
{
  //timestamp
  ntuple->column("runNb", _runNb );
  ntuple->column("BunchNb",  _modulus );
  ntuple->column("date", (float)_date );
  ntuple->column("platform", _platform );
  ntuple->column("partition", _partition );
  ntuple->column("upperID", (float)_upperID );
  ntuple->column("lowerID", (float)_lowerID );
  ntuple->column("s_fromDB", (float)_s );
  ntuple->column("UpsLab_x", (float)_UpsLab.x() );
  ntuple->column("UpsLab_y", (float)_UpsLab.y() );
  ntuple->column("UpsLab_z", (float)_UpsLab.z() );

  //Event shape
  ntuple->column("R2",_r2);
  ntuple->column("R2all",_r2all);
  ntuple->column("Sphericity",_sphericityall);
  ntuple->column( "TH_mag", (float)_Thrust_Mag );
  //ntuple->column( "TH_x", (float)_Thrust_Axis.x() );
  //ntuple->column( "TH_y", (float)_Thrust_Axis.y() );
  //ntuple->column( "TH_z", (float)_Thrust_Axis.z() );

  //tag bits
  //ntuple->column("isMultiHad",_ismultihadint);
  //ntuple->column("isBCMH", _isBCMH);

  //Trigger L3
  //ntuple->column("L3OutEmc",_L3OutEmc);
  //ntuple->column("L3OutDch",_L3OutDch);

  //leptons
  //ntuple->column("firstLepMaxTightCM",(float)_firstLepMax);
  //ntuple->column("secondLepMaxTightCM",(float)_secondLepMax);

  //EMC
  ntuple->column("eTotCal",_eTotCal);

  //Missing/Visible energy
  ntuple->column("mm2", (float)_PmissLab.m2() );
  //ntuple->column("eMissLab",(float)_PmissLab.e() );
  //ntuple->column("phiMissLab", (float)_PmissLab.phi() );
  //ntuple->column("thMissLab", (float)_PmissLab.theta() );
  //ntuple->column("pMissLab", (float)_PmissLab.vect().mag() );
  ntuple->column("eTotLab",  (float)_PtotLab.e() );
  ntuple->column("pxTotLab", (float)_PtotLab.vect().x() );
  ntuple->column("pyTotLab", (float)_PtotLab.vect().y() );
  ntuple->column("pzTotLab", (float)_PtotLab.vect().z() );
  //ntuple->column("eMissUps",(float)_PmissUps.e() );
  //ntuple->column("phiMissUps", (float)_PmissUps.phi() );
  //ntuple->column("thMissUps", (float)_PmissUps.theta() );
  //ntuple->column("pMissUps", (float)_PmissUps.vect().mag() );
  ntuple->column("eTotUps",  (float)_PtotUps.e() );
  ntuple->column("pxTotUps", (float)_PtotUps.vect().x() );
  ntuple->column("pyTotUps", (float)_PtotUps.vect().y() );
  ntuple->column("pzTotUps", (float)_PtotUps.vect().z() );

  return;
}


void
XSLReader::FillNtupleWithIfrListBlock( AbsEvent* anEvent, HepTuple* ntuple )
{
  HTValOrderedVector<float> multip,likeType,likeValue,thLab,phiLab,eLab;
  HTValOrderedVector<int> firstHit,lastHit,uid,nMatch,TruInd;
  
  //Most infos are copyed from BtaGoodNeutralSelector
  BtaCandidate* candid;
  HepAListIterator<BtaCandidate> iter_ifr(*BasicIfrList);  
  
  while ( 0 != ( candid = iter_ifr()) ) {
    
    const BtaIfrQual* IfrQual = candid->getMicroAdapter()->getIfrQual(); 

    /*
    //attempt to migrate the Quals out of mu code... Abandonned for the moment!
    cout<<"IfrQual->firstHit(): "<<IfrQual->firstHit()<<endl;
    if( candid->pidInfoSummary()!=0 ){
      cout<<"allo"<<endl;
      if( candid->pidInfoSummary()->ifrPidInfo()!=0 ){
	float first = candid->pidInfoSummary()->ifrPidInfo()->firstLayer();
	cout<<"first: "<<first<<endl;
    
	cout<<"IfrQual->lastHit(): "<<IfrQual->lastHit()<<endl;
	float last = candid->pidInfoSummary()->ifrPidInfo()->lastLayer();
	cout<<"last: "<<last<<endl;
    
	cout<<"IfrQual->IfrNStrips(): "<<IfrQual->IfrNStrips()<<endl;
	float strips = candid->pidInfoSummary()->ifrPidInfo()->numberOfStrips();
	cout<<"strips: "<<strips<<endl;    

	cout<<"IfrQual->IfrLayHits(): "<<IfrQual->IfrLayHits()<<endl;
	float layHits = candid->pidInfoSummary()->ifrPidInfo()->numberOfHits();
	cout<<"layHits: "<<layHits<<endl;    
      }
    }
    */

    //uid.append( (int)candid->uid() );
    thLab.append(candid->p3().theta());
    phiLab.append(candid->p3().phi());
    //eLab.append(candid->energy());
    firstHit.append(IfrQual->firstHit());
    lastHit.append(IfrQual->lastHit());    
    float mul = float(IfrQual->IfrNStrips())/float(IfrQual->IfrLayHits());
    multip.append( mul );

    /*
    double likelihood[2];
    int flag = 2; //flag=1 : EMC cand, flag=2 : IFR cand 
    KLselector::Selector(anEvent, _runNb, candid, likelihood, flag);
    likeType.append(likelihood[0]);
    likeValue.append(likelihood[1]);
    */

    if ( _lookMCTruth.value()==true ) 
      { 
	int ind[4],nMatch;
	float cons[4];
	IndicesOfTruthMatched(candid,ind,cons,nMatch);
	TruInd.append(ind[0]);
      }
  }//while ( 0 != ( candid = iter_ifr()) )

  ntuple->column("ifr_Nb",BasicIfrList->length(),0 ,"Ifr",HTRange<int>(0,100));
  //ntuple->column("ifr_eLab",eLab,"ifr_Nb",0 , "Ifr");
  ntuple->column("ifr_thLab",thLab,"ifr_Nb",0 , "Ifr");
  ntuple->column("ifr_phLab",phiLab,"ifr_Nb",0 , "Ifr");
  //ntuple->column("ifr_uid",uid,"ifr_Nb",0 , "Ifr");
  ntuple->column("ifr_firstHit",firstHit,"ifr_Nb",0 , "Ifr");
  ntuple->column("ifr_lastHit",lastHit,"ifr_Nb",0 , "Ifr");
  ntuple->column("ifr_multip",multip,"ifr_Nb",0 , "Ifr");
  //ntuple->column("ifr_likeValue",likeValue,"ifr_Nb",0 , "Ifr");
  //ntuple->column("ifr_likeType",likeType,"ifr_Nb",0 , "Ifr");

  if ( _lookMCTruth.value()==true ) 
    {
      ntuple->column("ifr_TruInd",TruInd, "ifr_Nb",0 ,"Ifr");
    }

}


void
XSLReader::FillNtupleWithBumpListBlock( AbsEvent* anEvent, HepTuple* ntuple )
{
  HTValOrderedVector<float> theta,phi,e,lat,uid,s9s25,thUps,phUps,eUps,secMom,zernike42;
  HTValOrderedVector<int> nCrys,TruInd,TruIndOne,TruIndTwo,TruIndThree,isOK,Bit,nearTrk;
  HTValOrderedVector<float> TruCons,TruConsOne,TruConsTwo,TruConsThree;

  //Most infos are copyed from BtaGoodNeutralSelector
  BtaCandidate* candid;
  HepAListIterator<BtaCandidate> iter_bump(*BasicBumpList);  

  BtaBooster* CMS = new BtaBooster(_UpsLab);
  
  while ( 0 != ( candid = iter_bump()) ) {
    
    const BtaCalQual* CalQual = candid->getMicroAdapter()->getCalQual();

    //uid.append( (int)candid->uid() );
    theta.append((float) candid->p3().theta() );
    phi.append((float) candid->p3().phi() );
    e.append((float) candid->energy() );
    
    lat.append((float) CalQual->lateralMoment() );	
    s9s25.append((float) CalQual->s9s25() );
    nCrys.append( CalQual->nCrystals() );
    secMom.append((float) CalQual->secondMomentTP() );    
    zernike42.append( CalQual->absZernike42() );
    //isOK.append((int)(CalQual->isOk()) ? 1:0 );
    
    /*
    //attempt to migrate the Quals out of mu code... Abandonned for the moment!
    if(candid->recoCalo()!=0){
      cout<<"bonjour"<<endl;
      if(candid->recoCalo()->dynamic_cast_EmcCluster()!=0) {
	float latB = candid->recoCalo()->dynamic_cast_EmcCluster()->Xmoments().lat();
	cout<<"lat: "<<CalQual->lateralMoment()<<endl;
	cout<<"latB: "<<latB<<endl;
	float zer42B = candid->recoCalo()->dynamic_cast_EmcCluster()->Xmoments().absZernikeMoment(4,2);
	cout<<"CalQual->absZernike42(): "<<CalQual->absZernike42()<<endl;
	cout<<"zer42B: "<<zer42B<<endl;
	float s9s25B = candid->recoCalo()->dynamic_cast_EmcCluster()->esums().e9e25();
	cout<<"s9s25: "<<CalQual->s9s25() <<endl;
	cout<<"s9s25B: "<<s9s25B<<endl;
	float secMomB = candid->recoCalo()->dynamic_cast_EmcCluster()->moments().secondMomentTP();
	cout<<"CalQual->secondMomentTP(): "<<CalQual->secondMomentTP()<<endl;
	cout<<"secMomB: "<<secMomB<<endl;
	int nBump = candid->recoCalo()->dynamic_cast_EmcCluster()->nBumps();
	int nCryB = candid->recoCalo()->dynamic_cast_EmcCluster()->numberOfDigis();
	cout<<" CalQual->nCrystals(): "<< CalQual->nCrystals()<<endl;
	cout<<"nCryB: "<<nCryB<<endl;
      }
    }
    */

    BtaCandidate BumpUps = CMS->boostTo(*candid);
    eUps.append((float)BumpUps.energy());
    thUps.append((float)BumpUps.p3().theta());
    phUps.append((float)BumpUps.p3().phi());

    /*
    double likelihood[2];
    int flag = 1; //flag=1 : EMC cand, flag=2 : IFR cand 
    KLselector::Selector(anEvent, _runNb, candid, likelihood, flag);
    
    bool noIfrMatch = (likelihood[0]==1.3);
    //likeType.append(likelihood[0]);  //this info moved to emc_Bit (see line 400)
    likeValue.append(likelihood[1]);
    */
    
    //Get the nearest trk to this photon... We now use the "BasedOnAlpha" method
    int NearestTrkInd=NearestTrkIndBasedOnAlpha(candid);    
    //double distToClosestTrk=999;
    //BtaCandidate* nearestTrk = GetNearestTrk(anEvent,candid,distToClosestTrk);
    //if(nearestTrk!=0){ NearestTrkInd=CandIndexInList(nearestTrk,BasicTrkList); }//if(nearestTrk!=0)
    nearTrk.append(NearestTrkInd);        

    //Looking if this photon part of GPVE and storing the answer in emc_Bit
    static const IfdStrKey GPVEKey("GoodPhotonsVisibleE");
    HepAList< BtaCandidate >* GPVE = Ifd<HepAList<BtaCandidate> >::get(anEvent,GPVEKey);
    if(GPVE == 0 ) ErrMsg(fatal)<<"couldn't find GoodPhotonsVisibleE"<<endl;
    bool is_GPVE=CandInList( candid, GPVE );

    //The final EMC bit!
    int emcBit=0;
    //if( isEl==1 )          emcBit = emcBit | 1; 
    //if( noIfrMatch )       emcBit = emcBit | 2; 
    if( is_GPVE )          emcBit = emcBit | 4; 
    Bit.append(emcBit);

    
    if ( _lookMCTruth.value()==true ) 
      { 
	int ind[4],nMatch;
	float cons[4];
	bool QualityMatch = IndicesOfTruthMatched(candid,ind,cons,nMatch);
	if(!QualityMatch){ _nBadRecoed+=1;  _BadRecoed+=candid->p4(); }
	
	TruIndOne.append(ind[0]);
	TruIndTwo.append(ind[1]);
	TruIndThree.append(ind[2]);
	TruInd.append(ind[3]);
	TruConsOne.append(cons[0]);
	TruConsTwo.append(cons[1]);
	TruConsThree.append(cons[2]);
	TruCons.append(cons[3]);
      }
    
    
  }// while ( 0 != ( candid = iter_bump()) ) {
  
  
  //now filling the ntuple...
  ntuple->column("emc_Nb",BasicBumpList->length(),0 ,"Photon",HTRange<int>(0,100));
  ntuple->column("emc_eLab",e,"emc_Nb",0 , "Photon");  //note: this is identical to CalQual->ecalEnergy()
  ntuple->column("emc_thLab",theta,"emc_Nb",0 , "Photon");
  ntuple->column("emc_phLab",phi,"emc_Nb",0 , "Photon");
  ntuple->column("emc_eUps",eUps,"emc_Nb",0 , "Photon");
  ntuple->column("emc_thUps",thUps,"emc_Nb",0 , "Photon");
  ntuple->column("emc_phUps",phUps,"emc_Nb",0 , "Photon");
  ntuple->column("emc_lat",lat,"emc_Nb",0 , "Photon");
  //ntuple->column("emc_uid",uid,"emc_Nb",0 , "Photon");
  ntuple->column("emc_s9s25",s9s25,"emc_Nb",0 , "Photon");
  ntuple->column("emc_nCrys",nCrys,"emc_Nb",0 , "Photon");
  //ntuple->column("emc_isOK",isOK,"emc_Nb",0 , "Photon");
  ntuple->column("emc_nearTrkInd",nearTrk,"emc_Nb",0 , "Photon");
  ntuple->column("emc_Bit",Bit,"emc_Nb",0 , "Photon");
  ntuple->column("emc_secMom",secMom,"emc_Nb",0 , "Photon");
  ntuple->column("emc_zern42",zernike42,"emc_Nb",0 , "Photon");
  //ntuple->column("emc_likeValue",likeValue,"emc_Nb",0 , "Photon");
  //ntuple->column("emc_likeType",likeType,"emc_Nb",0 , "Photon");  ->info moved to emc_Bit

  if ( _lookMCTruth.value()==true ) 
    {
      ntuple->column("emc_TruInd",TruInd, "emc_Nb",0 ,"Photon");
      ntuple->column("emc_TruIndOne",TruIndOne, "emc_Nb",0 ,"Photon");
      ntuple->column("emc_TruIndTwo",TruIndTwo, "emc_Nb",0 ,"Photon");
      ntuple->column("emc_TruIndThree",TruIndThree, "emc_Nb",0 ,"Photon");
      ntuple->column("emc_TruCons",TruCons, "emc_Nb",0 ,"Photon");
      ntuple->column("emc_TruConsOne",TruConsOne, "emc_Nb",0 ,"Photon");
      ntuple->column("emc_TruConsTwo",TruConsTwo, "emc_Nb",0 ,"Photon");
      ntuple->column("emc_TruConsThree",TruConsThree, "emc_Nb",0 ,"Photon");
    }

  delete CMS;

  //voila! :-)  
}

int
XSLReader::NearestTrkIndBasedOnAlpha(BtaCandidate* neut){
  
  int indice=-66,ind=0;
  float cosDeltaAlphaMax=-2;
  BtaCandidate* trk;
  HepAListIterator<BtaCandidate> iterTrk(*BasicTrkList);
  while( 0!=(trk = iterTrk()) )
    {
      const BtaPidQual* trkPidQual=trk->getMicroAdapter()->getPidQual();
      float thEmc=trkPidQual->thetaAtEMC();
      float phiEmc=trkPidQual->phiAtEMC();

      float phiLabNeut=neut->p4().phi();
      float thetaLabNeut=neut->p4().theta();
		  
      float cosDeltaAlpha = cos(thEmc)*cos(thetaLabNeut) + sin(thEmc)*sin(thetaLabNeut)*cos(phiEmc-phiLabNeut);
      
      if(cosDeltaAlpha > cosDeltaAlphaMax){ indice=ind; }
      ind+=1;
    }
  
  return indice;
}


BtaCandidate*
XSLReader::GetNearestTrk( AbsEvent* anEvent, BtaCandidate* c, double &distToClosestTrk )
{
  //This method returns the nearest trk unmatched to any photon
  BtaCandidate* theNearestTrk=0;

  const BtaCalQual* calQual = c->getMicroAdapter()->getCalQual();
  HepPoint xyzCal(calQual->centroid()); //cluster position
  double remc(sqrt(xyzCal.x()*xyzCal.x()+xyzCal.y()*xyzCal.y()));

  distToClosestTrk=999;
  BtaCandidate* Cand;
  HepAListIterator<BtaCandidate> iterTrk(*BasicTrkList);
  while( 0!=(Cand = iterTrk()) )
    {
      if ( Cand->p4().perp()<0.15 ) continue;
      
      const BtaPidQual* trkPidQual=Cand->getMicroAdapter()->getPidQual();
      if(!trkPidQual) continue;
      double theta(trkPidQual->thetaAtEMC());
      if(theta<0.1) continue;
      double phi(trkPidQual->phiAtEMC());
      const HepPoint pointOnTraj(remc*sin(theta)*cos(phi),
				 remc*sin(theta)*sin(phi),
				 remc*cos(theta));
      const double missDistance((pointOnTraj-xyzCal).mag());
      
      if( missDistance<distToClosestTrk ) 
	{ 
	  theNearestTrk=Cand; 
	  distToClosestTrk=missDistance;
	}
    }

  return theNearestTrk;
}



void
XSLReader::GetPidBit( AbsEvent* anEvent, BtaCandidate* c, int *bita, int *bitb )
{
  //Note: PidBit does NOT only contains information about the PID lists, but also about GTVE and other charged tracks lists

  //Electrons
  static const IfdStrKey eMicroVeryLooseKey("eMicroVeryLoose");
  HepAList<BtaCandidate>* eMicroVeryLoose = Ifd<HepAList<BtaCandidate> >::get(anEvent,eMicroVeryLooseKey);
  if(eMicroVeryLoose == 0 ) ErrMsg(fatal)<<"couldn't find eMicroVeryLoose"<<endl;
  bool is_eMicroVeryLoose=CandInList( c, eMicroVeryLoose );
 
  static const IfdStrKey eMicroLooseKey("eMicroLoose");
  HepAList<BtaCandidate>* eMicroLoose = Ifd<HepAList<BtaCandidate> >::get(anEvent,eMicroLooseKey);
  if(eMicroLoose == 0 ) ErrMsg(fatal)<<"couldn't find eMicroLoose"<<endl;
  bool is_eMicroLoose=CandInList( c, eMicroLoose );
 
  static const IfdStrKey eMicroTightKey("eMicroTight");
  HepAList<BtaCandidate>* eMicroTight = Ifd<HepAList<BtaCandidate> >::get(anEvent,eMicroTightKey);
  if(eMicroTight == 0 ) ErrMsg(fatal)<<"couldn't find eMicroTight"<<endl;
  bool is_eMicroTight=CandInList( c, eMicroTight );
 
  static const IfdStrKey eMicroVeryTightKey("eMicroVeryTight");
  HepAList<BtaCandidate>* eMicroVeryTight = Ifd<HepAList<BtaCandidate> >::get(anEvent,eMicroVeryTightKey);
  if(eMicroVeryTight == 0 ) ErrMsg(fatal)<<"couldn't find eMicroVeryTight"<<endl;
  bool is_eMicroVeryTight=CandInList( c, eMicroVeryTight );
 
  static const IfdStrKey eTightLHKey("PidLHElectrons");
  HepAList<BtaCandidate>* eTightLH = Ifd<HepAList<BtaCandidate> >::get(anEvent,eTightLHKey);
  if(eTightLH == 0 ) ErrMsg(fatal)<<"couldn't find PidLHElectrons"<<endl;
  bool is_PidLHElectrons=CandInList( c, eTightLH );
  
 
  //muons
  static const IfdStrKey muMicroTightKey("muMicroTight");
  HepAList<BtaCandidate>* muMicroTight = Ifd<HepAList<BtaCandidate> >::get(anEvent,muMicroTightKey);
  if(muMicroTight == 0 ) ErrMsg(fatal)<<"couldn't find muMicroTight"<<endl;
  bool is_muMicroTight=CandInList( c, muMicroTight );
  
  static const IfdStrKey muMicroVeryTightKey("muMicroVeryTight");
  HepAList<BtaCandidate>* muMicroVeryTight = Ifd<HepAList<BtaCandidate> >::get(anEvent,muMicroVeryTightKey);
  if(muMicroVeryTight == 0 ) ErrMsg(fatal)<<"couldn't find muMicroVeryTight"<<endl;
  bool is_muMicroVeryTight=CandInList( c, muMicroVeryTight );
  

  static const IfdStrKey muNNVeryLooseKey("muNNVeryLoose");
  HepAList<BtaCandidate>* muNNVeryLoose = Ifd<HepAList<BtaCandidate> >::get(anEvent,muNNVeryLooseKey);
  if(muNNVeryLoose == 0 ) ErrMsg(fatal)<<"couldn't find muNNVeryLoose"<<endl;
  bool is_muNNVeryLoose=CandInList( c, muNNVeryLoose );
  
  bool is_muNNLoose=false,is_muNNTight=false,is_muNNVeryTight=false;
  static const IfdStrKey muNNLooseKey("muNNLoose");
  HepAList<BtaCandidate>* muNNLoose = Ifd<HepAList<BtaCandidate> >::get(anEvent,muNNLooseKey);
  if(muNNLoose == 0 ) ErrMsg(fatal)<<"couldn't find muNNLoose"<<endl;
  is_muNNLoose=CandInList( c, muNNLoose );
  
  static const IfdStrKey muNNTightKey("muNNTight");
  HepAList<BtaCandidate>* muNNTight = Ifd<HepAList<BtaCandidate> >::get(anEvent,muNNTightKey);
  if(muNNTight == 0 ) ErrMsg(fatal)<<"couldn't find muNNTight"<<endl;
  is_muNNTight=CandInList( c, muNNTight );
  
  static const IfdStrKey muNNVeryTightKey("muNNVeryTight");
  HepAList<BtaCandidate>* muNNVeryTight = Ifd<HepAList<BtaCandidate> >::get(anEvent,muNNVeryTightKey);
  if(muNNVeryTight == 0 ) ErrMsg(fatal)<<"couldn't find muNNVeryTight"<<endl;
  is_muNNVeryTight=CandInList( c, muNNVeryTight );
  
  static const IfdStrKey muNNVeryLooseFakeRateKey("muNNVeryLooseFakeRate");
  HepAList<BtaCandidate>* muNNVeryLooseFakeRate = Ifd<HepAList<BtaCandidate> >::get(anEvent,muNNVeryLooseFakeRateKey);
  if(muNNVeryLooseFakeRate == 0 ) ErrMsg(fatal)<<"couldn't find muNNVeryLooseFakeRate"<<endl;
  bool is_muNNVeryLooseFakeRate=CandInList( c, muNNVeryLooseFakeRate );
  
  bool is_muNNLooseFakeRate=false,is_muNNTightFakeRate=false,is_muNNVeryTightFakeRate=false;
  static const IfdStrKey muNNLooseFakeRateKey("muNNLooseFakeRate");
  HepAList<BtaCandidate>* muNNLooseFakeRate = Ifd<HepAList<BtaCandidate> >::get(anEvent,muNNLooseFakeRateKey);
  if(muNNLooseFakeRate == 0 ) ErrMsg(fatal)<<"couldn't find muNNLooseFakeRate"<<endl;
  is_muNNLooseFakeRate=CandInList( c, muNNLooseFakeRate );
  
  static const IfdStrKey muNNTightFakeRateKey("muNNTightFakeRate");
  HepAList<BtaCandidate>* muNNTightFakeRate = Ifd<HepAList<BtaCandidate> >::get(anEvent,muNNTightFakeRateKey);
  if(muNNTightFakeRate == 0 ) ErrMsg(fatal)<<"couldn't find muNNTightFakeRate"<<endl;
  is_muNNTightFakeRate=CandInList( c, muNNTightFakeRate );
  
  static const IfdStrKey muNNVeryTightFakeRateKey("muNNVeryTightFakeRate");
  HepAList<BtaCandidate>* muNNVeryTightFakeRate = Ifd<HepAList<BtaCandidate> >::get(anEvent,muNNVeryTightFakeRateKey);
  if(muNNVeryTightFakeRate == 0 ) ErrMsg(fatal)<<"couldn't find muNNVeryTightFakeRate"<<endl;
  is_muNNVeryTightFakeRate=CandInList( c, muNNVeryTightFakeRate );
  
  //Kaons
  //HepAList<BtaCandidate>* KLHVeryLoose = Ifd<HepAList<BtaCandidate> >::get(anEvent,IfdStrKey(HepString("KLHVeryLoose")));
  //if(KLHVeryLoose == 0 ) ErrMsg(fatal)<<"couldn't find KLHVeryLoose"<<endl;
  //bool is_KLHVeryLoose=CandInList( c, KLHVeryLoose );

  static const IfdStrKey KLHLooseKey("KLHLoose");
  HepAList<BtaCandidate>* KLHLoose = Ifd<HepAList<BtaCandidate> >::get(anEvent,KLHLooseKey);
  if(KLHLoose == 0 ) ErrMsg(fatal)<<"couldn't find KLHLoose"<<endl;
  bool is_KLHLoose=CandInList( c, KLHLoose );
  
  bool is_KLHTight=false,is_KLHVeryTight=false;
  static const IfdStrKey KLHTightKey("KLHTight");
  HepAList<BtaCandidate>* KLHTight = Ifd<HepAList<BtaCandidate> >::get(anEvent,KLHTightKey);
  if(KLHTight == 0 ) ErrMsg(fatal)<<"couldn't find KLHTight"<<endl;
  is_KLHTight=CandInList( c, KLHTight );

  static const IfdStrKey KLHVeryTightKey("KLHVeryTight");
  HepAList<BtaCandidate>* KLHVeryTight = Ifd<HepAList<BtaCandidate> >::get(anEvent,KLHVeryTightKey);
  if(KLHVeryTight == 0 ) ErrMsg(fatal)<<"couldn't find KLHVeryTight"<<endl;
  is_KLHVeryTight=CandInList( c, KLHVeryTight );
  
  //protons
  //HepAList<BtaCandidate>* pLHVeryLoose = Ifd<HepAList<BtaCandidate> >::get(anEvent,IfdStrKey(HepString("pLHVeryLoose")));
  //if(pLHVeryLoose == 0 ) ErrMsg(fatal)<<"couldn't find pLHVeryLoose"<<endl;
  //bool is_pLHVeryLoose=CandInList( c, pLHVeryLoose );

  static const IfdStrKey pLHLooseKey("pLHLoose");
  HepAList<BtaCandidate>* pLHLoose = Ifd<HepAList<BtaCandidate> >::get(anEvent,pLHLooseKey);
  if(pLHLoose == 0 ) ErrMsg(fatal)<<"couldn't find pLHLoose"<<endl;
  bool is_pLHLoose=CandInList( c, pLHLoose );

  bool is_pLHTight=false,is_pLHVeryTight=false;
  static const IfdStrKey pLHTightKey("pLHTight");
  HepAList<BtaCandidate>* pLHTight = Ifd<HepAList<BtaCandidate> >::get(anEvent,pLHTightKey);
  if(pLHTight == 0 ) ErrMsg(fatal)<<"couldn't find pLHTight"<<endl;
  is_pLHTight=CandInList( c, pLHTight );
  
  static const IfdStrKey pLHVeryTightKey("pLHVeryTight");
  HepAList<BtaCandidate>* pLHVeryTight = Ifd<HepAList<BtaCandidate> >::get(anEvent,pLHVeryTightKey);
  if(pLHVeryTight == 0 ) ErrMsg(fatal)<<"couldn't find pLHVeryTight"<<endl;
  is_pLHVeryTight=CandInList( c, pLHVeryTight );
  
  //pions
  static const IfdStrKey piLHLooseKey("piLHLoose");
  HepAList<BtaCandidate>* piLHLoose = Ifd<HepAList<BtaCandidate> >::get(anEvent,piLHLooseKey);
  if(piLHLoose == 0 ) ErrMsg(fatal)<<"couldn't find piLHLoose"<<endl;
  bool is_piLHLoose=CandInList( c, piLHLoose );

  bool is_piLHTight=false,is_piLHVeryTight=false;
  static const IfdStrKey piLHTightKey("piLHTight");
  HepAList<BtaCandidate>* piLHTight = Ifd<HepAList<BtaCandidate> >::get(anEvent,piLHTightKey);
  if(piLHTight == 0 ) ErrMsg(fatal)<<"couldn't find piLHTight"<<endl;
  is_piLHTight=CandInList( c, piLHTight );
  
  static const IfdStrKey piLHVeryTightKey("piLHVeryTight");
  HepAList<BtaCandidate>* piLHVeryTight = Ifd<HepAList<BtaCandidate> >::get(anEvent,piLHVeryTightKey);
  if(piLHVeryTight == 0 ) ErrMsg(fatal)<<"couldn't find piLHVeryTight"<<endl;
  is_piLHVeryTight=CandInList( c, piLHVeryTight );
  
  //Tracking lists
  //These lists ain't nothing to do with PID... but are still useful!
  static const IfdStrKey GTVEKey("GoodTracksVisibleE");
  HepAList<BtaCandidate>* GTVE = Ifd<HepAList<BtaCandidate> >::get(anEvent,GTVEKey);
  if(GTVE == 0 ) ErrMsg(fatal)<<"couldn't find GoodTracksVisibleE"<<endl;
  bool is_GTVE=CandInList( c, GTVE );
  

  static const IfdStrKey GTVLKey("GoodTracksVeryLoose");
  HepAList< BtaCandidate >* GTVL = Ifd<HepAList<BtaCandidate> >::get(anEvent,GTVLKey);
  if(GTVL == 0 ) ErrMsg(fatal)<<"couldn't find GoodTracksVeryLoose"<<endl;
  bool is_GTVL=CandInList( c, GTVL );

  bool is_GTL=false;
  if(is_GTVL) {
    static const IfdStrKey GTLKey("GoodTracksLoose");
    HepAList< BtaCandidate >* GTL = Ifd<HepAList<BtaCandidate> >::get(anEvent,GTLKey);
    if(GTL == 0 ) ErrMsg(fatal)<<"couldn't find GoodTracksLoose"<<endl;
    is_GTL=CandInList( c, GTL );
  }//if(is_GTVL)


  static const IfdStrKey GamConvKey("gammaConversionDefault");
  HepAList<BtaCandidate>* GamConv = Ifd<HepAList<BtaCandidate> >::get(anEvent,GamConvKey);
  if(GamConv == 0 ) ErrMsg(fatal)<<"couldn't find GoodTracksVisibleE"<<endl;
  bool is_Gam=CandIsDaughterInList( c, GamConv );


  *bita=0;
  //if ( is_XXX )       *bita = *bita | 1; 
  if ( is_PidLHElectrons )    *bita = *bita | 2;         //corrected with weight
  if ( is_muMicroTight )      *bita = *bita | 4;         //corrected with weight
  if ( is_muMicroVeryTight )  *bita = *bita | 8;         //corrected with weight

  //if ( is_pLHVeryLoose )      *bita = *bita | 16;        //corrected with tweaking
  if ( is_pLHLoose )          *bita = *bita | 32;        //corrected with tweaking
  if ( is_pLHTight )          *bita = *bita | 64;        //corrected with tweaking
  if ( is_pLHVeryTight )      *bita = *bita | 128;       //corrected with tweaking

  if ( is_piLHLoose )         *bita = *bita | 256;       //corrected with weight
  if ( is_piLHTight )         *bita = *bita | 512;       //corrected with tweaking
  if ( is_piLHVeryTight )     *bita = *bita | 1024;      //corrected with tweaking

  if ( is_GTVE )              *bita = *bita | 2048;      //not corrected for data/MC
  if ( is_GTVL )              *bita = *bita | 4096;      //not corrected for data/MC
  //The GTL bit value is used elsewhere. Do NOT modify!
  if ( is_GTL )               *bita = *bita | 8192;      //not corrected for data/MC
  if ( is_Gam )               *bita = *bita | 16384;     //not corrected for data/MC


  *bitb=0;
  if ( is_eMicroVeryLoose )       *bitb = *bitb | 1;      //corrected with tweaking
  if ( is_eMicroLoose )           *bitb = *bitb | 2;      //corrected with tweaking
  if ( is_eMicroTight )           *bitb = *bitb | 4;      //corrected with tweaking
  if ( is_eMicroVeryTight )       *bitb = *bitb | 8;      //corrected with tweaking
  if( is_muNNVeryLooseFakeRate )  *bitb = *bitb | 16;     //corrected with tweaking
  if( is_muNNVeryLoose )          *bitb = *bitb | 32;     //corrected with tweaking
  if( is_muNNLooseFakeRate )      *bitb = *bitb | 64;     //corrected with tweaking
  if( is_muNNLoose )              *bitb = *bitb | 128;    //corrected with tweaking
  if( is_muNNTightFakeRate )      *bitb = *bitb | 256;    //corrected with tweaking
  if( is_muNNTight )              *bitb = *bitb | 512;    //corrected with tweaking
  if( is_muNNVeryTightFakeRate )  *bitb = *bitb | 1024;   //corrected with tweaking
  if( is_muNNVeryTight )          *bitb = *bitb | 2048;   //corrected with tweaking
  if( is_KLHLoose )               *bitb = *bitb | 4096;   //corrected with tweaking
  if( is_KLHTight )               *bitb = *bitb | 8192;   //corrected with tweaking
  if( is_KLHVeryTight )           *bitb = *bitb | 16384;  //corrected with tweaking

  return;
}

float
XSLReader::EncodeValAndErrInFloat(float val, float err)
{
  //note: for w = 1.2346 +/- 0.5615, this returns code=1235.562  
  //To read back:
  /*
  int iRead = (int)code;
  float rVal = (float)iRead/1000.0;
  cout<<"rVal: "<<rVal<<endl;

  float W= code*1000;  //
  iRead = (int)W;     //important to do this in two steps!
  iRead = iRead%1000;
  float rErr = (float)iRead/1000.0;
  cout<<"rErr: "<<rErr<<endl;
  */

  float code=0.0;
  int roundup=0;
  
  float fVal;
  if(val>9.999) fVal=9999;
  else if(val<0) fVal=1000;
  else {
    long int tmp = (long int)(val*10000);
    roundup = ( (tmp)%10 >4 ) ? 1:0;  
    long int iVal = (long int)(val*1000)+roundup;
    fVal=(float)iVal;
  }
  
  float fErr;
  if(err>0.999 || err<0) fErr=0.999;  //there is occasional Err=-1 cases...
  else {
    long int tmp=(long int)(err*10000);
    roundup = ( (tmp)%10 >4 ) ? 1:0;  
    long int iErr = (long int)(err*1000)+roundup;  
    fErr = (float)iErr/1000.0;
  }

  code = fVal + fErr;
  return code;
}

void
XSLReader::FillNtupleWithTrkListBlock( AbsEvent* anEvent, HepTuple* ntuple )
{
  //Trk Eff
  EffTable* table=NULL;
  //Thomas: "The multiplicity is the number of GTVL tracks
  //in the event ( GTVL is our reference even as the tables are derived for GTL tracks )"
  int Multiplicity=0;
  //Thomas: "One minor issue on the angles - I had to change from rad to degree to
  //follow the tables of the PID group so you would need a factor angleScale = 180.0/3.14159"
  double angleScale = 180.0/3.14159;

  //PID weights
  std::map< int,PidWeight> *eTightLH_Map=NULL; 
  std::map< int,PidWeight> *muMicroTight_Map=NULL; 
  std::map< int,PidWeight> *muMicroVeryTight_Map=NULL; 
  std::map< int,PidWeight> *piLooseLH_Map=NULL; 

  if(_sigModeString!="data") {
    
    ////////////////////////
    //  TrkEff Tables
    //See: http://www.slac.stanford.edu/BFROOT/www/Organization/CollabMtgs/2004/detMay04/Wed4i/tracking.pdf
    static const IfdStrKey tK("/physicstools/trk/svt_tables/all");  //also possible to use all,pos,neg (all used to be named "average")
    table = Ifd<EffTable>::get(gblPEnv, tK);

    //Thomas: "The multiplicity is the number of GTVL tracks
    //in the event ( GTVL is our reference even as the tables are derived for GTL tracks )"
    HepAList< BtaCandidate >* GTVL;
    getTmpAList( anEvent, GTVL, IfdStrKey(HepString("GoodTracksVeryLoose")));
    Multiplicity = GTVL->length();
    
    
    //////////////////////////
    //  PID weights
    //See: http://www.slac.stanford.edu/BFROOT/www/Physics/Tools/Pid/PidOnMc/pidonmc.html#weight
    eTightLH_Map = Ifd< map< int,PidWeight> >::get( anEvent, "PidLHElectrons"); 
    muMicroTight_Map = Ifd< map< int,PidWeight> >::get( anEvent, "muNNTight");  //switched from cut-based to NN selector!  
    muMicroVeryTight_Map = Ifd< map< int,PidWeight> >::get( anEvent, "muNNVeryTight"); 
    piLooseLH_Map = Ifd< map< int,PidWeight> >::get( anEvent, "piLHLoose"); 

  }//if(_sigModeString!="data")


  HTValOrderedVector<float> e,p,th,Phi,pt,chi2,docaXY,docaZ,uid,pUps,eUps,thUps,phUps,eCalo;//trkEff,eEff,muTEff,muVTEff,piEff;
  HTValOrderedVector<float> delThGhost,delPhGhost,delPtGhost,delThLoopSame,delPhLoopSame,delThLoopOpp,delPhLoopOpp,delPtLoop;
  HTValOrderedVector<int> TruInd,TruIndOne,TruIndTwo,TruIndThree,DataLund,Ch,PidBitA,PidBitB,nSVT,nDCH;
  HTValOrderedVector<float> thAtEmc,phAtEmc,trkW,trkWE,eW,eWE,muTW,muTWE,muVTW,muVTWE,piW,piWE;
  HTValOrderedVector<float> TruCons,TruConsOne,TruConsTwo,TruConsThree;

  BtaBooster* CMS = new BtaBooster(_UpsLab);

  //Most infos are copyed from BtaGoodTrkSelector
  BtaCandidate* trk;
  HepAListIterator<BtaCandidate> iter_trk(*BasicTrkList);  

  while ( 0 != ( trk = iter_trk()) ) {

    //DataLund.append((int)trk->pdtEntry()->lundId());

    //Note: PieBit does NOT only contains information about the PID lists, but also about GTVE
    int bita=0,bitb=0;
    GetPidBit(anEvent,trk,&bita,&bitb);
    PidBitA.append( bita );
    PidBitB.append( bitb );

    e.append((float)trk->energy());
    p.append((float)trk->p());
    double ch = trk->charge();
    Ch.append((int)ch);

    //  information from fit 
    const TrkAbsFit * aTrackFit = trk->trkAbsFit();
    double pT = aTrackFit->pt(); 
    //pt.append((float)pT);
    double phi = aTrackFit->momentum(0).phi();
    Phi.append((float)phi);    
    double theta = aTrackFit->momentum(0).theta();
    th.append((float)theta);

    const BtaPidQual* trkPidQual=trk->getMicroAdapter()->getPidQual();
    float thEmc=-666,phEmc=-666;
    if(trkPidQual!=NULL) {
      thEmc=(float)trkPidQual->thetaAtEMC();
      phEmc=(float)trkPidQual->phiAtEMC();
    }
    thAtEmc.append( thEmc );
    phAtEmc.append( phEmc );    


    if(_sigModeString!="data")
      {
	double eff(1.0),err(0);
	std::vector<double> vect(4,0); 
	vect[0]= pT;  //pT LAB
	vect[1]= theta*angleScale; //theta LAB in degree
	vect[2]= phi*angleScale; //phi LAB in degree
	vect[3]= Multiplicity; //Nb of GTVL in the event 
	if(table!=0) {
	  if((bita&8192)==8192) table->get_efficiency_and_error( eff, err, vect);  //weights only for GTL trks
	} else cout<<"table is a NULL pointer!!"<<endl;
	//trkEff.append( EncodeValAndErrInFloat(eff,err) );	
	trkW.append((float)eff);
	trkWE.append((float)err);
	
	float eWW=1.0,eWWE=0.0;
	if(eTightLH_Map!=NULL){
	  PidWeight w = (*eTightLH_Map)[trk->uid()]; 
	  if(w._status==PidWeight::ok) { eWW=w._value; eWWE=w._error; }
	}
	//eEff.append(EncodeValAndErrInFloat(eWW,eWWE));
	eW.append(eWW);
	eWE.append(eWWE);

	float muTWW=1.0,muTWWE=0.0;
	if(muMicroTight_Map!=NULL){
	  PidWeight w = (*muMicroTight_Map)[trk->uid()]; 
	  if(w._status==PidWeight::ok) { muTWW=w._value; muTWWE=w._error; }
	} 
	//muTEff.append(EncodeValAndErrInFloat(muTWW,muTWWE));
	muTW.append(muTWW);
	muTWE.append(muTWWE);

	float muVTWW=1.0,muVTWWE=0.0;
	if(muMicroVeryTight_Map!=NULL){
	  PidWeight w = (*muMicroVeryTight_Map)[trk->uid()]; 
	  if(w._status==PidWeight::ok) { muVTWW=w._value; muVTWWE=w._error; }
	}
	//muVTEff.append(EncodeValAndErrInFloat(muVTWW,muVTWWE));
	muVTW.append(muVTWW);
	muVTWE.append(muVTWWE);


	float piWW=1.0,piWWE=0.0;
	if(piLooseLH_Map!=NULL){
	  PidWeight w = (*piLooseLH_Map)[trk->uid()]; 
	  if(w._status==PidWeight::ok) { piWW=w._value; piWWE=w._error; }
	}
	//piEff.append(EncodeValAndErrInFloat(piWW,piWWE));
	piW.append(piWW);
	piWE.append(piWWE);

      }

    //info from micro DB
    const BtaTrkQual* curTrkQual = trk->getMicroAdapter()->getTrkQual();
    int nHitsSvtTest = curTrkQual->nSvtHits();
    int nHitsDchTest = curTrkQual->nDchHits();
    nSVT.append(nHitsSvtTest);
    nDCH.append(nHitsDchTest);
    chi2.append((float)curTrkQual->prob());

    //attempt to migrate the Quals out of mu code... Abandonned for the moment!
    //int nSVTB = trk->recoTrk()->hots()->nSvt();
    //int nDCHB = trk->recoTrk()->hots()->nDch();
    //float chi2B = trk->trkAbsFit()->chisq(); also nDof()

    const BtaCalQual* CalQual = trk->getMicroAdapter()->getCalQual();
    float eCal=-666.0;
    if(CalQual!=NULL) eCal = CalQual->ecalEnergy();
    eCalo.append(eCal);

    //  BeamSpot (for DOCA)
    const BbrLorentzVectorErr & ip = gblEnv->getBta()->pepBeams()->interactionPoint();
    HepPoint beamSpotPoint(ip.x(), ip.y(), ip.z());
    
    //DOCA
    HepPoint closest= trk->recoObject()->position(beamSpotPoint,BtaAbsRecoObject::XY);
    double docaTest = sqrt(pow((closest.x()-beamSpotPoint.x()),2)+pow((closest.y()-beamSpotPoint.y()),2));
    docaXY.append((float)docaTest);
    double zTest = closest.z()-beamSpotPoint.z(); 
    docaZ.append((float)zTest);

    int UID = (int)trk->uid();  //used later on in the function... (by the way)
    //uid.append(UID);

    BtaCandidate TrkUps = CMS->boostTo(*trk);
    eUps.append((float)TrkUps.energy());
    pUps.append((float)TrkUps.p());
    thUps.append((float)TrkUps.p3().theta());
    phUps.append((float)TrkUps.p3().phi());
    

    //////////
    //Ghost - Following Bob's algorithm, we store the deltaTheta/Phi/Pt of the "ghostest" trk 
    double smallDelThGh=1000,smallDelPhiGh=1000,smallDelPtGh=1000;
    if( nHitsDchTest>=_nDchMinGhost.value() && pT<_ptMaxGhost.value() )
      {
	if(phi<0) phi+=2.*3.1415926;
		
	BtaCandidate* Cand;
	HepAListIterator<BtaCandidate> trkIter(*BasicTrkList);
	while( 0!=(Cand = trkIter()) )
	  {
	    if(UID!=(int)Cand->uid() && ch==Cand->charge() )
	      {
		
		const BtaTrkQual* curTrkQual2 = Cand->getMicroAdapter()->getTrkQual();
		int nDch2 = curTrkQual2->nDchHits();
		
		//Even though a ghost trk is found, this trk won't be rejected if it has more DCH hits
		if(nDch2>nHitsDchTest)
		  {
		    double pT2 = Cand->p4().perp();
		    double theta2= Cand->p4().theta();
		    double phi2= Cand->p4().phi();
		    if(phi2<0) phi2+=2.*3.1415926;
		    
		    double deltaPhi = phi2-phi;
		    if(deltaPhi<(-2.)*3.1415926) deltaPhi+=2*3.1415926; 
		    else if(deltaPhi>2*3.1415926) deltaPhi-=2*3.1415926; 
		    double deltaTheta=theta2-theta;
		    double deltaPt=pT-pT2;
		    
		    int nDchCriteria = 45-nDch2;
		    
		    //The cut!
		    //This trk will be rejected if a ghost candidate with more DCH hits is found
		    if( nHitsDchTest<nDchCriteria && fabs(deltaPhi)<smallDelPhiGh
			&& fabs(deltaTheta)<smallDelThGh
			&& fabs(deltaPt)<smallDelPtGh )
		      {
			smallDelPhiGh = deltaPhi;
			smallDelThGh = deltaTheta;
			smallDelPtGh = deltaPt; 
		      }
		  }//if(nDch2>nHitsDchTest)
	      }//if(UID!=(int)Cand->uid())
	  }//while 	  
      }//Ghost
    delThGhost.append((float)smallDelThGh);
    delPhGhost.append((float)smallDelPhiGh);
    delPtGhost.append((float)smallDelPtGh);
    
    ////////////
    //Loopers
    double cosTh = cos(theta);
    double smallDelThLoopSame=1000,smallDelPhLoopSame=1000,smallDelThLoopOpp=1000,smallDelPhLoopOpp=1000,smallDelPtLoop=1000;
    if( nHitsSvtTest>=_nSvtMinLooper.value() && fabs(cosTh)<_cosThetaMaxLooper.value() && pT<_ptMaxLooper.value() )
      {
	if(phi<0) phi+=2.*3.1415926;
	
	
	BtaCandidate* Cand;
	HepAListIterator<BtaCandidate> trkIter(*BasicTrkList);
	while( 0!=(Cand = trkIter()) )
	  {
	    double cosTh2=Cand->p3().cosTheta();
	    double pT2=Cand->pt();
	    double deltaPT=pT-pT2;
	    if(UID!=(int)Cand->uid()
	       && fabs(cosTh2)<_cosThetaMaxLooper.value() && pT2<_ptMaxLooper.value() )
	      {
		double phi2= Cand->p4().phi();
		if(phi2<0) phi2+=2.*3.1415926;
		double deltaPhi = phi2-phi;
		if(deltaPhi<(-2.)*3.1415926) deltaPhi+=2*3.1415926; 
		if(deltaPhi>2*3.1415926) deltaPhi-=2*3.1415926; 
		if(deltaPhi>3.1415926/2) deltaPhi = 3.1415926 - deltaPhi;
		if(deltaPhi<(-3.1415926/2) ) deltaPhi = 3.1415926 + deltaPhi;		  
		
		int ch2= (int)Cand->charge();
		double theta2= Cand->p4().theta();
		double deltaTheta=0,deltaThetaMax=1000,deltaPhiMax=1000;
		if(ch*ch2>0) 
		  {
		    deltaTheta = theta2-theta;
		    deltaThetaMax = smallDelThLoopSame;
		    deltaPhiMax = smallDelPhLoopSame;
		  }
		else if(ch*ch2<0)
		  {
		    deltaTheta = 3.1415926-theta2-theta;
		    deltaThetaMax = smallDelThLoopOpp;
		    deltaPhiMax = smallDelPhLoopOpp;
		  }
		else ErrMsg(error) << "Error in BetaMicro/BtaGoodTrkSelector!! A trk with charge 0!!" << endmsg ;
		
		if( fabs(deltaTheta)<deltaThetaMax && fabs(deltaPhi)<deltaPhiMax  && fabs(deltaPT)<smallDelPtLoop ) 
		  { 
		    //This trk is part of a looper!
		    HepPoint closest2 = Cand->recoObject()->position(beamSpotPoint,BtaAbsRecoObject::XY);
		    double z0 = fabs(closest2.z()-beamSpotPoint.z()); 
		    
		    //Only the trk with smallest z0 is kept...
		    if( z0<fabs(zTest) ) 
		      {
			smallDelPtLoop=fabs(deltaPT);
			if(ch*ch2>0){ smallDelThLoopSame=deltaTheta;  smallDelPhLoopSame=deltaPhi;   }
			else{ smallDelThLoopOpp=deltaTheta;  smallDelPhLoopOpp=deltaPhi;   }
					  
		      }
		  }// if(fabs(deltaTheta)<0.16 && fabs(deltaPhi)<0.18)
	      }//if(UID!=(int)Cand->uid())
	  }//while
      }// if( nHitsSvtTest!=0 && fabs(cosTh)<0.2 && pT<0.25 )
    delThLoopSame.append((float)smallDelThLoopSame);
    delPhLoopSame.append((float)smallDelPhLoopSame);
    delThLoopOpp.append((float)smallDelThLoopOpp);
    delPhLoopOpp.append((float)smallDelPhLoopOpp);
    delPtLoop.append((float)smallDelPtLoop);


    if ( _lookMCTruth.value()==true ) 
      { 
	int ind[4],nMatch;
	float cons[4];
	bool QualityMatch = IndicesOfTruthMatched(trk,ind,cons,nMatch);
	if(!QualityMatch){ _nBadRecoed+=1;  _BadRecoed+=trk->p4(); }

	TruIndOne.append(ind[0]);
	TruIndTwo.append(ind[1]);
	TruIndThree.append(ind[2]);
	TruInd.append(ind[3]);
	TruConsOne.append(cons[0]);
	TruConsTwo.append(cons[1]);
	TruConsThree.append(cons[2]);
	TruCons.append(cons[3]);
      }
    
  }//end of track list
  
  
  //Now filling the ntuple...
  ntuple->column("trk_Nb",BasicTrkList->length(),0 ,"Track",HTRange<int>(0,100));
  ntuple->column("trk_Bit",PidBitA, "trk_Nb",0 ,"Track");
  ntuple->column("trk_BitB",PidBitB, "trk_Nb",0 ,"Track");
  ntuple->column("trk_eLab",e, "trk_Nb",0 ,"Track");
  ntuple->column("trk_pLab",p, "trk_Nb",0 ,"Track");
  ntuple->column("trk_thAtEmc",thAtEmc, "trk_Nb",0 ,"Track");
  ntuple->column("trk_phiAtEmc",phAtEmc, "trk_Nb",0 ,"Track");
  ntuple->column("trk_thLab",th, "trk_Nb",0 ,"Track");
  ntuple->column("trk_phiLab",Phi, "trk_Nb",0 ,"Track");
  //ntuple->column("trk_ptLab",pt, "trk_Nb",0 ,"Track");
  ntuple->column("trk_eUps",eUps, "trk_Nb",0 ,"Track");
  ntuple->column("trk_pUps",pUps, "trk_Nb",0 ,"Track");
  ntuple->column("trk_thUps",thUps, "trk_Nb",0 ,"Track");
  ntuple->column("trk_phiUps",phUps, "trk_Nb",0 ,"Track");
  ntuple->column("trk_ch",Ch, "trk_Nb",0 ,"Track");
  ntuple->column("trk_nSVT",nSVT, "trk_Nb",0 ,"Track");
  ntuple->column("trk_nDCH",nDCH, "trk_Nb",0 ,"Track");
  ntuple->column("trk_chi2",chi2, "trk_Nb",0 ,"Track");
  ntuple->column("trk_eCal",eCalo, "trk_Nb",0 ,"Track");
  ntuple->column("trk_docaXY",docaXY, "trk_Nb",0 ,"Track");
  ntuple->column("trk_docaZ",docaZ, "trk_Nb",0 ,"Track");
  //ntuple->column("trk_uid",uid, "trk_Nb",0 ,"Track");
  ntuple->column("trk_delThGhost",delThGhost, "trk_Nb",0 ,"Track");
  ntuple->column("trk_delPhGhost",delPhGhost, "trk_Nb",0 ,"Track");
  ntuple->column("trk_delPtGhost",delPtGhost, "trk_Nb",0 ,"Track");
  ntuple->column("trk_delPtLoop",delPtLoop, "trk_Nb",0 ,"Track");
  ntuple->column("trk_delPhLoopOpp",delPhLoopOpp, "trk_Nb",0 ,"Track");
  ntuple->column("trk_delThLoopOpp",delThLoopOpp, "trk_Nb",0 ,"Track");
  ntuple->column("trk_delPhLoopSame",delPhLoopSame, "trk_Nb",0 ,"Track");
  ntuple->column("trk_delThLoopSame",delThLoopSame, "trk_Nb",0 ,"Track");
  //ntuple->column("trk_Lund",DataLund, "trk_Nb",0 ,"Track");
  if(_sigModeString!="data")
    {
      //note: for w = 1.234 +/- 0.567, we store 1234.567  
      //ntuple->column("trk_eWeight",eEff, "trk_Nb",0 ,"Track");
      //ntuple->column("trk_muTWeight",muTEff, "trk_Nb",0 ,"Track");
      //ntuple->column("trk_muVTWeight",muVTEff, "trk_Nb",0 ,"Track");
      //ntuple->column("trk_piWeight",piEff, "trk_Nb",0 ,"Track");
      //ntuple->column("trk_Eff",trkEff, "trk_Nb",0 ,"Track");

      ntuple->column("trk_eW",eW, "trk_Nb",0 ,"Track");
      ntuple->column("trk_eWE",eWE, "trk_Nb",0 ,"Track");
      ntuple->column("trk_muTW",muTW, "trk_Nb",0 ,"Track");
      ntuple->column("trk_muTWE",muTWE, "trk_Nb",0 ,"Track");
      ntuple->column("trk_muVTW",muVTW, "trk_Nb",0 ,"Track");
      ntuple->column("trk_muVTWE",muVTWE, "trk_Nb",0 ,"Track");
      ntuple->column("trk_piW",piW, "trk_Nb",0 ,"Track");
      ntuple->column("trk_piWE",piWE, "trk_Nb",0 ,"Track");
      ntuple->column("trk_W",trkW, "trk_Nb",0 ,"Track");
      ntuple->column("trk_WE",trkWE, "trk_Nb",0 ,"Track");      
    }

  if ( _lookMCTruth.value()==true ) 
    {
      ntuple->column("trk_TruInd",TruInd, "trk_Nb",0 ,"Track");
      ntuple->column("trk_TruIndOne",TruIndOne, "trk_Nb",0 ,"Track");
      ntuple->column("trk_TruIndTwo",TruIndTwo, "trk_Nb",0 ,"Track");
      ntuple->column("trk_TruIndThree",TruIndThree, "trk_Nb",0 ,"Track");
      ntuple->column("trk_TruCons",TruCons, "trk_Nb",0 ,"Track");
      ntuple->column("trk_TruConsOne",TruConsOne, "trk_Nb",0 ,"Track");
      ntuple->column("trk_TruConsTwo",TruConsTwo, "trk_Nb",0 ,"Track");
      ntuple->column("trk_TruConsThree",TruConsThree, "trk_Nb",0 ,"Track");
    }

  delete CMS;

  //voila! :-P
}




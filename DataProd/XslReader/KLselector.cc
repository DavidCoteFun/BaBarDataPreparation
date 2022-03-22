//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: KLselector.cc,v 1.1 2004/10/05 23:41:39 cote Exp $
//
// Description:
//      Class KLselector
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      PAPPAGALLO MARCO               
//
// INPUT: Flag = 1 ---> Neutral Candidate from CalorClusterNeutral List (EMC)
//        Flag = 2 ---> Neutral Candidate from NeutralHad List (IFR)
//
// OUTPUT:Likelihood[0]={1.1 EMC candidate with energy < 0.2 GeV
//                       1.2 EMC candidate with a matched IFR neutral object
//                       1.3 EMC candidate without a match IFR neutral object
//                       2.1 IFR candidate with hit layers < 2
//                       2.2 IFR candidate with a matched EMC neutral object
//                       2.3 IFR candidate without a matched EMC neutral object}
//
//        Likelihood[1]=likelihood value
//
//
// NOTE: Use BetaCoreTools V00-01-22 in case of incompatibility with BtaCandMap.hh
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "XslReader/KLselector.hh"

//-------------
// C Headers --
//-------------
#include <assert.h>

//---------------
// C++ Headers --
//---------------
#include <iostream>
#include <math.h>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

#include "AbsEnv/AbsEnv.hh"
#include "AbsEvent/AbsEvent.hh"

#include "Beta/BtaCandidate.hh"
#include "BetaMicroAdapter/BtaCalQual.hh"
#include "BetaMicroAdapter/BtaIfrQual.hh"
#include "BetaCoreTools/BtaCandMap2.hh"

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

//----------------
// Constructors --
//----------------

// in general, a module constructor should not do much.  The begin(job) or
// begin(run) members are better places to put initialization
//

void KLselector::Selector(AbsEvent* anEvent,int cRun, const BtaCandidate* Candidate, double likelihood[], int flag){

  BtaCandMap2 *map=Ifd<BtaCandMap2>::get(anEvent, IfdStrKey("NeutralHadEMCNeutralMap"));

  likelihood[0]=0.;
  likelihood[1]=-100.;

  if(flag==1){

    const BtaCalQual* CalQualNHEMC = Candidate->getMicroAdapter()->getCalQual(); 
    if(CalQualNHEMC == 0) {return;} 
    if(CalQualNHEMC->ecalEnergy() < 0.2){
      //cout << "Energy<0.2 GeV"<<endl;
      likelihood[0]=1.1;
      likelihood[1]=-100.;
      return;
    }
    
    //*************** NHEMC coming/not coming from NHIFR ******************
    bool match_ifr= 0;
    if (map !=0) {
      match_ifr = map->containsImage(*Candidate);
      if(match_ifr == 1) {
        //cout << "EMC candidate matches with IFR "<<endl;
        likelihood[0]=1.2;
      }
      else likelihood[0]=1.3;
    }
    else likelihood[0]=1.3;

    if(cRun<18190){
      likelihood[1]=like2000_emc(Candidate);
    }
    if(cRun>18190 && cRun<25281){
      likelihood[1]=like2001_emc(Candidate);
    }
    if(cRun>25281){
      likelihood[1]=like2002_emc(Candidate);
    }
    return;
  }

  if(flag==2){
    const BtaIfrQual* IfrQualNHIFR = Candidate->getMicroAdapter()->getIfrQual();
    if(IfrQualNHIFR == 0) {
      return;
    }
    if(IfrQualNHIFR->IfrLayHits() < 2 ) {
      //cout << "Hit Layers < 2 "<<endl;
      likelihood[0]=2.1;
      likelihood[1]=-100.;
      return;
    }

    //*************** NHIFR coming from NHEMC ******************
    bool match_emc= 0;
    if (map !=0) {
      match_emc = map->containsClonable(*Candidate);
      if(match_emc == 1) { 
        //cout << "IFR candidate match with EMC "<<endl;
        likelihood[0]=2.2;
        likelihood[1]=-100.;
        return;
      }
    }
    //*********************************************************

    if(cRun<18190){
      likelihood[0]=2.3;
      likelihood[1]=like2000_ifr(Candidate);
    }
    if(cRun>18190 && cRun<25281){
      likelihood[0]=2.3;
      likelihood[1]=like2001_ifr(Candidate);
    }
    if(cRun>25281){
      likelihood[0]=2.3;
      likelihood[1]=like2002_ifr(Candidate);
    }
    return;
  } 
}


double KLselector::like2000_emc(const BtaCandidate* btaNHEMC) {

  const BtaCalQual* CalQual = btaNHEMC->getMicroAdapter()->getCalQual();
   
  double bin_ecal[12]={0.2, 0.23, 0.265, 0.31, 0.36 ,0.415, 0.48, 0.56, 0.66, 0.8, 1.05, 4.};
   
  double bin_secmom[11] = {0., 0.0005, 0.0007, 0.0009, 0.0011, 0.0014, 0.0018, 0.0024, 0.0035, 0.006, 0.05};
   
  double bin_zern42[12]={0., 0.01 , 0.02, 0.026, 0.035, 0.043, 0.05, 0.062, 0.081, 0.118, 0.175, 0.5};

  //---------------- PDFs(signal)/PDFs(background) -----------------------    
  double rapp_ecal[11] = {0.301264, 0.603186, 0.454526, 0.623907, 0.898884,0.880019, 1.10583, 1.548, 1.55347, 1.68106, 1.33072 };

  double rapp_secmom[10] = {0.647, 0.41775, 0.480883, 0.502786, 0.523571, 1.15033, 1.6998, 2.30794, 2.19139, 1.39902};

  double rapp_zern42[11] = {0.660619, 0.382701, 0.625601, 0.578903, 0.857912, 0.839228, 0.88757, 1.35823, 2.02477, 1.90629, 1.45819};
  

  double like_ecal = -10.;
  double like_secmom = -10.;
  double like_zern42 = -10.;
  double likelihood = 0.;
  
  for(int i = 0; i<=10; i++) {
    if(CalQual->ecalEnergy()>= bin_ecal[i] && CalQual->ecalEnergy()< bin_ecal[i+1]){
      like_ecal = log(rapp_ecal[i]);
      break;
    }
  }
  if(CalQual->ecalEnergy()>=bin_ecal[11]){
    like_ecal = -10.;
  }

  for(int j = 0; j<=9; j++) {
    if(CalQual->secondMomentTP()>= bin_secmom[j] && CalQual->secondMomentTP() < bin_secmom[j+1]){
      like_secmom = log(rapp_secmom[j]);
      break;
    }
  }
  if(CalQual->secondMomentTP()>=bin_secmom[10]){
    like_secmom = -10.;
  }
  
  for(int k = 0; k<=10; k++) {
    if(CalQual->absZernike42() >= bin_zern42[k] && CalQual->absZernike42() < bin_zern42[k+1]){
      like_zern42 = log(rapp_zern42[k]);
      break;
    }
  }
  if(CalQual->absZernike42()>=bin_zern42[11]){
    like_zern42 = -10.;
  }

  likelihood = like_ecal + like_secmom + like_zern42;
  return likelihood;
  
}

double KLselector::like2001_emc(const BtaCandidate* btaNHEMC) {

  const BtaCalQual* CalQual = btaNHEMC->getMicroAdapter()->getCalQual();
   
  double bin_ecal[14]={0.2, 0.225, 0.255, 0.29, 0.33 ,0.375, 0.425, 0.485, 0.555, 0.64, 0.75, 0.92, 1.25, 4.};
   
  double bin_secmom[11] ={0., 0.0005, 0.0007, 0.0009, 0.0011, 0.0014, 0.0018, 0.0024, 0.0035, 0.006, 0.05};

   
  double bin_zern42[13]= {0., 0.01 , 0.02, 0.026, 0.035, 0.043, 0.05, 0.062, 0.081, 0.118, 0.175, 0.3, 0.5};

  //---------------- PDFs(signal)/PDFs(background) -----------------------    
  double rapp_ecal[13] = {0.397696, 0.549548, 0.297421, 0.597827, 0.74877, 0.839402, 1.00334, 1.32575, 1.61193, 1.59421, 1.61394, 1.59477, 0.867589};

  double rapp_secmom[10] = {0.585188, 0.453729, 0.446209, 0.532921, 0.666205, 0.690908, 1.70983, 2.13953, 2.30114, 1.57608};

  double rapp_zern42[12] = {0.393179, 0.569172, 0.629292, 0.611394, 0.662729, 0.825844, 1.20914, 1.45163, 1.88207, 1.94868, 1.46432, 1.23173};
  

  double like_ecal = -10.;
  double like_secmom = -10.;
  double like_zern42 = -10.;
  double likelihood = 0.;
  
  for(int i = 0; i<=12; i++) {
    if(CalQual->ecalEnergy()>= bin_ecal[i] && CalQual->ecalEnergy()< bin_ecal[i+1]){
      like_ecal = log(rapp_ecal[i]);
      break;
    }
  }
  if(CalQual->ecalEnergy()>=bin_ecal[13]){
    like_ecal = -10.;
  }

  for(int j = 0; j<=9; j++) {
    if(CalQual->secondMomentTP()>= bin_secmom[j] && CalQual->secondMomentTP() < bin_secmom[j+1]){
      like_secmom = log(rapp_secmom[j]);
      break;
    }
  }
  if(CalQual->secondMomentTP()>=bin_secmom[10]){
    like_secmom = -10.;
  }

  for(int k = 0; k<=11; k++) {
    if(CalQual->absZernike42() >= bin_zern42[k] && CalQual->absZernike42() < bin_zern42[k+1]){
      like_zern42 = log(rapp_zern42[k]);
      break;
    }
  }

  if(CalQual->absZernike42()>=bin_zern42[12]){
    like_zern42 = -10.;
  }
  
  likelihood = like_ecal + like_secmom + like_zern42;
  return likelihood;
      
}

double KLselector::like2002_emc(const BtaCandidate* btaNHEMC) {

  const BtaCalQual* CalQual = btaNHEMC->getMicroAdapter()->getCalQual();
 
  double bin_ecal[12]={0.2, 0.23, 0.265, 0.31, 0.36 ,0.415, 0.48, 0.56, 0.66, 0.8, 1.05, 4.};
   
  double bin_secmom[11] ={0., 0.0005, 0.0007, 0.0009, 0.0011, 0.0014, 0.0018, 0.0024,0.0035,0.006,0.05};
   
  double bin_zern42[12]= {0., 0.01 , 0.02, 0.026, 0.035, 0.043, 0.05, 0.062, 0.081, 0.118, 0.175, 0.5};

  //---------------- PDFs(signal)/PDFs(background) -----------------------    
  double rapp_ecal[11] = {0.433207, 0.43654, 0.679982, 0.577463, 0.893009, 1.19631, 1.16281,1.53605, 1.61174, 1.45411, 1.11085};

  double rapp_secmom[10] = {0.474654, 0.354145, 0.546827, 0.420432, 0.656281, 1.07523, 1.57887, 2.17165, 2.25741, 1.38324 };

  double rapp_zern42[11] = {0.310454, 0.514843, 0.438853, 0.634539, 0.78664, 0.939736, 1.18344,  1.53788, 2.1836, 1.79884, 1.43289};
  

  double like_ecal = -10.;
  double like_secmom = -10.;
  double like_zern42 = -10.;
  double likelihood = 0.;
  
  for(int i = 0; i<=10; i++) {
    if(CalQual->ecalEnergy() >= bin_ecal[i] && CalQual->ecalEnergy()< bin_ecal[i+1]){
      like_ecal = log(rapp_ecal[i]);
      break;
    }
  }
  if(CalQual->ecalEnergy()>=bin_ecal[11]){
    like_ecal = -10.;
  }

  for(int j = 0; j<=9; j++) {
    if(CalQual->secondMomentTP()>= bin_secmom[j] && CalQual->secondMomentTP() < bin_secmom[j+1]){
      like_secmom = log(rapp_secmom[j]);
      break;
    }
  } 
  if(CalQual->secondMomentTP()>=bin_secmom[10]){
    like_secmom = -10.;
  }

  for(int k = 0; k<=10; k++) {
    if(CalQual->absZernike42() >= bin_zern42[k] && CalQual->absZernike42() < bin_zern42[k+1]){
      like_zern42 = log(rapp_zern42[k]);
      break;
    }
  }
  if(CalQual->absZernike42()>=bin_zern42[11]){
    like_zern42 = -10.;
  }

  likelihood = like_ecal + like_secmom + like_zern42;
  return likelihood;
      
}




double KLselector::like2000_ifr(const BtaCandidate* btaNHIFR) {
  const BtaIfrQual* IfrQual = btaNHIFR->getMicroAdapter()->getIfrQual(); 
 
  double bin_first[7]={-1.5, 0.5, 1.5, 3.5, 6.5, 14.5, 19.5};
   
  double bin_last[7] ={0.5, 3.5, 5.5, 7.5, 10.5, 16.5, 19.5};
   
  double bin_multp[8]= {0., 4., 4.8, 5.5, 6.2, 7.5, 10., 17.};

  //---------------- PDFs(signal)/PDFs(background) ----------------------- 
  double rapp_first[6] = {0.434161, 0.547672, 0.979651, 1.76972, 2.28175, 0.617311};
  
  double rapp_last[6] = {0.447262, 0.877661, 0.824655, 1.12196, 2.19836, 0.581326};
  
  double rapp_multp[7] = {0.247768, 1.1735, 0.911682, 1.07593, 1.01256, 1.6752, 1.15174};


  double like_first = -10.;
  double like_last = -10.;
  double like_multp = -10.;
  double likelihood = 0.;
  
  for(int i = 0; i<=5; i++) {
    if(IfrQual->firstHit() > bin_first[i] && IfrQual->firstHit() < bin_first[i+1]){
      like_first = log(rapp_first[i]);
    }
  }

  for(int j = 0; j<=5; j++) {
    if(IfrQual->lastHit()  > bin_last[j] && IfrQual->lastHit()  < bin_last[j+1]){
      like_last = log(rapp_last[j]);
      break;
    }
  }
  if(IfrQual->lastHit()<bin_last[0]){
    like_last = -10.;
  }

  double multip = float(IfrQual->IfrNStrips())/float(IfrQual->IfrLayHits());
  for(int k = 0; k<=6; k++) {
    if(multip >= bin_multp[k] && multip < bin_multp[k+1]){
      like_multp = log(rapp_multp[k]);
      break;
    }
  }
  if(multip>=bin_multp[7]){
    like_multp = -10.;
  }
  
  likelihood = like_first + like_last + like_multp;
  return likelihood;
      
}

double KLselector::like2001_ifr(const BtaCandidate* btaNHIFR) {

  const BtaIfrQual* IfrQual = btaNHIFR->getMicroAdapter()->getIfrQual(); 

  double bin_first[9]={-1.5, 0.5, 1.5, 2.5, 3.5, 5.5 ,8.5, 13.5, 19.5};
   
  double bin_last[9] ={0.5, 3.5, 5.5, 7.5, 9.5, 11.5, 14.5, 17.5, 19.5};
   
  double bin_multp[9]= {0., 4., 4.6, 5.3, 6. ,6.6, 8., 10., 17};

  //---------------- PDFs(signal)/PDFs(background) ----------------------- 
  double rapp_first[8] = {1.13473, 0.342949, 0.788372, 1.45693, 2.05458, 1.651, 2.33409, 0.494942};
  
  double rapp_last[8] = {0.624435, 0.865978, 1.30984,  0.909373, 1.57373, 2.0794, 1.4868, 0.298958};
  
  double rapp_multp[8] = {0.246419, 0.86089, 1.0322, 1.14574, 0.963061, 1.42706, 1.67481, 1.54033};


  double like_first = -10.;
  double like_last = -10.;
  double like_multp = -10.;
  double likelihood = 0.;
  
  for(int i = 0; i<=7; i++) {
    if(IfrQual->firstHit() > bin_first[i] && IfrQual->firstHit() < bin_first[i+1]){
      like_first = log(rapp_first[i]);
      break;
    }
  }


  for(int j = 0; j<=7; j++) {
    if(IfrQual->lastHit()  > bin_last[j] && IfrQual->lastHit()  < bin_last[j+1]){
      like_last = log(rapp_last[j]);
      break;
    }
  }
  if(IfrQual->lastHit()<bin_last[0]){
    like_last = -10.;
  }

  double multip = float(IfrQual->IfrNStrips())/float(IfrQual->IfrLayHits());
  for(int k = 0; k<=7; k++) {
    if(multip >= bin_multp[k] && multip < bin_multp[k+1]){
      like_multp = log(rapp_multp[k]);
      break;
    }
  }
  if(multip>=bin_multp[8]){
    like_multp = -10.;
  }
  
  likelihood = like_first + like_last + like_multp;
  return likelihood;
      
}

double KLselector::like2002_ifr(const BtaCandidate* btaNHIFR) {

  const BtaIfrQual* IfrQual = btaNHIFR->getMicroAdapter()->getIfrQual(); 

  double bin_first[9]={-1.5, 0.5, 1.5, 2.5, 3.5, 5.5 ,8.5, 13.5, 19.5};
   
  double bin_last[9] ={0.5, 3.5, 5.5, 7.5, 9.5, 11.5, 14.5, 17.5, 19.5};
   
  double bin_multp[9]= {0., 4., 4.6, 5.3, 6. ,6.6, 8., 10., 17};

  //---------------- PDFs(signal)/PDFs(background) ----------------------- 
  double rapp_first[8] = {1.17399, 0.532221, 0.968757, 1.45615, 1.40816, 1.69632, 1.82307, 0.716372};

  double rapp_last[8] = {0.731041, 0.9189, 1.27461, 0.917349, 1.49564, 1.91641, 1.50251, 0.281553};

  double rapp_multp[8] = {0.320662, 1.25842, 1.20846, 1.69136, 1.39783, 1.34666, 1.52331, 0.840678};


  double like_first = -10.;
  double like_last = -10.;
  double like_multp = -10.;
  double likelihood = 0.;
  
  for(int i = 0; i<=7; i++) {
    if(IfrQual->firstHit() > bin_first[i] && IfrQual->firstHit() < bin_first[i+1]){
      like_first = log(rapp_first[i]);
      break;
    }
  }


  for(int j = 0; j<=7; j++) {
    if(IfrQual->lastHit()  > bin_last[j] && IfrQual->lastHit()  < bin_last[j+1]){
      like_last = log(rapp_last[j]);
      break;
    }
  }
  if(IfrQual->lastHit()<bin_last[0]){
    like_last = -10.;
  }

  double multip = float(IfrQual->IfrNStrips())/float(IfrQual->IfrLayHits());
  for(int k = 0; k<=7; k++) {
    if(multip >= bin_multp[k] && multip < bin_multp[k+1]){
      like_multp = log(rapp_multp[k]);
      break;
    }
  }
  if(multip>=bin_multp[8]){
    like_multp = -10.;
  }
  
  likelihood = like_first + like_last + like_multp;
  return likelihood;
      
}






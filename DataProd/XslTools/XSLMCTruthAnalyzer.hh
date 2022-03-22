//*****************************************************************
//   
//   Creation: David Cote, Universite de Montreal, 12/11/03
//

#ifndef XSLMCTRUTHANALYZER
#define XSLMCTRUTHANALYZER

#include "CLHEP/Alist/AList.h"
#include "CLHEP/Alist/AListBase.h"
#include "Beta/BtaCandidate.hh"



class XSLMCTruthAnalyzer {

public:

  //This dummy constructor is mandatory to use XSLMCTruthAnalyzer objects as AppModule class variables
  XSLMCTruthAnalyzer(HepAList<BtaCandidate> inputMcList);
  XSLMCTruthAnalyzer();
  ~XSLMCTruthAnalyzer();
  void AnalyzeEvent();

  int eventType(){ return _eventType; }
  int typeB1(){ return _typeB1; }
  int typeB2(){ return _typeB2; }
  bool IsBadXuDecayMode( BtaCandidate* Xu, int XuLund );
  std::string modeB1(){ return _modeB1; }
  std::string modeB2(){ return _modeB2; }

  void AnalyzeMissingMomentum();
  HepLorentzVector p4Neutrinos() { return _Nus; }
  int nbNeutrinos() { return _nNus; }
  HepLorentzVector p4KLongs() { return _KLongs; }
  int nbKLongs() { return _nKLongs; }
  HepLorentzVector p4Neutrons() { return _Neutrons; }
  int nbNeutrons() { return _nNeutrons; }
  void SetB1Family();
  void SetB2Family();

  BtaCandidate Ups4S(){ return _Ups4S; }
  BtaCandidate B1(){ return _B1; }
  BtaCandidate Lep1_Lab(){ return _Lep1_Lab; }
  BtaCandidate Nu1_Lab(){ return _Nu1_Lab; }
  BtaCandidate Xu1_Lab(){ return _Xu1_Lab; }
  BtaCandidate VDaug1_Lab(){ return _VDaug1_Lab; }
  BtaCandidate B2(){ return _B2; }
  BtaCandidate Lep2_Lab(){ return _Lep2_Lab; }
  BtaCandidate Nu2_Lab(){ return _Nu2_Lab; }
  BtaCandidate Xu2_Lab(){ return _Xu2_Lab; }
  BtaCandidate VDaug2_Lab(){ return _VDaug2_Lab; }

  //In case of 2 signal B (rare!), this is done only for B2, which is the signal one in signal collections
  BtaCandidate fille1Lab(){ return _fille1Lab; } //pi+ or Gam1
  BtaCandidate fille2Lab(){ return _fille2Lab; }//pi- or Gam2
  BtaCandidate filleX0Lab(){ return _filleX0Lab; }//pi0/eta/rho0

  BtaCandidate pfille1Lab(){ return _pfille1Lab; } //pi+ or Gam1
  BtaCandidate pfille2Lab(){ return _pfille2Lab; }//pi- or Gam2
  BtaCandidate pfillePi0Lab(){ return _pfillePi0Lab; }

  BtaCandidate ppfilleGam1Lab(){ return _ppfilleGam1Lab; }
  BtaCandidate ppfilleGam2Lab(){ return _ppfilleGam2Lab; }


private:
  int BType( BtaCandidate* B, BtaCandidate &XuLab, BtaCandidate &LepLab, BtaCandidate &NuLab );
  int AnalyzeTau( BtaCandidate* tau );
  int AnalyzeD( BtaCandidate* D );
  int AnalyzeD0( BtaCandidate* D0 );
  int AnalyzeDplus( BtaCandidate* Dplus );
  int AnalyzeBtoHadrons( BtaCandidate* B );
  int AnalyzeBtoXulnu( BtaCandidate* B, BtaCandidate &XuLab, BtaCandidate &LepLab, BtaCandidate &NuLab, int type);
  int AnalyzeBtoXclnu( BtaCandidate* B, int type );
  std::string FigureOutXuFamily(BtaCandidate Xu, int type);     
  BtaCandidate GetXuDaughter(BtaCandidate XuLab, std::string mode); //GetXuDaughter needs FigureOutXuFamily first

  HepAList<BtaCandidate> _mcList;
  BtaCandidate _Ups4S;
  BtaCandidate _B1;
  BtaCandidate _Lep1_Lab;
  BtaCandidate _Nu1_Lab;
  BtaCandidate _Xu1_Lab;
  BtaCandidate _VDaug1_Lab;
  BtaCandidate _B2;
  BtaCandidate _Lep2_Lab;
  BtaCandidate _Nu2_Lab;
  BtaCandidate _Xu2_Lab;
  BtaCandidate _VDaug2_Lab;

  BtaCandidate _fille1Lab;
  BtaCandidate _fille2Lab;
  BtaCandidate _filleX0Lab;

  BtaCandidate _pfille1Lab;
  BtaCandidate _pfille2Lab;
  BtaCandidate _pfillePi0Lab;

  BtaCandidate _ppfilleGam1Lab;
  BtaCandidate _ppfilleGam2Lab;

  int _eventType;
  int _typeB1;
  int _typeB2;
  
  std::string _modeB1;
  std::string _modeB2;

  //From AnalyzeMissingMomentum
  HepLorentzVector _Nus;
  HepLorentzVector _KLongs;
  HepLorentzVector _Neutrons;
  int _nNus;
  int _nKLongs; 
  int _nNeutrons; 

};

#endif





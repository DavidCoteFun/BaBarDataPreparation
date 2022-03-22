//*****************************************************************
//   
//   Creation: David Cote, Universite de Montreal, 02/16/04
//

#ifndef XSLYRECOYANALYZER
#define XSLYRECOYANALYZER

#include "CLHEP/Alist/AList.h"
#include "Beta/BtaCandidate.hh"

class HepLorentzVector;
class BtaBooster;
class EventInfo;

class XSLRecoYAnalyzer {

public:

  enum WhichDaughter{NotApplicable, HigherEnergy, LowerEnergy};

  XSLRecoYAnalyzer();
  ~XSLRecoYAnalyzer();

  void Init( BtaCandidate* YLab );

  BtaCandidate YLab(){ return _YLab; }
  BtaCandidate XuLab(){ return _XuLab; }
  BtaCandidate LepLab(){ return _LepLab; }
  BtaCandidate YUps(){ return _YUps; }
  BtaCandidate XuUps(){ return _XuUps; }
  BtaCandidate LepUps(){ return _LepUps; }
  BtaCandidate VDLab();
  BtaCandidate VDUps();

  std::string mode(){ return _mode; }

  void FigureOutXuFamily();
  BtaCandidate fille1Lab(){ return _fille1Lab; } //pi+ or Gam1
  BtaCandidate fille2Lab(){ return _fille2Lab; } //pi- or Gam2
  BtaCandidate filleX0Lab(){ return _filleX0Lab; } //pi0/eta/rho0
  BtaCandidate filleX0Ups(){ return _filleX0Ups; } //pi0/eta/rho0

  BtaCandidate pfille1Lab(){ return _pfille1Lab; } //pi+ or Gam1
  BtaCandidate pfille2Lab(){ return _pfille2Lab; } //pi- or Gam2
  BtaCandidate pfillePi0Lab(){ return _pfillePi0Lab; }
  BtaCandidate pfillePi0Ups(){ return _pfillePi0Ups; }

  BtaCandidate ppfilleGam1Lab(){ return _ppfilleGam1Lab; }
  BtaCandidate ppfilleGam2Lab(){ return _ppfilleGam2Lab; }


  double cosBY();
  double pBCM();

  double openAngleL0(); //This means Xu vs Lep opening angle in Ups frame

  double HelicityAngleL1(); //Helicity angle of the composite daughter of Xu in Xu frame
  double HelicityAngleL2(); //Helicity angle of the composite daughter of Xu's daughter in Xu's daughter frame
  double HelicityAngleL3(); //Helicity angle of the composite daughter of Xu's grand-daughter in Xu's grand-daughter frame

  double mXuLep2(); //For Dalitz plot (with q2)
  void DalitzL1( double &m12, double &m23 ); //For Dalitz plot of eta/eta'/omega -> 3 bodys

  BtaCandidate fittedY_Tree(EventInfo* evtInfo, BtaCandidate &fittedXu, BtaCandidate &fittedLep, 
			    BtaCandidate &f1, BtaCandidate &f2,  BtaCandidate &fX0);
  BtaCandidate fittedB_Tree(EventInfo* evtInfo, Hep3Vector p3NuLab,  BtaCandidate &fittedNu, BtaCandidate &fittedXu, BtaCandidate &fittedLep,
			    BtaCandidate &f1, BtaCandidate &f2,  BtaCandidate &fX0);
  BtaCandidate fittedXu(EventInfo* evtInfo,  BtaCandidate &f1, BtaCandidate &f2,  BtaCandidate &fX0);
  BtaCandidate fittedEtaL1(EventInfo* evtInfo);
  double ProbChi2TreeToCutOn(EventInfo* evtInfo);
  double ProbChi2Y_Fast();

  BtaCandidate VDLabFromFittedDaughters(BtaCandidate f1, BtaCandidate f2,  BtaCandidate fX0);

  BtaCandidate* XuToFit(EventInfo* evtInfo);
  BtaCandidate* PiToFit(EventInfo* evtInfo);
  BtaCandidate* Pi0ToFit(EventInfo* evtInfo);
  BtaCandidate* Eta2ToFit(EventInfo* evtInfo);
  BtaCandidate* Eta3ToFit(EventInfo* evtInfo);
  BtaCandidate* EtapRGToFit(EventInfo* evtInfo);
  BtaCandidate* EtapE3PPToFit(EventInfo* evtInfo);
  BtaCandidate* EtapE2PPToFit(EventInfo* evtInfo);
  BtaCandidate* RhoCToFit(EventInfo* evtInfo);
  BtaCandidate* Rho0ToFit(EventInfo* evtInfo);
  BtaCandidate* OmegaToFit(EventInfo* evtInfo);

  double deltaE( HepLorentzVector PmissUps );
  double deltaE2( HepLorentzVector PmissUps );
  double delThetaNu( HepLorentzVector PmissUps );
  double mES( HepLorentzVector PmissUps );
  

private:
  void GetVD();
  BtaCandidate GetXu( BtaCandidate Y );
  BtaCandidate GetLep( BtaCandidate Y );

  HepLorentzVector _Ups;
  bool _FamilyFiguredOut;

  double _mB;
  double _eBeamCM;
  double _f;
  BtaCandidate _YLab;
  BtaCandidate _LepLab;
  BtaCandidate _XuLab;
  BtaCandidate _YUps;
  BtaCandidate _LepUps;
  BtaCandidate _XuUps;

  BtaCandidate _VDLab;
  BtaCandidate _VDUps;

  BtaCandidate _fille1Lab;
  BtaCandidate _fille2Lab;
  BtaCandidate _filleX0Lab;
  BtaCandidate _filleX0Ups;

  BtaCandidate _filleX0LabForFit;

  BtaCandidate _pfille1Lab;
  BtaCandidate _pfille2Lab;
  BtaCandidate _pfillePi0Lab;
  BtaCandidate _pfillePi0Ups;

  BtaCandidate _ppfilleGam1Lab;
  BtaCandidate _ppfilleGam2Lab;

  std::string _mode;
  XSLRecoYAnalyzer::WhichDaughter _DaughterToPick;

};

#endif





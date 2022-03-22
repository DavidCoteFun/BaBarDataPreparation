//*****************************************************************
//   
//   Creation: David Cote, Universite de Montreal, 12/11/03
//   see comments in the XSLKin.hh file
#include "BaBar/BaBar.hh"

#include "XslFFReweighting/XSLKin.hh"


XSLKin::XSLKin( XSLKin* Kin ) :
  _BLab(Kin->BLab()),
  _LepLab(Kin->LepLab()),
  _XuLab(Kin->XuLab()),
  _XuDaugLab(Kin->XuDLab()),
  _isVector(Kin->isVector())
{ 
  Init(); 
  Compute();
}

XSLKin::XSLKin(HepLorentzVector BLab, HepLorentzVector LepLab, HepLorentzVector XuLab) :
  _BLab(BLab),
  _LepLab(LepLab),
  _XuLab(XuLab)
{ 
  _isVector=false;
  _XuDaugLab=HepLorentzVector(0,0,0,0);
  Init(); 
  Compute();
}


XSLKin::XSLKin(HepLorentzVector BLab, HepLorentzVector LepLab, HepLorentzVector XuLab, HepLorentzVector XuDaughterLab) :
  _BLab(BLab),
  _LepLab(LepLab),
  _XuLab(XuLab),
  _XuDaugLab(XuDaughterLab)
{ 
  _isVector=true;
  Init(); 
  Compute();
}

void
XSLKin::Init()
{
  _ctl=-666;
  _theta_l=-666;
  _q2=-666;
  _w=-666;
  _ctv=-666;
  _theta_v=-666;
  _chi=-666;
}


void 
XSLKin::Compute()
{
  //Xu meson in the B frame
  HepLorentzVector XuB = _XuLab;
  XuB.boost(-_BLab.boostVector());

  //q2 computed as (B-Xu)^2 in the B frame
  _q2=_BLab.m2() + _XuLab.m2() - 2*_BLab.m()*XuB.e();
  _w=(_BLab.m2()+ _XuLab.m2()-_q2)/(2*_BLab.m()*_XuLab.m());

  //W boson in the Lab frame
  HepLorentzVector WLab = _BLab - _XuLab;
  if(WLab.mag2()<=0)
    {
      //this non-physical case is identified in the ntuple by q2<0 
      //To compute theta_l and chi however, we'll arbitrarly set the W mass to "almost zero" to avoid nan in the ntuple
      WLab.setVectM(WLab.vect(),0.000001);   
    }
  HepLorentzVector W_B = WLab;
  W_B.boost(-_BLab.boostVector());

  //Lepton in the W frame
  HepLorentzVector LepW = _LepLab;
  LepW.boost(-_BLab.boostVector());
  LepW.boost(-W_B.boostVector());

  _ctl=LepW.vect().unit()*W_B.vect().unit();
  _theta_l=LepW.vect().angle(W_B.vect());


  if(_isVector){

    //Daugther of the Xu meson in Xu meson's frame
    HepLorentzVector VD_Xu(_XuDaugLab);
    VD_Xu.boost(-_BLab.boostVector());
    VD_Xu.boost(-XuB.boostVector());
    
    //ctv
    _ctv=VD_Xu.vect().unit()*XuB.vect().unit();
    _theta_v=VD_Xu.vect().angle(XuB.vect());
    
    //chi      
    //Art Snyder's (correct) 
    Hep3Vector p3LepTrans=LepW.vect()-_ctl*LepW.rho()*W_B.vect().unit();
    Hep3Vector p3VDTrans=VD_Xu.vect()-_ctv*VD_Xu.rho()*XuB.vect().unit();
    _chi=p3VDTrans.angle(p3LepTrans);
    
    //Leif Wilden's  (wrong...)
      
    Hep3Vector LepPlane = LepW.vect().cross(W_B.vect());
    Hep3Vector HadPlane = VD_Xu.vect().cross(XuB.vect());
    _chib = LepPlane.angle(HadPlane);
    

  }//if(_isVector)

  return;
}







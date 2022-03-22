//*****************************************************************
//   
//   Creation: David Cote, Universite de Montreal, 12/11/03
//   see comments in the XSLEvtFFWeight.hh file or BAD #809
//
//   Calculation from EvtGenModels/EvtISGW2.cc 
#include "BaBar/BaBar.hh"

#include "XslFFReweighting/XSLVectorISGW2.hh"
#include "CLHEP/Vector/LorentzVector.h"
#include <math.h>
using std::cout;
using std::endl;
using std::string;

XSLVectorISGW2::XSLVectorISGW2( double mB, double mXu, double q2, double theta_l, double theta_V, double chi, string mode) : 
  XSLVectorFF(mB,mXu,q2,theta_l,theta_V,chi) 
{  SetNormalizations(mode);  Compute(); }

XSLVectorISGW2::XSLVectorISGW2( XSLKin* DecayKin, string mode) : XSLVectorFF(DecayKin) 
{  SetNormalizations(mode);  Compute(); }

XSLVectorISGW2::XSLVectorISGW2(HepLorentzVector BLab, HepLorentzVector LepLab, HepLorentzVector XuLab, HepLorentzVector XuDaughterLab, string mode) : 
  XSLVectorFF(BLab, LepLab, XuLab, XuDaughterLab) 
{  SetNormalizations(mode);  Compute(); }


void XSLVectorISGW2::GetAllFF(double *A1, double *A2, double *V)
{ 
  GetAllISGW2FF(A1,A2,V,_q2); 
  return; 
}

double XSLVectorISGW2::GetA1(double q2)
{
  double A1=0,A2=0,V=0;
  GetAllISGW2FF(&A1,&A2,&V,q2); 

  return A1;
}

double XSLVectorISGW2::GetA2(double q2)
{
  double A1=0,A2=0,V=0;
  GetAllISGW2FF(&A1,&A2,&V,q2); 

  return A2;
}

double XSLVectorISGW2::GetV(double q2)
{
  double A1=0,A2=0,V=0;
  GetAllISGW2FF(&A1,&A2,&V,q2); 

  return V;
}

void XSLVectorISGW2::GetAllISGW2FF(double *A1, double *A2, double *V, double q2)
{
	double f,g,ap;
	double t=q2;
	double mass=_mV;
	double mb=_mB;
	///////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//	   if (daught==DST0||daught==DSTP||daught==DSTM||daught==DSTB||
//       daught==OMEG||daught==RHO0||daught==RHOM||daught==RHOP||
//       daught==KSTP||daught==KSTM||daught==KST0||daught==KSTB||
//       daught==PHI||daught==DSSTP||daught==DSSTM) {
//     EvtISGW2FF3S1(parent,daught,t,mass,&ff,&gf,&apf,&amf);

	
	//void  EvtISGW2FF::EvtISGW2FF3S1(EvtId parent,EvtId daugt,double t,double mass,
      //double *f,double *g,double *ap,double *am){


  double cf,mtb,wt,msd,mup,f3f,msq,bb2,mum,mtx,bbx2,f3g;
  double cji,bx2,f3appam,msb,tm,mbb,mbx;
  double f3apmam,appam,apmam,mx,f3;
  double r_f,r_g,r_appam,r_apmam, betaji_f,betaji_g;
  double betaji_appam, betaji_apmam;
  double w,mqm,r2,chiji,zji,ai,al,rcji,nf,nfp,gammaji;

  //For B0/B+ -> Rho/Omega
    msb=5.2;
    msd=0.33;
    bb2=0.431*0.431;
    mbb=5.31;
    nf = 4.0;
	cf=0.905;
	msq=0.33;
	bx2=0.299*0.299;
	mbx=0.75*0.770+0.25*0.14;
	nfp = 0.0;

	//OK, here we go!

  mtb=msb+msd;
  mtx=msq+msd;
  
  mup=1.0/(1.0/msq+1.0/msb);
  mum=1.0/(1.0/msq-1.0/msb);
  bbx2=0.5*(bb2+bx2);
  mx=mass;
  tm=(mb-mx)*(mb-mx);
  if ( t > tm ) t = 0.99*tm;

  wt=1.0+(tm-t)/(2.0*mbb*mbx);
  mqm = 0.1;
  
  r2=3.0/(4.0*msb*msq)+3*msd*msd/(2*mbb*mbx*bbx2) + 
    (16.0/(mbb*mbx*(33.0-2.0*nfp)))*
    log(EvtGetas(mqm,mqm)/EvtGetas(msq,msq));

  w = 1.0 + (( tm - t ) / ( 2.0* mb * mx ));

  rcji = ( 1/sqrt(w*w -1 ))*log( w + sqrt( w*w -1 ));
  
  al = (8.0 / ( 33.0 - 2.0*nfp ))*(w*rcji -1.0 );

  ai = -1.0* ( 6.0/( 33.0 - 2.0*nf));  
  
  cji = pow(( EvtGetas( msb,msb ) / EvtGetas( msq,msq ) ),ai);
  zji = msq / msb;

  gammaji = EvtGetGammaji( zji );

  chiji = -1.0 - ( gammaji / ( 1- zji ));
  
  betaji_g = (2.0/3.0)+gammaji;
  betaji_f = (-2.0/3.0)+gammaji;
  betaji_appam = -1.0-chiji+(4.0/(3.0*(1.0-zji)))+
                 (2.0*(1+zji)*gammaji/(3.0*(1.0-zji)*(1.0-zji)));
  
  betaji_apmam = (1.0/3.0)-chiji-(4.0/(3.0*(1.0-zji)))-
                 (2.0*(1+zji)*gammaji/(3.0*(1.0-zji)*(1.0-zji)))+
                 gammaji;

  r_g = cji*(1+(betaji_g*EvtGetas( msq,sqrt(mb*msq) )/PI));
  r_f = cji*(1+(betaji_f*EvtGetas( msq,sqrt(mb*msq) )/PI));
  r_appam = cji*(1+(betaji_appam*EvtGetas( msq,sqrt(mb*msq) )/PI));
  r_apmam = cji*(1+(betaji_apmam*EvtGetas( msq,sqrt(mb*msq) )/PI));

  
  f3=sqrt(mtx/mtb)*pow(sqrt(bx2*bb2)/bbx2,1.5)/
    ((1.0+r2*(tm-t)/12.0)*(1.0+r2*(tm-t)/12.0));
  
  f3f=sqrt(mbx*mbb/(mtx*mtb))*f3;
  f3g=sqrt(mtx*mtb/(mbx*mbb))*f3;
  f3appam=sqrt(mtb*mtb*mtb*mbx/(mbb*mbb*mbb*mtx))*f3;
  f3apmam=sqrt(mtx*mtb/(mbx*mbb))*f3;
  f=cf*mtb*(1+wt+msd*(wt-1)/(2*mup))*f3f*r_f;
  g=0.5*(1/msq-msd*bb2/(2*mum*mtx*bbx2))*f3g*r_g;
  
  appam=cji*(msd*bx2*(1-msd*bx2/(2*mtb*bbx2))/ 
	     ((1+wt)*msq*msb*bbx2)-
	     betaji_appam*EvtGetas( msq,sqrt(msq*mb) )/
	     (mtb*PI))*f3appam;
  
  apmam=-1.0*(mtb/msb-msd*bx2/(2*mup*bbx2)+wt*msd*mtb*bx2*
	      (1-msd*bx2/(2*mtb*bbx2))/((wt+1)*msq*msb*bbx2))*
    f3apmam*r_apmam/mtx;
  
  ap=0.5*(appam+apmam);

	///////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////

  *V = (g)*(mb+mass);
  *A1 = (f)/(mb+mass);
  *A2 = -1.0*(ap)*(mb+mass);
  

}


//////////////////////////////////////////////////////
//Utilities from EvtGen
double XSLVectorISGW2::EvtGetGammaji ( double z )
{
double temp;

   temp = 2+((2.0*z)/(1-z))*log(z);
   temp = -1.0*temp;

   return temp;

} //EvtGetGammaji



double XSLVectorISGW2::EvtGetas ( double massq, double massx )
{
  double lqcd2 = 0.04;
  double nflav = 4;
  double temp = 0.6;
  
  if ( massx > 0.6 ) {
    if ( massq < 1.85 ) {
      nflav = 3.0;}
    
    temp = 12.0*PI / ( 33.0 - 2.0*nflav) /
      log( massx*massx/lqcd2);
  }
  return temp;
  
} //EvtGetas

double XSLVectorISGW2::EvtGetas ( double mass )
     
{
  double lqcd2 = 0.04;
  double nflav = 4;
  double temp = 0.6;

  if ( mass > 0.6 ) {
    if ( mass < 1.85 ) {
      nflav = 3.0;}
    
    temp = 12.0*PI / ( 33.0 - 2.0*nflav) /
      log( mass*mass/lqcd2);
  }
  return temp;
  
} //EvtGetas


void
XSLVectorISGW2::SetNormalizations(string mode)
{
  //Note: These normalization constants were determined empirically from MC samples of 1000k events.

  if(mode=="Ignore") { _FLATQ2Normalization=1.0; _PHSPNormalization=1.0; _ISGW2Normalization=1.0; }
  else if(mode=="rhoClnu") {
    _FLATQ2Normalization=0.1098110339;  //determined from 1000k events: 0.1% error stat
    _PHSPNormalization=0.1650473868;   //determined from 1000k events: 0.1% error stat
    _ISGW2Normalization=1.0; }  //1.0 by construction
  else if(mode=="rho0lnu") {
    _FLATQ2Normalization=0.1100517287;  //determined from 1000k events: 0.1% error stat
    _PHSPNormalization=0.1655613652;  //determined from 1000k events: 0.1% error stat
    _ISGW2Normalization=1.0; }  //1.0 by construction
  else if(mode=="omegalnu" || mode=="omegalnu_o") { 
    _FLATQ2Normalization=0.1123092638; //determined from 1000k events: 0.1% error stat
    _PHSPNormalization=0.1697422693; //determined from 1000k events: 0.1% error stat
    _ISGW2Normalization=1.0; }  //1.0 by construction
  else if(mode=="pilnu"||mode=="pi0lnu"||mode=="etalnu"||mode=="eta2lnu"||mode=="eta3lnu"||
	  mode=="etaplnu"||mode=="etaplnuE2PP"||mode=="etaplnuE3PP"||mode=="etaplnuRG") { 
    cout<<"XSLVectorISGW2::SetNormalizations  -- This model is not intended to be used with mode: "<<mode<<" !!!!"<<endl;
    cout<<"You probably want to use XSLPseudoScalarISGW2 instead... We still fix the normalizations to 1.0, but these results might be meaningless."<<endl;
    _FLATQ2Normalization=1.0; _PHSPNormalization=1.0; _ISGW2Normalization=1.0;  }
  else {
    cout<<"XSLVectorISGW2::SetNormalizations  -- Unknown mode: "<<mode<<" !!!    Normalizations set to 1..."<<endl;
    _FLATQ2Normalization=1.0; _PHSPNormalization=1.0; _ISGW2Normalization=1.0;  }
  
  return;
}

















// su2_hb_update.cc, 2004/11/16 velytsky
#include "chromabase.h"
#include "util/gauge/su2extract.h"
#include "util/gauge/sunfill.h"

namespace Chroma 
{

  void print_field(const LatticeReal& a0);
  void su2_a_0(const LatticeReal&, LatticeReal& ,
	       const Subset& sub, const int NmaxHB,
	       LatticeBoolean& lAccept);
  void su2_a_0_kp(const LatticeReal&, LatticeReal& ,
		  const Subset& sub, const int NmaxHB,
		  LatticeBoolean& lAccept);

  void su2_hb_update(LatticeColorMatrix& u_mu, const LatticeColorMatrix& 
		     u_mu_staple, Double BetaMC, const int su2_index,
		     const Subset& sub, const int NmaxHB) 
  {
    // ****************************************** 
    //              Parameters
    Double fuzz=1e-16; //the smallest number for devision (det)
    // ******************************************
    //              Temp Storages
    LatticeBoolean lbtmp; // storage for lattice booleans
    //LatticeReal lftmp; // storage for lattice float tmps
    // ******************************************
    //V=U*U_staple
    LatticeColorMatrix v;
    v[sub]=u_mu*u_mu_staple;
    // convert to Pauli matrices parametrization
    LatticeReal a_0;
    multi1d<LatticeReal> r(4);
    su2Extract(r,v,su2_index,sub);
    //compensate for extra 2 !!! su(2)
    Real half(0.5);
    r[0][sub] *= half;
    r[1][sub] *= half;
    r[2][sub] *= half;
    r[3][sub] *= half;
    LatticeReal SqDet;
    //SqDet=-1;
    SqDet[sub]=sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]+r[3]*r[3]);
    //print_field(SqDet);
    // Normalize
    //lftmp[sub]=1.0/where(lbtmp,SqDet,LatticeReal(1)); //not needed???
    // Inverse matrix u^-1
    r[0][sub]=r[0]/SqDet;
    r[1][sub]=-r[1]/SqDet;
    r[2][sub]=-r[2]/SqDet;
    r[3][sub]=-r[3]/SqDet;
		
    LatticeBoolean lAccept;
    su2_a_0(BetaMC*SqDet, a_0, sub, NmaxHB,lAccept);
    //su2_a_0_kp(BetaMC*SqDet, a_0, sub, NmaxHB,lAccept);
    //print_field(a_0);
    multi1d<LatticeReal> a(4);
    a[0][sub]=a_0;
    //other a components
    LatticeReal CosTheta, Phi;
    LatticeReal a_abs;
    LatticeInt lWarning;
    int iWarning;
    a_abs[sub] = 1.0-a[0]*a[0];
    lbtmp = (1>0);
    lbtmp[sub] = (a_abs >= 0);
    lWarning=where(lbtmp,0,1);
    iWarning=toInt(sum(lWarning));
    if(iWarning>0) QDPIO::cerr <<"wrong a_0!!!"<<endl;
    lbtmp = (1>0);
    lbtmp[sub] = (a_abs > fuzz);
    lWarning=where(lbtmp,0,1);
    iWarning=toInt(sum(lWarning));
    if(iWarning>0) QDPIO::cerr <<"large a_0!!!"<<endl;
    LatticeReal a_r;
    a_r[sub]=sqrt(a_abs);
    random(CosTheta,sub);
    Real RDummy;
    random(RDummy);
    CosTheta[sub]=1.0-2.0*CosTheta;
    a[3][sub]=a_r*CosTheta;
    //LatticeReal pr_a;
    //pr_a[sub]=a[3];
    //print_field(pr_a);
    CosTheta[sub]=(1-CosTheta*CosTheta);
    CosTheta[sub]=sqrt(CosTheta);//SinTheta
    random(Phi,sub);
    random(RDummy);
    Phi[sub]*=8.0*atan(1.0);
    a_r[sub]*=CosTheta; //a_r*SinTheta
    a[1][sub]=a_r*cos(Phi);
    a[2][sub]=a_r*sin(Phi);
    //u'=uu^-1 -> b = a*r
    multi1d<LatticeReal> b(4);
    b[0][sub]=a[0]*r[0]-a[1]*r[1]-a[2]*r[2]-a[3]*r[3];
    b[1][sub]=a[0]*r[1]+a[1]*r[0]-a[2]*r[3]+a[3]*r[2];
    b[2][sub]=a[0]*r[2]+a[2]*r[0]-a[3]*r[1]+a[1]*r[3];
    b[3][sub]=a[0]*r[3]+a[3]*r[0]-a[1]*r[2]+a[2]*r[1];
    sunFill(v,b,su2_index,sub);
    u_mu[sub]=where(lAccept,v*u_mu,u_mu);

  }

  void su2_a_0(const LatticeReal& weight, LatticeReal& a_0, 
	       const Subset& sub,const int NmaxHB,
	       LatticeBoolean& lAccept) 
  {
    // received weight=SqDet*BetaMC
    int vol_cb = Layout::vol() >> 1; //volume of cb sublattice
    int vol_accept;
    LatticeInt ilbtmp=0;
    lAccept = (1 < 0);
    LatticeReal w_exp;
    w_exp[sub]=exp(-2.0*weight); //too small, need to avoid
    LatticeReal x; //container for random numbers
    int n_runs=0;
    Real RDummy;
    do {
      n_runs++;
      //random(x[sub]);
      random(x,sub);
      random(RDummy);
      //a_0[sub]=where(lAccept,a_0,1.0+log(x*(1.0-w_exp)+w_exp)/weight);
      a_0[sub]=where(lAccept,a_0,1.0+log(w_exp*(1-x)+x)/weight);
      //print_field(a_0);exit(1);
      //random(x[sub]);
      random(x,sub);
      random(RDummy);
      //lAccept[sub] = where(lAccept,(1 > 0),((x*x) < (1.0-a_0*a_0)));
      //x=1.0l-x;
      x=x*x;
      LatticeReal a_0trial;
      a_0trial=a_0*a_0;
      a_0trial=1-a_0trial;
      lAccept[sub] = where(lAccept,(1 > 0),(x < a_0trial));
      //lAccept[sub] = where(lAccept,(1 > 0),(x > (1.0-sqrt(1.0-a_0*a_0))));
      ilbtmp[sub]=where(lAccept,1,0); //convert to 1/0
      vol_accept = toInt(sum(ilbtmp));
    } while ((vol_accept < vol_cb) && ((NmaxHB <= 0) || (n_runs<NmaxHB)));
  }

  void su2_a_0_kp(const LatticeReal& weight, LatticeReal& a_0,
		  const Subset& sub, const int NmaxHB,
		  LatticeBoolean& lAccept) 
  {
    // received weight=SqDet*BetaMC
    int vol_cb = Layout::vol() >> 1; //volume of cb sublattice
    lAccept = (1 < 0);
    LatticeInt ilbtmp=0;
    int vol_accept;
    LatticeReal xr1,xr2,xr3,xr4;
    int n_runs=0;
    do {
      n_runs++;
      random(xr1);
      random(xr2);
      random(xr3);
      random(xr4);
      xr1=-(log(xr1)/weight);
      xr3 = cos(2.0l*M_PI*xr3);
      xr3=xr3*xr3;
      xr2=-(log(xr2)/weight);
      //a_0=xr2+xr1*xr3;
      a_0[sub]=where(lAccept,a_0,xr2+xr1*xr3);
      lAccept[sub] = where(lAccept,(1 > 0),(xr4*xr4) < (1.0l-a_0/2.0l));
      ilbtmp[sub]=where(lAccept,1,0); //convert to 1/0
      vol_accept = toInt(sum(ilbtmp));
      //QDPIO::cout<<vol_accept<<endl;

    } while ((vol_accept < vol_cb) && ((NmaxHB <= 0) || (n_runs<NmaxHB)));
    a_0=1.0l-a_0;
  }

  void print_field(const LatticeReal& a0) 
  {
    // single node only !!!
    double temp,fuzz=1e-44;
    const int nodeSites = Layout::sitesOnNode();
    const int nodeNumber = Layout::nodeNumber();
    for(int i=0;i<nodeSites;i++) {
      temp=toDouble(peekSite(a0,Layout::siteCoords(nodeNumber,i)));
      if(temp > fuzz || temp<-fuzz) cout<<temp<<" ";
      //if(temp > -0.1 && temp <0.01) cout<<temp<<" ";
      //cout<<temp<<" ";
    }
    cout<<endl;

  }

}  // end namespace Chroma

// $Id: overlapbu_linop_w.cc,v 1.3 2003-10-09 21:06:29 edwards Exp $
/*! \file
 *  \brief Overlap operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/overlapbu_linop_w.h"

void zol_const(Double& qq1, Double qq[], Double cq[]);
static void wldirac(const LatticeFermion& psi, 
		    LatticeFermion& chi,
		    const multi1d<LatticeColorMatrix>& u, 
		    const Real& Kappa,
		    int   isign);



//! Creation routine
/*! \ingroup fermact
 *
 * \param _u 	      gauge field     	       (Read)
 * \param _OverMass   Mass within Wilson operator kernel (Read)
 * \param _m_q        fermion kappa   	       (Read)
 */
void OverlapBULinOp::create(const multi1d<LatticeColorMatrix>& _u, const Real& _OverMass, const Real& _m_q)
{
  m_q = _m_q;
  OverMass = _OverMass;
  u = _u;
}


//! Apply Overlap operator
/*!
 * \ingroup linop
 *
 * The operator acts on the entire lattice
 *
 * \param psi 	  Pseudofermion field     	       (Read)
 * \param isign   Flag ( PLUS | MINUS )   	       (Read)
 */
LatticeFermion OverlapBULinOp::operator() (const LatticeFermion& psi, enum LinOpSign isign) const
{
  LatticeFermion chi;

  START_CODE("OverlapBULinOp");

  const int n_s=14;
  int  n_count, MaxCG;
  Double  qq[n_s],cq[n_s],qq1;
  Real    kappa = 0.5 / (Nd - OverMass);
  LatticeFermion v0,v1,v2,v3,sol,x[n_s],p[n_s],r[n_s],r_aux;
  DComplex gamma[n_s],omega[n_s],sigma[n_s],delta[n_s];
  Double  beta,eps1,eps2;
  Double  rsd_sq;
  Double  tmp_eps;
  int     G5=Ns*Ns-1;

  n_count=0;
  rsd_sq=1.0e-9;
  MaxCG =1000;

  zol_const(qq1,qq,cq);
  
  qq1 *= (2.1/2.5)*(2.1/2.5);
  
  for(int k=1; k< n_s; ++k)
    { 
      qq[k] *= (2.1/2.5)*(2.1/2.5);
      cq[k] *= (2.1/2.5);
    }
  
  for(int k = 1; k < n_s; ++k)
    {
      qq[k]-=qq1;
    }
   
  x[0]=zero;
  r[0]=where(isign==PLUS,psi,Gamma(G5)*psi);

  for(int k = 1; k < n_s; ++k)
    {
      x[k]=zero;
      r[k]=r[0];
    }
  
  
  beta=sqrt(norm2(r[0]));
  v1=Real(1.0/beta)*r[0];
  

  for(int k = 0; k < n_s; ++k)
    {
      p[k]=r[k];
      sigma[k]=beta;
      omega[k]=0.0;
      gamma[k]=1.0;
    }
  
  v0=zero;
  
  for(int l = 1; l <= MaxCG; ++l)
    {
      
      wldirac(v1,r_aux,u,kappa,PLUS);
      wldirac(r_aux,v2,u,kappa,MINUS);
      
      delta[0]=innerProduct(v2,v1);
      delta[0]+=qq1;
      
      for(int k = 1; k < n_s; ++k)
	{
	  delta[k]=delta[0]+DComplex(qq[k]);
	}

      v2-=Complex(delta[0])*v1+Real(beta)*v0;
      beta=sqrt(norm2(v2));
      
      v0=v1;
      v1=Real(1.0/beta)*v2;
      
      for(int k = 0; k < n_s; ++k)
	{
	  gamma[k]=1.0/(delta[k]-omega[k]/gamma[k]);
	  omega[k]=(beta*gamma[k])*(beta*gamma[k]);
	  sigma[k]=-beta*gamma[k]*sigma[k];
	  
	  x[k]+=Complex(gamma[k])*p[k];
	  r[k]=Complex(sigma[k])*v1;
	  p[k]=r[k]+Complex(omega[k])*p[k];
	}
      
      eps1=norm2(r[0]);
      eps2=norm2(x[0]);
      tmp_eps=sqrt(eps1/eps2);
      
//      cout << "iter = " << l << "  " << tmp_eps << endl;
      
      if ( toBool( tmp_eps <=  rsd_sq) )  //why?
	{
	  n_count = l;
	  break;
	}
      
    }

  for(int k = 0; k < n_s; ++k)
    {
      wldirac(x[k],r_aux,u,kappa,PLUS);
      wldirac(r_aux,v1,u,kappa,MINUS);
      v1+=Real(qq1)*x[k];

      v1+=Real(qq[k])*x[k];
      v1-=psi;
      eps1=norm2(v1);
      eps2=norm2(x[k]);
      tmp_eps=sqrt(eps1/eps2);
//      cout << "tmp_eps = " << tmp_eps << endl;
    }
  
  v1=zero;
  for(int k = 0; k < n_s; ++k)
    {
      v1+=Real(cq[k])*x[k];
    }
  
  if (isign == PLUS)
  {
    wldirac(v1,chi,u,kappa,PLUS);
    chi+=psi;
  }
  else
  {
    wldirac(v1,v2,u,kappa,PLUS);
    chi=psi+Gamma(G5)*v2;
  }


  QDPIO::cout << "Inner CG iters = " << n_count << endl;

  END_CODE("OverlapBULinOp");

  return chi;
}


static void wldirac(const LatticeFermion& psi, 
		    LatticeFermion& chi,
		    const multi1d<LatticeColorMatrix>& u, 
		    const Real& Kappa,
		    int   isign)

{
  chi = psi;
  
  
  for(int mu = 0; mu < Nd; ++mu)
  {
    chi -= Kappa*(spinReconstruct(LatticeHalfFermion(u[mu] * shift(spinProject(psi,mu,-isign), FORWARD, mu)),mu,-isign)
		  + spinReconstruct(LatticeHalfFermion(shift(adj(u[mu]) * spinProject(psi,mu,+isign), BACKWARD, mu)),mu,+isign));
  }
}



void zol_const(Double& qq1, Double qq[], Double cq[])
{
  qq1   =0.4203852824188996e-04;
  qq[ 0]=0.0;
  qq[ 1]=0.4230339571009471e-03;
  qq[ 2]=0.1458829357287322e-02;
  qq[ 3]=0.3894118671815802e-02;
  qq[ 4]=0.9481580995276351e-02;
  qq[ 5]=0.2225210303473774e-01;
  qq[ 6]=0.5146873615634705e-01;
  qq[ 7]=0.1185868564259925e+00;
  qq[ 8]=0.27428938359092e+00;
  qq[ 9]=0.6437234073136943e+00;
  qq[10]=0.1567367648339766e+01;
  qq[11]=0.4183844803034055e+01;
  qq[12]=0.1442795672202632e+02;
  qq[13]=0.1451886134043600e+03;
  cq[ 0]=0.8371803911435546e-02;
  cq[ 1]=0.9813674433769438e-02;
  cq[ 2]=0.1294644813006436e-01;
  cq[ 3]=0.1831241561050564e-01;
  cq[ 4]=0.2684553134920763e-01;
  cq[ 5]=0.4004971283185919e-01;
  cq[ 6]=0.6031838397347246e-01;
  cq[ 7]=0.9155798451341547e-01;
  cq[ 8]=0.1406107013842408e+00;
  cq[ 9]=0.2211980404641135e+00;
  cq[10]=0.3673892837229880e+00;
  cq[11]=0.6933239005020095e+00;
  cq[12]=0.1812368256185377e+01;
  cq[13]=0.1555827970041009e+02;
  
}



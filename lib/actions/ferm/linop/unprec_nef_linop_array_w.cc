// $Id: unprec_nef_linop_array_w.cc,v 1.3 2004-09-01 03:32:59 kostas Exp $
/*! \file
 *  \brief Unpreconditioned NEF domain-wall linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/unprec_nef_linop_array_w.h"

//! Creation routine
/*! \ingroup fermact
 *
 * \param u_            gauge field   (Read)
 * \param WilsonMass_   DWF height    (Read)
 * \param b5_           NEF parameter (Read)
 * \param c5_           NEF parameter (Read)
 * \param m_q_          quark mass    (Read)
 */
void UnprecNEFDWLinOpArray::create(const multi1d<LatticeColorMatrix>& u_, 
				   const Real& WilsonMass_, const Real& b5_,
				   const Real& c5_, const Real& m_q_, int N5_)
{
  WilsonMass = WilsonMass_;
  b5 = b5_ ;
  c5 = c5_ ;
  m_q = m_q_ ;
  N5  = N5_ ;

  D.create(u_);
//    CoeffWilsr_s = (AnisoP) ? Wilsr_s / xiF_0 : 1;
}



//! Apply unpreconditioned domain-wall fermion linear operator
/*!
 * \ingroup linop
 *
 * The operator acts on the entire lattice
 *
 * \param psi 	  Pseudofermion field     	       (Read)
 * \param isign   Flag ( PLUS | MINUS )   	       (Read)
 */
void UnprecNEFDWLinOpArray::operator() (multi1d<LatticeFermion>& chi, 
				     const multi1d<LatticeFermion>& psi, 
				     enum PlusMinus isign) const
{
//  multi1d<LatticeFermion> chi(N5);   // probably should check psi.size() == N5

  START_CODE();

  //
  //  Chi   =  D' Psi
  //
  LatticeFermion  tmp,c_tmp,cmb ;
  Real fb5 = b5*(Nd - WilsonMass) + 1;
  Real fc5 = c5*(Nd - WilsonMass) - 1;
  Real Dfb5(-.5*b5);
  Real Dfc5(-.5*c5);

  switch (isign)
  {
  case PLUS:
    for(int n=0; n < N5; ++n){
      if (n == 0)
	//c_tmp=0.5*(psi[1]-
	//   m_q*psi[N5-1] - GammaConst<Ns,Ns*Ns-1>()*(m_q*psi[N5-1]+psi[1])) ;
	c_tmp = chiralProjectMinus(psi[1]) -m_q*chiralProjectPlus(psi[N5-1]);
      else if (n == N5-1)
	//c_tmp = 0.5*(psi[N5-2]  - 
	//   m_q*psi[0] + GammaConst<Ns,Ns*Ns-1>()*(psi[N5-2] + m_q*psi[0])) ;
	c_tmp = chiralProjectPlus(psi[N5-2]) - m_q*chiralProjectMinus(psi[0]) ;
      else
	//c_tmp=0.5*(psi[n+1]+
	//	   psi[n-1]+GammaConst<Ns,Ns*Ns-1>()*(psi[n-1]-psi[n+1]));
	c_tmp = chiralProjectPlus(psi[n-1]) + chiralProjectMinus(psi[n+1]);
      
      cmb = Dfb5*psi[n] + Dfc5*c_tmp ;
      D(tmp, cmb, isign);
      chi[n] = tmp + fb5*psi[n] +fc5*c_tmp ;
    }
    break;

  case MINUS:
    multi1d<LatticeFermion>  Dpsi(N5), c5Dpsi(N5);
    for(int n=0; n < N5; ++n)
      {
	D(Dpsi[n], psi[n], isign);
	c5Dpsi[n] = fc5*psi[n] + Dfc5*Dpsi[n] ;
      }

    for(int n=0; n < N5; ++n){
      if (n == 0){
	c_tmp=chiralProjectPlus(c5Dpsi[1])-m_q*chiralProjectMinus(c5Dpsi[N5-1]);
      }
      else if (n == N5-1)
	c_tmp=chiralProjectMinus(c5Dpsi[N5-2])-m_q*chiralProjectPlus(c5Dpsi[0]);
      else
	c_tmp=chiralProjectMinus(c5Dpsi[n-1]) + chiralProjectPlus(c5Dpsi[n+1]);
      
      chi[n] = fb5*psi[n] +Dfb5*Dpsi[n] + c_tmp ;
    }
    break;
  }

  END_CODE();
}

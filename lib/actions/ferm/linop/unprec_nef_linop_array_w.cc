// $Id: unprec_nef_linop_array_w.cc,v 1.1 2004-08-08 11:12:19 kostas Exp $
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
  LatticeFermion  tmp;
  Real fact1 = b5*(Nd - WilsonMass) + 1;

  multi1d<LatticeFermion> cmb(N5) ;
  CompbineHoppingVectors(cmb, chi)
  switch (isign)
  {
  case PLUS:
    for(int n=0; n < N5; ++n)
    {
      D(tmp, cmb[n], isign);

      if (n == 0)
	chi[n] = tmp + fact1*psi[n] 
	       + m_q*chiralProjectPlus(psi[N5-1]) - chiralProjectMinus(psi[1]);
      else if (n == N5-1)
	chi[n] = tmp + fact1*psi[n] 
	       - chiralProjectPlus(psi[N5-2]) + m_q*chiralProjectMinus(psi[0]);
      else
	chi[n] = tmp + fact1*psi[n] 
	       - chiralProjectPlus(psi[n-1]) - chiralProjectMinus(psi[n+1]);
    }          
    break;

  case MINUS:
    for(int n=0; n < N5; ++n)
    {
      D(tmp, psi[n], isign);

      if (n == 0)
	chi[n] = tmp + fact1*psi[n] 
	       + m_q*chiralProjectMinus(psi[N5-1]) - chiralProjectPlus(psi[1]);
      else if (n == N5-1)
	chi[n] = tmp + fact1*psi[n] 
	       - chiralProjectMinus(psi[N5-2]) + m_q*chiralProjectPlus(psi[0]);
      else
	chi[n] = tmp + fact1*psi[n] 
	       - chiralProjectMinus(psi[n-1]) - chiralProjectPlus(psi[n+1]);
    }          
    break;
  }

  END_CODE();
}

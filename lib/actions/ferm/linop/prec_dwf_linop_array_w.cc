// $Id: prec_dwf_linop_array_w.cc,v 1.1 2003-11-22 21:34:01 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned domain-wall linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/prec_dwf_linop_array_w.h"

//! Creation routine
/*! \ingroup fermact
 *
 * \param u_            gauge field   (Read)
 * \param WilsonMass_   DWF height    (Read)
 * \param m_q_          quark mass    (Read)
 */
void 
EvenOddPrecDWLinOpArray::create(const multi1d<LatticeColorMatrix>& u_, 
				const Real& WilsonMass_, const Real& m_q_, int N5_)
{
  WilsonMass = WilsonMass_;
  m_q = m_q_;
  a5  = 1.0;
  N5  = N5_;

  D.create(u_);
}


//! Apply unpreconditioned domain-wall fermion linear operator
/*!
 * \ingroup linop
 *
 * The operator acts on the entire lattice
 *
 * \param chi 	  result     	                       (Write)
 * \param psi 	  source     	                       (Read)
 * \param isign   Flag ( PLUS | MINUS )   	       (Read)
 */
void 
EvenOddPrecDWLinOpArray::operator() (multi1d<LatticeFermion>& chi, 
				     const multi1d<LatticeFermion>& psi, 
				     enum PlusMinus isign) const
{
//  multi1d<LatticeFermion> chi(N5);   // probably should check psi.size() == N5

  multi1d<LatticeFermion> tmp1(N5);
  multi1d<LatticeFermion> tmp2(N5);

  /*  Tmp1   =  D     A^(-1)     D    Psi  */
  /*      O      O,E        E,E   E,O    O */
  evenOddLinOp(tmp1, psi, isign);
  evenEvenInvLinOp(tmp2, tmp1, isign);
  oddEvenLinOp(tmp1, tmp2, isign);

  /*  Chi   =  A    Psi  -  Tmp1  */
  /*     O      O,O    O        O */
  oddOddLinOp(chi, psi, isign);
  chi[rb[1]] -= tmp1;
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
void 
EvenOddPrecDWLinOpArray::operator() (multi1d<LatticeFermion>& chi, 
				     const multi1d<LatticeFermion>& psi, 
				     enum PlusMinus isign) const
{
//  multi1d<LatticeFermion> chi(N5);   // probably should check psi.size() == N5

  START_CODE("EvenOddPrecDWLinOpArray");

  //
  //  Chi   =  D' Psi
  //
  LatticeFermion  tmp;
  Real fact1 = a5*(Nd - WilsonMass) + 1;
  Real fact2 = -0.5*a5;

  switch (isign)
  {
  case PLUS:
    for(int n=0; n < N5; ++n)
    {
      D(tmp, psi[n], isign);

      if (n == 0)
	chi[n] = fact2*tmp + fact1*psi[n] 
	       + m_q*chiralProjectPlus(psi[N5-1]) - chiralProjectMinus(psi[1]);
      else if (n == N5-1)
	chi[n] = fact2*tmp + fact1*psi[n] 
	       - chiralProjectPlus(psi[N5-2]) + m_q*chiralProjectMinus(psi[0]);
      else
	chi[n] = fact2*tmp + fact1*psi[n] 
	       - chiralProjectPlus(psi[n-1]) - chiralProjectMinus(psi[n+1]);
    }          
    break;

  case MINUS:
    for(int n=0; n < N5; ++n)
    {
      D(tmp, psi[n], isign);

      if (n == 0)
	chi[n] = fact2*tmp + fact1*psi[n] 
	       + m_q*chiralProjectMinus(psi[N5-1]) - chiralProjectPlus(psi[1]);
      else if (n == N5-1)
	chi[n] = fact2*tmp + fact1*psi[n] 
	       - chiralProjectMinus(psi[N5-2]) + m_q*chiralProjectPlus(psi[0]);
      else
	chi[n] = fact2*tmp + fact1*psi[n] 
	       - chiralProjectMinus(psi[n-1]) - chiralProjectPlus(psi[n+1]);
    }          
    break;
  }

  END_CODE("EvenOddPrecDWLinOpArray");
}

// $Id: unprec_dwf_linop_w.cc,v 1.6 2003-11-20 05:43:41 edwards Exp $
/*! \file
 *  \brief Unpreconditioned domain-wall linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/unprec_dwf_linop_w.h"

//! Creation routine
/*! \ingroup fermact
 *
 * \param u_            gauge field   (Read)
 * \param WilsonMass_   DWF height    (Read)
 * \param m_q_          quark mass    (Read)
 */
void UnprecDWLinOp::create(const multi1d<LatticeColorMatrix>& u_, 
			   const Real& WilsonMass_, const Real& m_q_)
{
  u = u_;
  WilsonMass = WilsonMass_;
  m_q = m_q_;
  a5  = 1.0;

//    CoeffWilsr_s = (AnisoP) ? Wilsr_s / xiF_0 : 1;
}


//-----------------------------------------------------------------------------
// HACK
/* Terrible implementation just to get the silly thing going */
static inline LatticeDWFermion 
cycleUpDW(const LatticeDWFermion& l)
{
  LatticeDWFermion d = zero;

  for(int i=0; i < Ls-1; ++i)
    pokeDW(d, peekDW(l,i+1), i);

  return d;
}

// HACK
/* Terrible implementation just to get the silly thing going */
static inline LatticeDWFermion 
cycleDownDW(const LatticeDWFermion& l)
{
  LatticeDWFermion d = zero;

  for(int i=1; i < Ls; ++i)
    pokeDW(d, peekDW(l,i-1), i);

  return d;
}


// HACK
/* Terrible implementation just to get the silly thing going */
static inline LatticeDWFermion 
chiralProjectPlus(const LatticeDWFermion& l)
{
  return 0.5*(l + Gamma(15)*l);
}

static inline LatticeDWFermion 
chiralProjectMinus(const LatticeDWFermion& l)
{
  return 0.5*(l - Gamma(15)*l);
}


static inline LatticeFermion 
chiralProjectPlus(const LatticeFermion& l)
{
  return 0.5*(l + Gamma(15)*l);
}

static inline LatticeFermion 
chiralProjectMinus(const LatticeFermion& l)
{
  return 0.5*(l - Gamma(15)*l);
}


//-----------------------------------------------------------------------------



//! Apply unpreconditioned domain-wall fermion linear operator
/*!
 * \ingroup linop
 *
 * The operator acts on the entire lattice
 *
 * \param psi 	  Pseudofermion field     	       (Read)
 * \param isign   Flag ( PLUS | MINUS )   	       (Read)
 */
void UnprecDWLinOp::operator() (LatticeDWFermion& chi, const LatticeDWFermion& psi, 
				enum PlusMinus isign) const
{
  START_CODE("UnprecDWLinOp");

  if (Nd != 4)
    QDP_error_exit("UnprecDWLinOP: expects 4D");

  //
  //  Chi   =  D' Psi
  //
  LatticeDWFermion tmp1, tmp2;

  Real fact1 = a5*(Nd - WilsonMass) + 1;
  Real fact2 = -0.5*a5;

  /* Why are these lines split? An array syntax would help, but the problem is deeper.
   * The expression templates require NO variable args (like int's) to a function
   * and all args must be known at compile time. Hence, the function names carry
   * (as functions usually do) the meaning (and implicit args) to a function.
   */
  switch (isign)
  {
  case PLUS:
    tmp1 = chiralProjectMinus(cycleUpDW(psi));
    pokeDW(tmp1, -m_q*chiralProjectMinus(LatticeFermion(peekDW(psi,0))), Ls-1);
    tmp2 = chiralProjectPlus(cycleDownDW(psi));
    pokeDW(tmp2, -m_q*chiralProjectPlus(LatticeFermion(peekDW(psi,Ls-1))), 0);

    chi = fact1*psi - tmp1 - tmp2
        + fact2*(spinReconstructDir0Minus(u[0] * shift(spinProjectDir0Minus(psi), FORWARD, 0)) +
		 spinReconstructDir0Plus(shift(adj(u[0]) * spinProjectDir0Plus(psi), BACKWARD, 0)) +
		 spinReconstructDir1Minus(u[1] * shift(spinProjectDir1Minus(psi), FORWARD, 1)) +
		 spinReconstructDir1Plus(shift(adj(u[1]) * spinProjectDir1Plus(psi), BACKWARD, 1)) +
		 spinReconstructDir2Minus(u[2] * shift(spinProjectDir2Minus(psi), FORWARD, 2)) +
		 spinReconstructDir2Plus(shift(adj(u[2]) * spinProjectDir2Plus(psi), BACKWARD, 2)) +
		 spinReconstructDir3Minus(u[3] * shift(spinProjectDir3Minus(psi), FORWARD, 3)) +
		 spinReconstructDir3Plus(shift(adj(u[3]) * spinProjectDir3Plus(psi), BACKWARD, 3)));
          
    break;

  case MINUS:
    tmp1 = chiralProjectPlus(cycleUpDW(psi));
    pokeDW(tmp1, -m_q*chiralProjectPlus(LatticeFermion(peekDW(psi,0))), Ls-1);
    tmp2 = chiralProjectMinus(cycleDownDW(psi));
    pokeDW(tmp2, -m_q*chiralProjectMinus(LatticeFermion(peekDW(psi,Ls-1))), 0);

    chi = fact1*psi - tmp1 - tmp2
        + fact2*(spinReconstructDir0Plus(u[0] * shift(spinProjectDir0Plus(psi), FORWARD, 0)) +
		 spinReconstructDir0Minus(shift(adj(u[0]) * spinProjectDir0Minus(psi), BACKWARD, 0)) +
		 spinReconstructDir1Plus(u[1] * shift(spinProjectDir1Plus(psi), FORWARD, 1)) +
		 spinReconstructDir1Minus(shift(adj(u[1]) * spinProjectDir1Minus(psi), BACKWARD, 1)) +
		 spinReconstructDir2Plus(u[2] * shift(spinProjectDir2Plus(psi), FORWARD, 2)) +
		 spinReconstructDir2Minus(shift(adj(u[2]) * spinProjectDir2Minus(psi), BACKWARD, 2)) +
		 spinReconstructDir3Plus(u[3] * shift(spinProjectDir3Plus(psi), FORWARD, 3)) +
		 spinReconstructDir3Minus(shift(adj(u[3]) * spinProjectDir3Minus(psi), BACKWARD, 3)));
    break;
  }


  END_CODE("UnprecDWLinOp");
}

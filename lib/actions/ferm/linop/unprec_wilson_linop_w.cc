// $Id: unprec_wilson_linop_w.cc,v 1.4 2003-11-16 06:21:49 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Wilson linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/unprec_wilson_linop_w.h"

//! Creation routine
/*! \ingroup fermact
 *
 * \param u_ 	    gauge field     	       (Read)
 * \param Mass_   fermion kappa   	       (Read)
 */
void UnprecWilsonLinOp::create(const multi1d<LatticeColorMatrix>& u_, const Real& Mass_)
{
  Mass = Mass_;
  u = u_;
  D.create(u);

//    CoeffWilsr_s = (AnisoP) ? Wilsr_s / xiF_0 : 1;
}


//! Apply unpreconditioned Wilson fermion linear operator
/*!
 * \ingroup linop
 *
 * The operator acts on the entire lattice
 *
 * \param psi 	  Pseudofermion field     	       (Read)
 * \param isign   Flag ( PLUS | MINUS )   	       (Read)
 */
LatticeFermion UnprecWilsonLinOp::operator() (const LatticeFermion& psi, enum PlusMinus isign) const
{
  LatticeFermion chi;

  START_CODE("UnprecWilsonLinOp");

  //
  //  Chi   =  (Nd+Mass)*Psi  -  (1/2) * D' Psi
  //
  Real fact1 = Nd + Mass;
  Real fact2 = -0.5;

  chi = fact1*psi + fact2* D(psi, isign);
  
  END_CODE("UnprecWilsonLinOp");

  return chi;
}

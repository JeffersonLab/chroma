// $Id: unprec_wilson_linop_w.cc,v 1.3 2003-11-09 22:35:19 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Wilson linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/unprec_wilson_linop_w.h"

//! Creation routine
/*! \ingroup fermact
 *
 * \param _u 	    gauge field     	       (Read)
 * \param _Kappa   fermion kappa   	       (Read)
 */
void UnprecWilsonLinOp::create(const multi1d<LatticeColorMatrix>& _u, const Real& _Kappa)
{
  Kappa = _Kappa;
  u = _u;
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
  //  Chi   =  Psi  -  k   D' Psi
  //
  chi = psi - Kappa * D(psi, isign);
  
  END_CODE("UnprecWilsonLinOp");

  return chi;
}

// $Id: lmpsim_w.cc,v 1.2 2003-04-09 19:43:22 edwards Exp $
/*! \file
 *  \brief Preconditioned Wilson linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/lmpsim_w.h"

//! Creation routine
/*!
 * \ingroup linop
 *
 * \param _u 	    gauge field     	       (Read)
 * \param _Kappa   fermion kappa   	       (Read)
 */
void PreconditionedWilson::create(const multi1d<LatticeColorMatrix>& _u, const Real& _Kappa)
{
  Kappa = _Kappa;
  u = _u;
  D.create(u);

//    CoeffWilsr_s = (AnisoP) ? Wilsr_s / xiF_0 : 1;
}

//! Apply fermion linear operator
/*!
 * \ingroup linop
 *
 * \param psi 	  Pseudofermion field     	       (Read)
 * \param isign   Flag ( PLUS | MINUS )   	       (Read)
 */
LatticeFermion PreconditionedWilson::operator() (const LatticeFermion& psi, enum LinOpSign isign) const
{
  LatticeFermion chi;

  START_CODE("lmpsi");

  Real Kappa_sq = Kappa * Kappa;

  /*      	         2                 */
  /*  Chi   =  Psi  -  k   D'    D'    Psi */
  /*     O        O         O,E   E,O      */
  chi[rb[1]] = psi - Kappa_sq * D.apply(D.apply(psi, isign, 0), isign, 1);
  
  END_CODE("lmpsi");

  return chi;
}

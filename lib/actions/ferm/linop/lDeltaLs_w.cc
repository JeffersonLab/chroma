// $Id: lDeltaLs_w.cc,v 1.1 2004-10-15 16:37:19 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Wilson linear operator
 */

#include "chromabase.h"
#include "lDeltaLs_w.h"

using namespace QDP;
using namespace Chroma;

namespace Chroma { 
//! Apply unpreconditioned Wilson fermion linear operator
/*!
 * \ingroup linop
 *
 * The operator acts on the entire lattice
 *
 * \param chi 	  Pseudofermion field     	       (Read)
 * \param psi 	  Pseudofermion field     	       (Read)
 * \param isign   Flag ( PLUS | MINUS )   	       (Read)
 */
void lDeltaLs::operator() (LatticeFermion& chi, const LatticeFermion& psi, 
			   enum PlusMinus isign) const
{
  START_CODE();

  int G5 = Ns*Ns-1;

  LatticeFermion tmp1, tmp2;

  // Apply g5 epsilon psi
  (*g5eps)(tmp1, psi, isign);

  // Multiply by g5 -> tmp2 = epsilon psi
  tmp2 = Gamma(G5)*tmp1;

  // Multiply by g5eps -> tmp1 = gamma_5 epsilon tmp2 
  //                           = gamma_5 epsilon epsilon psi
  (*g5eps)(tmp1, tmp2, isign);

  // Multiply by gamma5 -> tmp2 = gamma_5*tmp1 = epsilon^2 psi
  tmp2 = Gamma(G5)*tmp1;

  // Chi = psi-tmp2 = psi - epsilon^2 psi
  //     =           ( 1 - epsilon^2) psi
  chi = psi - tmp2;

  // Multiply in factor:
  Real factor = Real(0.25);

  // chi *= (1/4) = (1/4)(1-epsilon^2) psi = DeltaLs psi
  chi *=factor;

  END_CODE();
}

};

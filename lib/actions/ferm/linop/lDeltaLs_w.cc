// $Id: lDeltaLs_w.cc,v 1.2 2004-12-09 03:58:03 edwards Exp $
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
  LatticeFermion tmp1, tmp2, tmp3;

  // Construct  eps(H)*psi
  (*D)(tmp1, psi, PLUS);
  tmp2 = Gamma(G5)*(2*tmp1 - psi);

  // Construct  eps(H)*eps(H)*psi
  (*D)(tmp1, tmp2, PLUS);
  tmp3 = Gamma(G5)*(2*tmp1 - tmp2);

  // Construct  (1/4)(1 - eps(H)*eps(H))*psi
  chi = 0.25*(psi - tmp3);

  END_CODE();
}

};

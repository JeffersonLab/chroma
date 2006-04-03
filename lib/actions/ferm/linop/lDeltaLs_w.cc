// $Id: lDeltaLs_w.cc,v 3.0 2006-04-03 04:58:50 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Wilson linear operator
 */

#include "chromabase.h"
#include "lDeltaLs_w.h"

using namespace QDP::Hints;

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
  moveToFastMemoryHint(tmp1);
  moveToFastMemoryHint(tmp2);
  moveToFastMemoryHint(tmp3);

  // Construct  eps(H)*psi
  (*D)(tmp1, psi, PLUS);
  tmp2 = GammaConst<Ns,Ns*Ns-1>()*(Real(2)*tmp1 - Real(1)*psi);

  // Construct  eps(H)*eps(H)*psi
  (*D)(tmp1, tmp2, PLUS);
  tmp3 = GammaConst<Ns,Ns*Ns-1>()*(Real(2)*tmp1 - Real(1)*tmp2);

  // Construct  (1/4)(1 - eps(H)*eps(H))*psi
  chi = Real(0.25)*(psi - tmp3);

  END_CODE();
}

};

// $Id: asqtad_linop_s.cc,v 1.6 2004-11-06 11:30:45 mcneile Exp $
/*! \file
 *  \brief Unpreconditioned Wilson linear operator
 */
// NEW $Id: asqtad_linop.cc 2003/11/13 steve

#include "chromabase.h"
#include "linop.h"


void AsqtadLinOp::evenOddLinOp(LatticeStaggeredFermion& chi, 
					  const LatticeStaggeredFermion& psi, 
					  enum PlusMinus isign) const
{
  D.apply(chi, psi, isign, 0);
}

void AsqtadLinOp::oddEvenLinOp(LatticeStaggeredFermion& chi, 
					  const LatticeStaggeredFermion& psi, 
					  enum PlusMinus isign) const
{
  D.apply(chi, psi, isign, 1);
}

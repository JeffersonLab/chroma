// $Id: asqtad_linop_s.cc,v 1.5 2004-01-07 13:50:07 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Wilson linear operator
 */
// NEW $Id: asqtad_linop.cc 2003/11/13 steve

#include "chromabase.h"
#include "linop.h"


void AsqtadLinOp::evenOddLinOp(LatticeFermion& chi, 
					  const LatticeFermion& psi, 
					  enum PlusMinus isign) const
{
  D.apply(chi, psi, isign, 0);
}

void AsqtadLinOp::oddEvenLinOp(LatticeFermion& chi, 
					  const LatticeFermion& psi, 
					  enum PlusMinus isign) const
{
  D.apply(chi, psi, isign, 1);
}

// $Id: asqtad_linop_s.cc,v 1.4 2003-12-12 14:28:27 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Wilson linear operator
 */
// NEW $Id: asqtad_linop.cc 2003/11/13 steve

#include "chromabase.h"
#include "actions/ferm/linop/asqtad_linop_s.h"


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

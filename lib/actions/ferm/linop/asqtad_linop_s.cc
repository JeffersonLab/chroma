// $Id: asqtad_linop_s.cc,v 1.3 2003-12-12 13:56:40 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Wilson linear operator
 */
// NEW $Id: asqtad_linop.cc 2003/11/13 steve

#include "chromabase.h"
#include "actions/ferm/linop/asqtad_linop_s.h"
#include "actions/ferm/linop/dslash_s.h"


void AsqtadLinOp::evenOddLinOp(LatticeFermion& chi, 
					  const LatticeFermion& psi, 
					  enum PlusMinus isign) const
{
  AsqtadDslash D(u_fat, u_triple);

  LatticeFermion tmp;

  D.apply(chi, psi, isign, 0);

}

void AsqtadLinOp::oddEvenLinOp(LatticeFermion& chi, 
					  const LatticeFermion& psi, 
					  enum PlusMinus isign) const
{
  AsqtadDslash D(u_fat, u_triple);

  LatticeFermion tmp;

  D.apply(chi, psi, isign, 1);

}

// $Id: prec_asqtad_linop_s.cc,v 1.3 2003-12-11 17:11:17 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Wilson linear operator
 */
// NEW $Id: asqtad_linop.cc 2003/11/13 steve

#include "chromabase.h"
#include "actions/ferm/linop/prec_asqtad_linop_s.h"
#include "actions/ferm/linop/dslash_s.h"


void EvenOddPrecAsqtadLinOp::evenEvenLinOp(LatticeFermion& chi, 
					   const LatticeFermion& psi, 
					   enum PlusMinus isign) const 
{
  chi[rb[0]] = 2*Mass*psi;
}


void EvenOddPrecAsqtadLinOp::oddOddLinOp(LatticeFermion& chi, 
					 const LatticeFermion& psi, 
					 enum PlusMinus isign) const
{
  chi[rb[1]] = 2*Mass*psi;
}

void EvenOddPrecAsqtadLinOp::evenOddLinOp(LatticeFermion& chi, 
					  const LatticeFermion& psi, 
					  enum PlusMinus isign) const
{
  AsqtadDslash D(u_fat, u_triple);

  LatticeFermion tmp;

  D.apply(chi, psi, isign, 0);

}

void EvenOddPrecAsqtadLinOp::oddEvenLinOp(LatticeFermion& chi, 
					  const LatticeFermion& psi, 
					  enum PlusMinus isign) const
{
  AsqtadDslash D(u_fat, u_triple);

  LatticeFermion tmp;

  D.apply(chi, psi, isign, 1);

}

// $Id: prec_asqtad_linop_s.cc,v 1.2 2003-12-10 16:21:00 bjoo Exp $
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
  chi[rb[0]] = zero;
}


void EvenOddPrecAsqtadLinOp::oddOddLinOp(LatticeFermion& chi, 
					 const LatticeFermion& psi, 
					 enum PlusMinus isign) const
{
  chi[rb[1]] = zero;
}

void EvenOddPrecAsqtadLinOp::evenOddLinOp(LatticeFermion& chi, 
					  const LatticeFermion& psi, 
					  enum PlusMinus isign) const
{
  AsqtadDslash D(u_fat, u_triple);

  LatticeFermion tmp;

  D.apply(tmp, psi, isign, 0);

  chi[rb[0]]  = 2*Mass*psi;
  chi[rb[0]] += tmp;
}

void EvenOddPrecAsqtadLinOp::oddEvenLinOp(LatticeFermion& chi, 
					  const LatticeFermion& psi, 
					  enum PlusMinus isign) const
{
  AsqtadDslash D(u_fat, u_triple);

  LatticeFermion tmp;

  PlusMinus isignNew = (isign == PLUS) ? MINUS : PLUS;
  D.apply(tmp, psi, isignNew, 1);

  chi[rb[1]] = 2*Mass*psi;
  chi[rb[1]] += tmp;
}

// $Id: asqtad_linop_s.cc,v 2.0 2005-09-25 21:04:28 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Asqtad linear operator
 */
// NEW $Id: asqtad_linop.cc 2003/11/13 steve

#include "chromabase.h"
#include "actions/ferm/linop/asqtad_linop_s.h"

namespace Chroma 
{ 
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

}; // End Namespace Chroma


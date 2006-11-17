// $Id: asqtad_linop_s.cc,v 3.1 2006-11-17 02:54:47 edwards Exp $
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

  //! Return flops performed by the operator()
  unsigned long AsqtadLinOp::nFlops() const
  {
    unsigned long cbsite_flops = 0;
    return cbsite_flops*(Layout::sitesOnNode()/2);
  }

} // End Namespace Chroma


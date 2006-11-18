// $Id: asqtad_linop_s.cc,v 3.2 2006-11-18 02:33:03 kostas Exp $
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
    unsigned long cbsite_flops = 1146; // I think this is correct... see flop count in asq_dsl_s.cc
    return cbsite_flops*(Layout::sitesOnNode()/2);
  }

} // End Namespace Chroma


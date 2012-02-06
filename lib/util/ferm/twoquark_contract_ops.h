// -*- C++ -*-
/*! \file
 * \brief Contraction operators for two quarks
 */

#ifndef __twoquark_contract_ops_h__
#define __twoquark_contract_ops_h__

#include "util/ferm/spin_rep.h"
#include <vector>

namespace Chroma
{
  //! Take transpose of a distilled object
  void transpose(multi2d<LatticeColorVector>& dist_rep, 
		 const multi2d<LatticeColorVector>& prop_rep);

  //! Use gamma_5 hermiticity on a prop
  void gamma5Herm(multi2d<LatticeColorVector>& prop);

  //----------------------------------------------------------------------------
  //! Dist(t2) = SpinMatrix*Prop(t2)
  void multiplyRep(multi2d<LatticeColorVector>& dist_rep, 
		   const std::vector<MatrixSpinRep_t>& spin, const multi2d<LatticeColorVector>& prop_rep);

  //! Dist(t2) = Prop(t2)*SpinMatrix
  void multiplyRep(multi2d<LatticeColorVector>& dist_rep, 
		   const multi2d<LatticeColorVector>& prop_rep, const std::vector<MatrixSpinRep_t>& spin);

} // namespace Chroma

#endif

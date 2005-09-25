// $Id: staggered_qprop.cc,v 2.0 2005-09-25 21:04:31 edwards Exp $
/*! \file
 *  \brief Propagator solver for an even-odd non-preconditioned fermion operator
 *
 *  Solve for the propagator of an even-odd non-preconditioned fermion operator
 */

#include "fermact.h"
#include "actions/ferm/qprop/prec_staggered_qprop.h"

namespace Chroma
{

  typedef LatticeStaggeredFermion LF;
  typedef multi1d<LatticeColorMatrix> LCM;

  /*
  template<>
  const SystemSolver<LF>* 
  EvenOddStaggeredTypeFermAct<LF,LCM>::qprop(Handle<const ConnectState> state,
					     const InvertParam_t& invParam) const
  {
    return new EvenOddFermActQprop<LF,LCM>(Handle< const EvenOddLinearOperator<LF,LCM> >(linOp(state)),
					   Handle< const LinearOperator<LF> >(lMdagM(state)),
					   getQuarkMass(),
					   invParam);
  }
  */

  template<>
  const SystemSolver<LF>* 
  EvenOddStaggeredTypeFermAct<LF,LCM>::qprop(Handle<const ConnectState> state,
					     const InvertParam_t& invParam) const
  {
    return new EvenOddFermActQprop<LF,LCM>(*this,
					   state,
					   invParam);
  }


}  

// $Id: staggered_qprop.cc,v 3.0 2006-04-03 04:58:53 edwards Exp $
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
  SystemSolver<LF>* 
  EvenOddStaggeredTypeFermAct<LF,LCM,LCM>::qprop(Handle< FermState<LF,LCM,LCM> > state,
					     const InvertParam_t& invParam) const
  {
    return new EvenOddFermActQprop<LF,LCM>(Handle< const EvenOddLinearOperator<LF,LCM> >(linOp(state)),
					   Handle< const LinearOperator<LF> >(lMdagM(state)),
					   getQuarkMass(),
					   invParam);
  }
  */

  template<>
  SystemSolver<LF>* 
  EvenOddStaggeredTypeFermAct<LF,LCM,LCM>::qprop(Handle< FermState<LF,LCM,LCM> > state,
						 const InvertParam_t& invParam) const
  {
    return new EvenOddFermActQprop<LF,LCM,LCM>(*this,
					       state,
					       invParam);
  }


}  

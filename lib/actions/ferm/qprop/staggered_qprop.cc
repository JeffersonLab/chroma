// $Id: staggered_qprop.cc,v 3.2 2006-10-19 17:36:07 edwards Exp $
/*! \file
 *  \brief Propagator solver for an even-odd non-preconditioned fermion operator
 *
 *  Solve for the propagator of an even-odd non-preconditioned fermion operator
 */

#include "fermact.h"
#include "actions/ferm/qprop/eoprec_staggered_qprop.h"
#include "actions/ferm/invert/syssolver_cg_params.h"

namespace Chroma
{

  typedef LatticeStaggeredFermion LF;
  typedef multi1d<LatticeColorMatrix> LCM;

  /*
  template<>
  SystemSolver<LF>* 
  EvenOddStaggeredTypeFermAct<LF,LCM,LCM>::qprop(Handle< FermState<LF,LCM,LCM> > state,
					     const GroupXML_t& invParam) const
  {
    return new EvenOddFermActQprop<LF,LCM>(Handle< const EvenOddLinearOperator<LF,LCM> >(linOp(state)),
					   Handle< const LinearOperator<LF> >(lMdagM(state)),
					   getQuarkMass(),
					   Handle< MdagMSystemSolverArray<LF> >(invMdagM(state,invParam)));
  }
  */


  //! Return a linear operator solver for this action to solve MdagM*psi=chi 
  /*! \ingroup qprop */
  template<>
  SystemSolver<LF>* 
  EvenOddStaggeredTypeFermAct<LF,LCM,LCM>::qprop(Handle< FermState<LF,LCM,LCM> > state,
						 const GroupXML_t& invParam) const
  {
    std::istringstream  is(invParam.xml);
    XMLReader  paramtop(is);
    SysSolverCGParams params(paramtop,invParam.path);

    return new EvenOddFermActQprop<LF,LCM,LCM>(*this,
					       state,
					       params);
  }

}  

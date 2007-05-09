// $Id: hisq_qprop.cc,v 1.1 2007-05-09 12:45:31 mcneile Exp $
/*! \file
 *  \brief Propagator solver for an even-odd non-preconditioned fermion operator
 *
 *  Solve for the propagator of an even-odd non-preconditioned fermion operator
 */

#include "chromabase.h"

#include "fermact.h"
#include "actions/ferm/invert/invcg1.h"

#include "actions/ferm/fermacts/hisq_fermact_s.h"
#include "actions/ferm/qprop/hisq_qprop.h"
#include "actions/ferm/invert/syssolver_cg_params.h"

namespace Chroma
{

  typedef LatticeStaggeredFermion LF;
  typedef multi1d<LatticeColorMatrix> LCM;

  SystemSolver<LF>* 
  HisqFermAct::qprop(Handle< FermState<LF,LCM,LCM> > state,
		       const GroupXML_t& invParam) const
  {
    const EvenOddStaggeredTypeFermAct<LF,LCM,LCM>& S_cast = 
      dynamic_cast< const EvenOddStaggeredTypeFermAct<LF,LCM,LCM>& >(*this);

    std::istringstream  is(invParam.xml);
    XMLReader  paramtop(is);
    SysSolverCGParams params(paramtop, invParam.path);

    return new HisqQprop(S_cast,
			   state,
			   params);
  }
  

}  

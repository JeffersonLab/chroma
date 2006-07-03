// $Id: asqtad_qprop.cc,v 3.1 2006-07-03 15:26:09 edwards Exp $
/*! \file
 *  \brief Propagator solver for an even-odd non-preconditioned fermion operator
 *
 *  Solve for the propagator of an even-odd non-preconditioned fermion operator
 */

#include "chromabase.h"

#include "fermact.h"
#include "actions/ferm/invert/invcg1.h"

#include "actions/ferm/fermacts/asqtad_fermact_s.h"
#include "actions/ferm/qprop/asqtad_qprop.h"
#include "actions/ferm/invert/syssolver_cg_params.h"

namespace Chroma
{

  typedef LatticeStaggeredFermion LF;
  typedef multi1d<LatticeColorMatrix> LCM;

  SystemSolver<LF>* 
  AsqtadFermAct::qprop(Handle< FermState<LF,LCM,LCM> > state,
		       const GroupXML_t& invParam) const
  {
    const EvenOddStaggeredTypeFermAct<LF,LCM,LCM>& S_cast = 
      dynamic_cast< const EvenOddStaggeredTypeFermAct<LF,LCM,LCM>& >(*this);

    std::istringstream  is(invParam.xml);
    XMLReader  paramtop(is);
    SysSolverCGParams params(paramtop, invParam.path);

    return new AsqtadQprop(S_cast,
			   state,
			   params);
  }
  

}  

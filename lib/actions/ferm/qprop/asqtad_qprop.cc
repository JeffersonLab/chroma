// $Id: asqtad_qprop.cc,v 2.0 2005-09-25 21:04:30 edwards Exp $
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
namespace Chroma
{

  typedef LatticeStaggeredFermion LF;
  typedef multi1d<LatticeColorMatrix> LCM;

  const SystemSolver<LF>* 
  AsqtadFermAct::qprop(Handle<const ConnectState> state,
		       const InvertParam_t& invParam) const
  {
    
    const EvenOddStaggeredTypeFermAct<LF, LCM>& S_cast = 
      dynamic_cast< const EvenOddStaggeredTypeFermAct<LF,LCM>& >(*this);

    return new AsqtadQprop(S_cast,
			   state,
			   invParam);
  }
  

}  

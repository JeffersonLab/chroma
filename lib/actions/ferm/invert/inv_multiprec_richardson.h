#ifndef _inv_multiprec_richardson_h_
#define _inv_multiprec_richardson_h_

#include <chromabase.h>
#include <syssolver.h>
#include <linearop.h>

namespace Chroma {

  
  void InvMultiPrecRichardson( const SystemSolver< LatticeFermionF >& Dinv,
			       const LinearOperator< LatticeFermionD >& D,
			       const LatticeFermionD& b, 
			       LatticeFermionD& x,
			       int MaxIter,
			       Real RsdTarget,
			       SystemSolverResults_t& res);


}


#endif

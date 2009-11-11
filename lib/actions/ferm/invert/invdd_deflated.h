// -*- C++ -*-
// $Id: invdd_deflated.h,v 1.1 2009/11/11 16:39:50 jbulava Exp $
/*! \file
 *  Domain-Decomposed Deflated inverter  
 */

#ifndef __invdd_deflated__
#define __invdd_deflated__

#include "linearop.h"
#include "syssolver.h"

namespace Chroma {


	namespace InvDDDeflatedEnv {
	
SystemSolverResults_t 
dd_def_inv(const LinearOperator<T>& A,
       const T& chi,
       T& psi,
       const Real& RsdCG, 
       int MaxCG, int MinCG=0);

}  // end namespace Chroma

#endif

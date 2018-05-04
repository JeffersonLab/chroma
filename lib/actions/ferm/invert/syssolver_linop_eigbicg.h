// -*- C++ -*-
/*! \file
 *  \brief Solve a M*psi=chi linear system by EigBiCG
 */

#ifndef __syssolver_eigbicg_h__
#define __syssolver_eigbicg_h__

#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "named_obj.h"
#include "meas/inline/io/named_objmap.h"

#include "actions/ferm/invert/syssolver_linop.h"
#include "actions/ferm/invert/syssolver_OPTeigbicg_params.h"
#include "actions/ferm/invert/containers.h"

#include "util/info/unique_id.h"

#ifdef BUILD_OPT_EIGCG

#include "actions/ferm/invert/syssolver_linop_OPTeigbicg.h"

namespace Chroma 
{
  //! Eigenstd::vector accelerated CG system solver namespace
  namespace LinOpSysSolverEigBiCGEnv
  {
    //! Register the syssolver
    inline bool registerAll() {return LinOpSysSolverOptEigBiCGEnv::registerAll();}
  }
}  // end namespace Chroma

#else

// Generic QDP code does not exist
// do nothing... and the thing will fail when somebody calls eigBiCG when it is not compiled in....
namespace Chroma 
{
  //! Eigenstd::vector accelerated CG system solver namespace
  namespace LinOpSysSolverEigBiCGEnv
  {
    //! Register the syssolver
    inline bool registerAll() {return true;}
  }
}  // end namespace Chroma

#endif 

#endif

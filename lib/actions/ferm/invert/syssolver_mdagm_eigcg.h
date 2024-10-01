// -*- C++ -*-
/*! \file
 *  \brief Solve a M^dag*M*psi=chi linear system by EigCG
 */

#ifndef __syssolver_mdagm_eigcg_h__
#define __syssolver_mdagm_eigcg_h__

#ifdef BUILD_OPT_EIGCG

#include "actions/ferm/invert/syssolver_mdagm_OPTeigcg.h"

namespace Chroma 
{
  //! Eigenstd::vector accelerated CG system solver namespace
  namespace MdagMSysSolverEigCGEnv
  {
    //! Register the syssolver
    inline bool registerAll() {
#if ! defined (QDP_IS_QDPJIT2)
      return MdagMSysSolverOptEigCGEnv::registerAll();
#endif
    }
  }
}  // end namespace Chroma

#else

// Generic QDP code

#include "actions/ferm/invert/syssolver_mdagm_eigcg_qdp.h"

namespace Chroma 
{
  //! Eigenstd::vector accelerated CG system solver namespace
  namespace MdagMSysSolverEigCGEnv
  {
    //! Register the syssolver
    inline bool registerAll() {
#if ! defined (QDP_IS_QDPJIT2)
      return MdagMSysSolverQDPEigCGEnv::registerAll();
#endif
    }
  }
}  // end namespace Chroma

#endif 

#endif



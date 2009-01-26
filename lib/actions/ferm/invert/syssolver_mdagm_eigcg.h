// -*- C++ -*-
// $Id: syssolver_mdagm_eigcg.h,v 3.1 2009-01-26 22:47:06 edwards Exp $
/*! \file
 *  \brief Solve a M^dag*M*psi=chi linear system by EigCG
 */

#ifndef __syssolver_mdagm_eigcg_h__
#define __syssolver_mdagm_eigcg_h__

#ifdef BUILD_OPT_EIGCG

#include "actions/ferm/invert/syssolver_mdagm_OPTeigcg.h"

namespace Chroma 
{
  //! Eigenvector accelerated CG system solver namespace
  namespace MdagMSysSolverEigCGEnv
  {
    //! Register the syssolver
    inline bool registerAll() {return MdagMSysSolverOptEigCGEnv::registerAll();}
  }
}  // end namespace Chroma

#else

// Generic QDP code

#include "actions/ferm/invert/syssolver_mdagm_eigcg_qdp.h"

namespace Chroma 
{
  //! Eigenvector accelerated CG system solver namespace
  namespace MdagMSysSolverEigCGEnv
  {
    //! Register the syssolver
    inline bool registerAll() {return MdagMSysSolverQDPEigCGEnv::registerAll();}
  }
}  // end namespace Chroma

#endif 

#endif


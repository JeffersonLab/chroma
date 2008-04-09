// -*- C++ -*-
// $Id: syssolver_linop_eigcg.h,v 1.12 2008-04-09 04:49:23 kostas Exp $
/*! \file
 *  \brief Solve a M*psi=chi linear system by CG2
 */

#ifndef __syssolver_linop_eigcg_h__
#define __syssolver_linop_eigcg_h__

#ifdef BUILD_OPT_EIGCG

#include "actions/ferm/invert/syssolver_linop_OPTeigcg.h"

namespace Chroma 
{
  //! Eigenvector accelerated CG system solver namespace
  namespace LinOpSysSolverEigCGEnv
  {
    //! Register the syssolver
    inline bool registerAll() {return LinOpSysSolverOptEigCGEnv::registerAll();}
  }
}  // end namespace Chroma

#else

// Generic QDP code

#include "actions/ferm/invert/syssolver_linop_eigcg_qdp.h"

namespace Chroma 
{
  //! Eigenvector accelerated CG system solver namespace
  namespace LinOpSysSolverEigCGEnv
  {
    //! Register the syssolver
    inline bool registerAll() {return LinOpSysSolverQDPEigCGEnv::registerAll();}
  }
}  // end namespace Chroma

#endif 

#endif


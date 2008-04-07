// -*- C++ -*-
// $Id: syssolver_linop_eigcg.h,v 1.11 2008-04-07 04:58:51 edwards Exp $
/*! \file
 *  \brief Solve a M*psi=chi linear system by CG2
 */

#ifndef __syssolver_linop_eigcg_h__
#define __syssolver_linop_eigcg_h__

#ifdef BUILD_LAPACK

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


// -*- C++ -*-
// $Id: syssolver_eigcg_params.h,v 1.3 2007-09-27 14:47:53 kostas Exp $
/*! \file
 *  \brief Solve a CG1 system
 */

#ifndef __syssolver_eigcg_params_h__
#define __syssolver_eigcg_params_h__

#include "chromabase.h"


namespace Chroma
{

  //! Params for EigCG inverter
  /*! \ingroup invert */
  struct SysSolverEigCGParams
  {
    SysSolverEigCGParams();
    SysSolverEigCGParams(XMLReader& in, const std::string& path);
    
    Real  RsdCG ;           /*!< CG residual */
    int   MaxCG ;           /*!< Maximum CG iterations */
    
    int   Neig     ; /*!< number of eigenvectors to compute  */
    int   Nmax     ; /*!< number of stored vectors */
    int   Neig_max ; /*!< maximum number of eigenvectors to be refined */

    Real  RsdCGRestart ; /*!< CG residual for restart. */
    
    int   vPrecCGvecs  ; /*!< number of vectors for preconditioned CG (if <=0 do regular CG) */
  std:string eigen_id ; /*!< named buffer holding the eigenvectors */

  };


  // Reader/writers
  /*! \ingroup invert */
  void read(XMLReader& xml, const string& path, SysSolverEigCGParams& param);

  /*! \ingroup invert */
  void write(XMLWriter& xml, const string& path, const SysSolverEigCGParams& param);

} // End namespace

#endif 


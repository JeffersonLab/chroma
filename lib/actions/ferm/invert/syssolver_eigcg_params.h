// -*- C++ -*-
// $Id: syssolver_eigcg_params.h,v 1.1 2007-09-25 21:17:12 edwards Exp $
/*! \file
 *  \brief Solve a CG1 system
 */

#ifndef __syssolver_cg_params_h__
#define __syssolver_cg_params_h__

#include "chromabase.h"


namespace Chroma
{

  //! Params for CG inverter
  /*! \ingroup invert */
  struct SysSolverCGParams
  {
    SysSolverCGParams();
    SysSolverCGParams(XMLReader& in, const std::string& path);
    
    Real          RsdCG;           /*!< CG residual */
    int           MaxCG;           /*!< Maximum CG iterations */

    Real          RsdCGRestart;    /*!< CG residual for a possibly double precision restart. Only valid for some solvers eg CG-DWF */

    int           MaxCGRestart;    /*!< Max no of CG iterations for a possibly double precision restart. Only valid for some solvers, eg CG-DWF */
  };


  // Reader/writers
  /*! \ingroup invert */
  void read(XMLReader& xml, const string& path, SysSolverCGParams& param);

  /*! \ingroup invert */
  void write(XMLWriter& xml, const string& path, const SysSolverCGParams& param);

} // End namespace

#endif 


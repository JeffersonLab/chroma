// -*- C++ -*-
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
    int           MinCG;           /*!< Minimum CG iterations (useful for charm) */

    Real          RsdCGRestart;    /*!< CG residual for a possibly double precision restart. Only valid for some solvers eg CG-DWF */

    int           MaxCGRestart;    /*!< Max no of CG iterations for a possibly double precision restart. Only valid for some solvers, eg CG-DWF */
  };


  // Reader/writers
  /*! \ingroup invert */
  void read(XMLReader& xml, const std::string& path, SysSolverCGParams& param);

  /*! \ingroup invert */
  void write(XMLWriter& xml, const std::string& path, const SysSolverCGParams& param);

} // End namespace

#endif 


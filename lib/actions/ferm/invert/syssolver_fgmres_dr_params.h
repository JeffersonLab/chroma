// -*- C++ -*-
/*! \file
 *  \brief Solve an FGMRESR-DR system
 */

#ifndef __syssolver_fgmres_dr_params_h__
#define __syssolver_fgmres_dr_params_h__

#include "chromabase.h"
#include "io/xml_group_reader.h"

namespace Chroma
{

  //! Params for FGMRESDR inverter
  /*! \ingroup invert */
  struct SysSolverFGMRESDRParams
  {
    SysSolverFGMRESDRParams();
    SysSolverFGMRESDRParams(XMLReader& in, const std::string& path);
    
    Real          RsdTarget;           /*!< Target Residuum */
    int           NKrylov;             /*!< Number of vectors before restart */
    int           NDefl;               /*!< Number of deflation vectors */
    int           MaxIter;             /*!< Total Number of Iterations */
    GroupXML_t    PrecondParams;       /*!< Parameters for a preconditioner */
  };


  // Reader/writers
  /*! \ingroup invert */
  void read(XMLReader& xml, const std::string& path, SysSolverFGMRESDRParams& param);

  /*! \ingroup invert */
  void write(XMLWriter& xml, const std::string& path, const SysSolverFGMRESDRParams& param);

} // End namespace

#endif 


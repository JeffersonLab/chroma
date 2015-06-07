// -*- C++ -*-
/*! \file
 *  \brief Params of MR inverter
 */

#ifndef __multi_syssolver_mr_params_h__
#define __multi_syssolver_mr_params_h__

#include "chromabase.h"

namespace Chroma
{

  //! Params for MR inverter
  /*! \ingroup invert */
  struct MultiSysSolverMRParams
  {
    MultiSysSolverMRParams();
    MultiSysSolverMRParams(XMLReader& in, const std::string& path);
    
    multi1d<Real> RsdCG;           /*!< MR residuals */
    int           MaxCG;           /*!< Maximum MR iterations */
  };


  // Reader/writers
  /*! \ingroup invert */
  void read(XMLReader& xml, const std::string& path, MultiSysSolverMRParams& param);

  /*! \ingroup invert */
  void write(XMLWriter& xml, const std::string& path, const MultiSysSolverMRParams& param);

} // End namespace

#endif 


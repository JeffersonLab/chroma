// -*- C++ -*-
// $Id: multi_syssolver_mr_params.h,v 1.1 2007-04-11 03:41:36 edwards Exp $
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
  void read(XMLReader& xml, const string& path, MultiSysSolverMRParams& param);

  /*! \ingroup invert */
  void write(XMLWriter& xml, const string& path, const MultiSysSolverMRParams& param);

} // End namespace

#endif 


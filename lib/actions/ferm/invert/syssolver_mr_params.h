// -*- C++ -*-
// $Id: syssolver_mr_params.h,v 1.1 2007-04-11 03:42:07 edwards Exp $
/*! \file
 *  \brief Solve a MR system
 */

#ifndef __syssolver_mr_params_h__
#define __syssolver_mr_params_h__

#include "chromabase.h"


namespace Chroma
{

  //! Params for MR inverter
  /*! \ingroup invert */
  struct SysSolverMRParams
  {
    SysSolverMRParams();
    SysSolverMRParams(XMLReader& in, const std::string& path);
    
    Real          MROver;          /*!< MR over-relaxation parameter */
    Real          RsdCG;           /*!< MR residual */
    int           MaxCG;           /*!< Maximum MR iterations */
  };


  // Reader/writers
  /*! \ingroup invert */
  void read(XMLReader& xml, const string& path, SysSolverMRParams& param);

  /*! \ingroup invert */
  void write(XMLWriter& xml, const string& path, const SysSolverMRParams& param);

} // End namespace

#endif 


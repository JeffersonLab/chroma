// -*- C++ -*-
// $Id: syssolver_mr_params.h,v 1.2 2007-05-01 14:39:13 bjoo Exp $
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
    Real          RsdMR;           /*!< MR residual */
    int           MaxMR;           /*!< Maximum MR iterations */
  };


  // Reader/writers
  /*! \ingroup invert */
  void read(XMLReader& xml, const string& path, SysSolverMRParams& param);

  /*! \ingroup invert */
  void write(XMLWriter& xml, const string& path, const SysSolverMRParams& param);

} // End namespace

#endif 


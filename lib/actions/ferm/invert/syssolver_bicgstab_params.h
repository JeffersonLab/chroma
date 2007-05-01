// -*- C++ -*-
// $Id: syssolver_bicgstab_params.h,v 3.1 2007-05-01 12:47:24 bjoo Exp $
/*! \file
 *  \brief Solve a BICGSTAB system
 */

#ifndef __syssolver_bicgstab_params_h__
#define __syssolver_bicgstab_params_h__

#include "chromabase.h"


namespace Chroma
{

  //! Params for BiCGStab inverter
  /*! \ingroup invert */
  struct SysSolverBiCGStabParams
  {
    SysSolverBiCGStabParams();
    SysSolverBiCGStabParams(XMLReader& in, const std::string& path);
    Real          RsdBiCGStab;           /*!< BiCGStab residual */
    int           MaxBiCGStab;           /*!< Maximum BiCGStab iterations */
  };


  // Reader/writers
  /*! \ingroup invert */
  void read(XMLReader& xml, const string& path, SysSolverBiCGStabParams& param);

  /*! \ingroup invert */
  void write(XMLWriter& xml, const string& path, const SysSolverBiCGStabParams& param);

} // End namespace

#endif 


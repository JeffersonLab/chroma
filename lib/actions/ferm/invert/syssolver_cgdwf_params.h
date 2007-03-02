// -*- C++ -*-
// $Id: syssolver_cgdwf_params.h,v 1.1 2007-03-02 20:59:34 bjoo Exp $
/*! \file
 *  \brief Solve a CG1 system
 */

#ifndef __syssolver_cgdwf_params_h__
#define __syssolver_cgdwf_params_h__

#include "chromabase.h"


namespace Chroma
{

  //! Params for CGDWF inverter
  /*! \ingroup invert */
  struct SysSolverCGDWFParams
  {
    SysSolverCGDWFParams();
    SysSolverCGDWFParams(XMLReader& in, const std::string& path);
    
    Real          RsdCGSingle;           /*!< CG residual for single precision solve */
    int           MaxCGSingle;           /*!< Maximum CG iterations for single preicsion solve */

    Real          RsdCGDouble;     /*!< CG residual for double precision solve */
    int           MaxCGDouble;     /*!< Maximum CG iterations for double precision solve */
  };


  // Reader/writers
  /*! \ingroup invert */
  void read(XMLReader& xml, const string& path, SysSolverCGDWFParams& param);

  /*! \ingroup invert */
  void write(XMLWriter& xml, const string& path, const SysSolverCGDWFParams& param);

} // End namespace

#endif 


// -*- C++ -*-
// $Id: syssolver_cg_params.h,v 3.1 2006-07-03 15:26:08 edwards Exp $
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
  };


  // Reader/writers
  /*! \ingroup invert */
  void read(XMLReader& xml, const string& path, SysSolverCGParams& param);

  /*! \ingroup invert */
  void write(XMLWriter& xml, const string& path, const SysSolverCGParams& param);

} // End namespace

#endif 


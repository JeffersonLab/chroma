// -*- C++ -*-
// $Id: multi_syssolver_cg_params.h,v 3.1 2006-07-03 15:26:08 edwards Exp $
/*! \file
 *  \brief Params of CG inverter
 */

#ifndef __multi_syssolver_cg_params_h__
#define __multi_syssolver_cg_params_h__

#include "chromabase.h"


namespace Chroma
{

  //! Params for CG inverter
  /*! \ingroup invert */
  struct MultiSysSolverCGParams
  {
    MultiSysSolverCGParams();
    MultiSysSolverCGParams(XMLReader& in, const std::string& path);
    
    multi1d<Real> RsdCG;           /*!< CG residuals */
    int           MaxCG;           /*!< Maximum CG iterations */
  };


  // Reader/writers
  /*! \ingroup invert */
  void read(XMLReader& xml, const string& path, MultiSysSolverCGParams& param);

  /*! \ingroup invert */
  void write(XMLWriter& xml, const string& path, const MultiSysSolverCGParams& param);

} // End namespace

#endif 


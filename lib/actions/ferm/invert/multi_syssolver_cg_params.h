// -*- C++ -*-
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
  void read(XMLReader& xml, const std::string& path, MultiSysSolverCGParams& param);

  /*! \ingroup invert */
  void write(XMLWriter& xml, const std::string& path, const MultiSysSolverCGParams& param);

} // End namespace

#endif 


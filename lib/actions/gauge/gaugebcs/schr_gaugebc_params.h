// -*- C++ -*-
// $Id: schr_gaugebc_params.h,v 2.1 2006-02-26 03:47:52 edwards Exp $
/*! \file
 *  \brief Schroedinger functional gauge boundary conditions
 */

#ifndef __gaugebc_schr_params_h__
#define __gaugebc_schr_params_h__

#include "chromabase.h"
#include "io/enum_io/enum_gaugebc_io.h"

namespace Chroma 
{ 

  /*! @ingroup gaugebcs */
  struct SchrGaugeBCParams {
    SchrGaugeBCParams();
    SchrGaugeBCParams(XMLReader& xml, const std::string& path);
    SchrFunType SchrFun;
    Real SchrPhiMult;
  };
  
  /*! @ingroup gaugebcs */
  void read(XMLReader& xml, const std::string& path, SchrGaugeBCParams& p);

  /*! @ingroup gaugebcs */
  void write(XMLWriter& xml, const std::string& path, const SchrGaugeBCParams& p);
  
}

#endif

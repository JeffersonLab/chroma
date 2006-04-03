// -*- C++ -*-
// $Id: schr_gaugebc_params.h,v 3.0 2006-04-03 04:58:54 edwards Exp $
/*! \file
 *  \brief Schroedinger functional gauge boundary conditions
 */

#ifndef __gaugebc_schr_params_h__
#define __gaugebc_schr_params_h__

#include "chromabase.h"

namespace Chroma 
{ 

  /*! @ingroup gaugebcs */
  struct SchrGaugeBCParams {
    SchrGaugeBCParams();
    SchrGaugeBCParams(XMLReader& xml, const std::string& path);

    int  loop_extent;    /*!< Maximum loop extent in decay direction */
    Real SchrPhiMult;    /*!< Multiplier of phases on boundaries */
    int  decay_dir;      /*!< Decay direction */
  };
  
  /*! @ingroup gaugebcs */
  void read(XMLReader& xml, const std::string& path, SchrGaugeBCParams& p);

  /*! @ingroup gaugebcs */
  void write(XMLWriter& xml, const std::string& path, const SchrGaugeBCParams& p);
  
}

#endif

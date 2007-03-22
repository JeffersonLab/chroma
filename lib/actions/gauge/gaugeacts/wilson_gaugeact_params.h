// -*- C++ -*-
// $Id: wilson_gaugeact_params.h,v 3.2 2007-03-22 19:52:04 bjoo Exp $
/*! \file
 *  \brief Params for Wilson gauge action
 */
#ifndef WILSON_GAUGEACT_PARAMS_H
#define WILSON_GAUGEACT_PARAMS_H
#include "chromabase.h"
#include "io/aniso_io.h"

namespace Chroma {
  //! Parameter structure
  /*! @ingroup gaugeacts */
  struct WilsonGaugeActParams 
  {
    // Base Constructor
    WilsonGaugeActParams() {} 
    
    // Read params from some root path
    WilsonGaugeActParams(XMLReader& xml_in, const std::string& path);

    Real beta;  
    AnisoParam_t aniso;
  };
  
  /*! @ingroup gaugeacts */
  void read(XMLReader& xml, const string& path, WilsonGaugeActParams& param);
}



#endif

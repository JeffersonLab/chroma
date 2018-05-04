// -*- C++ -*-
/*! \file
 *  \brief Wilson fermion action parameters
 */

#ifndef __wilson_fermact_params_w_h__
#define __wilson_fermact_params_w_h__

#include "chromabase.h"
#include "io/aniso_io.h"

namespace Chroma
{
  //! Params for wilson ferm acts
  /*! \ingroup fermacts */
  struct WilsonFermActParams
  {
    WilsonFermActParams();
    WilsonFermActParams(XMLReader& in, const std::string& path);
    
    Real Mass;
    AnisoParam_t anisoParam;
  };


  // Reader/writers
  /*! \ingroup fermacts */
  void read(XMLReader& xml, const std::string& path, WilsonFermActParams& param);

  /*! \ingroup fermacts */
  void write(XMLWriter& xml, const std::string& path, const WilsonFermActParams& param);

}

#endif

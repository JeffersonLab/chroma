// -*- C++ -*-
/*! \file
 *  \brief AsqTad fermion action parameters
 */

#ifndef __asqtad_fermact_params_w_h__
#define __asqtad_fermact_params_w_h__

#include "chromabase.h"
#include "io/aniso_io.h"

namespace Chroma
{
  //! Params for asqtad ferm acts
  /*! \ingroup fermacts */
  struct AsqtadFermActParams
  {
    AsqtadFermActParams();
    AsqtadFermActParams(XMLReader& in, const std::string& path);
    
    Real Mass;
    Real u0 ;
  };


  // Reader/writers
  /*! \ingroup fermacts */
  void read(XMLReader& xml, const std::string& path, AsqtadFermActParams& param);

  /*! \ingroup fermacts */
  void write(XMLWriter& xml, const std::string& path, const AsqtadFermActParams& param);

}

#endif

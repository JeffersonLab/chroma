// -*- C++ -*-
// $Id: asqtad_fermact_params_s.h,v 3.1 2006-11-16 19:49:33 kostas Exp $
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
  void read(XMLReader& xml, const string& path, AsqtadFermActParams& param);

  /*! \ingroup fermacts */
  void write(XMLWriter& xml, const string& path, const AsqtadFermActParams& param);

}

#endif

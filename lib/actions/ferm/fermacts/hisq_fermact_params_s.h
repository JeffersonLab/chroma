// -*- C++ -*-
// $Id: hisq_fermact_params_s.h,v 1.1 2007-05-09 12:43:20 mcneile Exp $
/*! \file
 *  \brief Hisq fermion action parameters
 */

#ifndef __hisq_fermact_params_w_h__
#define __hisq_fermact_params_w_h__

#include "chromabase.h"
#include "io/aniso_io.h"

namespace Chroma
{
  //! Params for hisq ferm acts
  /*! \ingroup fermacts */
  struct HisqFermActParams
  {
    HisqFermActParams();
    HisqFermActParams(XMLReader& in, const std::string& path);
    
    Real Mass;
    Real u0 ; // normally u0 = 1 for Hisq, leave in for now
  };


  // Reader/writers
  /*! \ingroup fermacts */
  void read(XMLReader& xml, const string& path, HisqFermActParams& param);

  /*! \ingroup fermacts */
  void write(XMLWriter& xml, const string& path, const HisqFermActParams& param);

}

#endif

// -*- C++ -*-
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
    Real epsilon ; // to improve the dispersion relation
  };


  // Reader/writers
  /*! \ingroup fermacts */
  void read(XMLReader& xml, const std::string& path, HisqFermActParams& param);

  /*! \ingroup fermacts */
  void write(XMLWriter& xml, const std::string& path, const HisqFermActParams& param);

}

#endif

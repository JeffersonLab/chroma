// -*- C++ -*-
/*! \file
 *  \brief Parameters for Clover fermion action
 */

#ifndef __clover_back_fermact_params_w_h__
#define __clover_back_fermact_params_w_h__

#include "io/aniso_io.h"
#include "clover_fermact_params_w.h"

namespace Chroma
{
  //! Params for clover ferm acts
  /*! \ingroup fermacts */
  class CloverBackFermActParams: public CloverFermActParams
  {
  public:
    CloverBackFermActParams();
    CloverBackFermActParams(XMLReader& in, const std::string& path) ;
    int gamma ;

  };


  // Reader/writers
  /*! \ingroup fermacts */
  void read(XMLReader& xml, const std::string& path, CloverBackFermActParams& param);

  /*! \ingroup fermacts */
  void write(XMLWriter& xml, const std::string& path, const CloverBackFermActParams& param);
}

#endif

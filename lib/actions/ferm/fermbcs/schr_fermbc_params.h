// -*- C++ -*-
// $Id: schr_fermbc_params.h,v 2.1 2006-02-26 03:47:52 edwards Exp $
/*! \file
 *  \brief Schroedinger fermion bc
 */

#ifndef __schr_fermbc_params_h__
#define __schr_fermbc_params_h__

#include "chromabase.h"
#include "io/enum_io/enum_gaugebc_io.h"

namespace Chroma 
{
  //! Params struct for Schroedinger functional bc
  /*! \ingroup fermbc */
  struct SchrFermBCParams 
  {
    SchrFermBCParams() {}
    SchrFermBCParams(XMLReader& in, const std::string& path);

    Real  theta;
    SchrFunType SchrFun;
    Real SchrPhiMult;
  };
  
  // Readers writers
  void read(XMLReader& xml, const std::string& path, SchrFermBCParams& param);
  void write(XMLWriter& xml, const std::string& path, const SchrFermBCParams& param);


} // Namespace Chroma 

// End of include guard 
#endif

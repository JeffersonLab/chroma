// -*- C++ -*-
// $Id: clover_fermact_params_w.h,v 2.1 2005-12-03 21:19:38 edwards Exp $
/*! \file
 *  \brief Parameters for Clover fermion action
 */

#ifndef __clover_fermact_params_w_h__
#define __clover_fermact_params_w_h__

#include "io/aniso_io.h"

namespace Chroma
{
  //! Params for clover ferm acts
  /*! \ingroup fermacts */
  struct CloverFermActParams
  {
    CloverFermActParams();
    CloverFermActParams(XMLReader& in, const std::string& path);
    
    Real Mass;
    Real clovCoeffR;
    Real clovCoeffT;
    Real u0;

    AnisoParam_t anisoParam;
  };


  // Reader/writers
  /*! \ingroup fermacts */
  void read(XMLReader& xml, const string& path, CloverFermActParams& param);

  /*! \ingroup fermacts */
  void write(XMLWriter& xml, const string& path, const CloverFermActParams& param);
}

#endif

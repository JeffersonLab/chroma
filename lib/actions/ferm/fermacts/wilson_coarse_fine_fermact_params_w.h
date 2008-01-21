// -*- C++ -*-
// $Id: wilson_coarse_fine_fermact_params_w.h,v 3.1 2008-01-21 20:18:50 edwards Exp $
/*! \file
 *  \brief Wilson coarse-fine 2+2 anisotropic lattice fermion action parameters
 */

#ifndef __wilson_fermact_coarse_fine_params_w_h__
#define __wilson_fermact_coarse_fine_params_w_h__

#include "chromabase.h"
#include "io/aniso_io.h"

namespace Chroma
{
  //! Params for wilson ferm acts
  /*! \ingroup fermacts */
  struct WilsonCoarseFineFermActParams
  {
    WilsonCoarseFineFermActParams();
    WilsonCoarseFineFermActParams(XMLReader& in, const std::string& path);
    
    multi1d<int> coarse_dirs;
    Real         Mass;
    Real         gamma_f;
  };


  // Reader/writers
  /*! \ingroup fermacts */
  void read(XMLReader& xml, const string& path, WilsonCoarseFineFermActParams& param);

  /*! \ingroup fermacts */
  void write(XMLWriter& xml, const string& path, const WilsonCoarseFineFermActParams& param);

}

#endif

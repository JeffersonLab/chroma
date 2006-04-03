// -*- C++ -*-
// $Id: schr_fermbc_params_w.h,v 3.0 2006-04-03 04:58:48 edwards Exp $
/*! \file
 *  \brief Schroedinger functional ferm boundary conditions
 */

#ifndef __fermbc_schr_params_w_h__
#define __fermbc_schr_params_w_h__

#include "chromabase.h"

namespace Chroma 
{ 

  /*! @ingroup fermbcs */
  struct SchrFermBCParams {
    SchrFermBCParams();
    SchrFermBCParams(XMLReader& xml, const std::string& path);

    int            loop_extent;    /*!< Maximum loop extent in decay direction */
    Real           SchrPhiMult;    /*!< Multiplier of phases on boundaries */
    int            decay_dir;      /*!< Decay direction */
    multi1d<Real>  theta;          /*!< Phase angles on boundaries */
  };
  
  /*! @ingroup fermbcs */
  void read(XMLReader& xml, const std::string& path, SchrFermBCParams& p);

  /*! @ingroup fermbcs */
  void write(XMLWriter& xml, const std::string& path, const SchrFermBCParams& p);
  
}

#endif

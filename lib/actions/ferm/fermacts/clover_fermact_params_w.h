// -*- C++ -*-
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

    // Optional Anisotropy
    AnisoParam_t anisoParam;
    Real max_norm;
    bool max_norm_usedP;


    // Optional mass twist...
    Real twisted_m;
    bool twisted_m_usedP;

    bool stabilized_wilson;

  };


  // Reader/writers
  /*! \ingroup fermacts */
  void read(XMLReader& xml, const std::string& path, CloverFermActParams& param);

  /*! \ingroup fermacts */
  void write(XMLWriter& xml, const std::string& path, const CloverFermActParams& param);
}

#endif

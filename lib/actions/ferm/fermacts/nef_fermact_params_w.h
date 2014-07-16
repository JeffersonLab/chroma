// -*- C++ -*-
// $Id: nef_fermact_params_w.h,v 3.3 2008-11-10 17:59:07 bjoo Exp $
/*! \file
 *  \brief Parameters for Clover fermion action
 */

#ifndef __nef_fermact_params_w_h__
#define __nef_fermact_params_w_h__

#include "io/aniso_io.h"

namespace Chroma
{
  //! Params for clover ferm acts
  /*! \ingroup fermacts */
  struct NEFFermActParams
  {
    NEFFermActParams();
    NEFFermActParams(XMLReader& in, const std::string& path);
    
    Real Mass;
	Real OverMass;
	multi1d<Real> b5;
	multi1d<Real> c5;
	int N5;
    Real clovCoeffR;
    Real clovCoeffT;
    Real u0;

    // Optional Anisotropy
    AnisoParam_t anisoParam;
    Real max_norm;
    bool max_norm_usedP;

    // Zero point energy
    Real sub_zero;
    bool sub_zero_usedP;

    // Optional mass twist...
    Real twisted_m;
    bool twisted_m_usedP;

  };


  // Reader/writers
  /*! \ingroup fermacts */
  void read(XMLReader& xml, const string& path, NEFFermActParams& param);

  /*! \ingroup fermacts */
  void write(XMLWriter& xml, const string& path, const NEFFermActParams& param);
}

#endif

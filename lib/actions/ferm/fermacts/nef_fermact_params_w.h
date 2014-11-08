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
    
    Real Mass;  /* The fermion Mass m_f */
    Real OverMass; /* The domain wall height. In our convention it is positive */

    multi1d<Real> b5; /* For general construction these need to be arrays */
    multi1d<Real> c5;

    int N5;  /* Length of 5th dim */
  };


  // Reader/writers
  /*! \ingroup fermacts */
  void read(XMLReader& xml, const std::string& path, NEFFermActParams& param);

  /*! \ingroup fermacts */
  void write(XMLWriter& xml, const std::string& path, const NEFFermActParams& param);
}

#endif

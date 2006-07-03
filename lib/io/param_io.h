// -*- C++ -*-
// $Id: param_io.h,v 3.1 2006-07-03 15:26:09 edwards Exp $
/*! \file
 *  \brief Various parameter structs and reader/writers
 */

#ifndef __param_io_h__
#define __param_io_h__

#include "chromabase.h"

// ! now needed

namespace Chroma 
{

  /*!
   * Types and structures
   *
   * \ingroup io
   *
   * @{
   */
  //! Convert a Kappa to a mass
  Real kappaToMass(const Real& Kappa);

  //! Convert a Kappa to a mass
  multi1d<Real> kappaToMass(const multi1d<Real>& Kappa);

  //! Convert a Kappa to a mass
  Real massToKappa(const Real& Mass);

  //! Convert a mass to a Kappa
  multi1d<Real> massToKappa(const multi1d<Real>& Mass);



  /*
   * Input 
   */
  struct IO_version_t
  {
    int version;
  };


  /*! @} */  // end of group io

} //end namespace chroma

#endif

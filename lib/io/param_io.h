// -*- C++ -*-
// $Id: param_io.h,v 2.2 2006-03-21 19:13:03 edwards Exp $
/*! \file
 *  \brief Various parameter structs and reader/writers
 */

#ifndef __param_io_h__
#define __param_io_h__

#include "chromabase.h"
#include "invtype.h"

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


  //---------------------------- Readers -----------------------------
  //! Read inverter parameters
  void read(XMLReader& xml, const string& path, InvertParam_t& param);

  //! Read inverter parameters
  void read(XMLReader& xml, const string& path, MultiInvertParam_t& param);


  //---------------------------- Writers -----------------------------
  //! Write inverter parameters
  void write(XMLWriter& xml, const string& path, const InvertParam_t& param);

  //! Write inverter parameters
  void write(XMLWriter& xml, const string& path, const MultiInvertParam_t& param);

  /*! @} */  // end of group io

} //end namespace chroma

#endif

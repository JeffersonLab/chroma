// -*- C++ -*-
// $Id: smearing_io.h,v 3.0 2006-04-03 04:58:56 edwards Exp $
/*! \file
 *  \brief Smearing parameters
 */

#ifndef __smearing_io_h__
#define __smearing_io_h__

#include "chromabase.h"

//! reading enums
#include "io/enum_io/enum_wvfkind_io.h"

namespace Chroma 
{

  /*!
   * Types and structures
   *
   * \ingroup io
   *
   * @{
   */


  //! Parameters for sources and sinks
  struct SmearingParam_t
  {
    SmearingParam_t();  // default constructor
    ~SmearingParam_t() {}
  
    WvfKind       wvf_kind;
    Real          wvf_param;
    int           wvfIntPar;
  };


  //! Read a smearing param struct
  void read(XMLReader& xml, const string& path, SmearingParam_t& param);

  //! Write a smearing param struct
  void write(XMLWriter& xml, const string& path, const SmearingParam_t& param);

  /*! @} */  // end of group io

} //end namespace chroma

#endif

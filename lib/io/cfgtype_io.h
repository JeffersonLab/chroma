// -*- C++ -*-
// $Id: cfgtype_io.h,v 3.0 2006-04-03 04:58:55 edwards Exp $
/*! \file
 *  \brief Configuration structure IO
 */

#ifndef __cfgtype_io_h__
#define __cfgtype_io_h__

#include "io/enum_io/enum_cfgtype_io.h"

namespace Chroma 
{
  //! Gauge configuration structure
  /*! \ingroup io */
  struct Cfg_t
  {
    CfgType      cfg_type;   // storage order for stored gauge configuration
    string       cfg_file;
  };


  //! Configuration input
  /*! \ingroup io */
  void read(XMLReader& xml, const string& path, Cfg_t& input);

  //! Configuration input
  /*! \ingroup io */
  void write(XMLWriter& xml, const string& path, const Cfg_t& input);

} //end namespace chroma

#endif

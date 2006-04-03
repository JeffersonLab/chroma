// -*- C++ -*-
// $Id: gauge_startup.h,v 3.0 2006-04-03 04:59:12 edwards Exp $
/*! \file
 *  \brief Initialize the gauge fields
 */

#ifndef GAUGE_STARTUP_H
#define GAUGE_STARTUP_H

#include "io/cfgtype_io.h"

namespace Chroma {

  //! Initialize the gauge fields
  /*!
   * \ingroup gauge
   *
   * \param gauge_file_xml  File xml
   * \param gauge_xml       Record xml
   * \param u               Gauge fields
   * \param cfg             Configuration structure
   */
  void gaugeStartup(XMLReader& gauge_file_xml,
		    XMLReader& gauge_xml,
		    multi1d<LatticeColorMatrix>& u,
		    Cfg_t& cfg);

}  // end namespace Chroma

#endif

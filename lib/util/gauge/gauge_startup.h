// -*- C++ -*-
// $Id: gauge_startup.h,v 1.4 2004-04-28 14:31:59 edwards Exp $
/*! \file
 *  \brief Initialize the gauge fields
 */

#ifndef GAUGE_STARTUP_H
#define GAUGE_STARTUP_H

#include "io/param_io.h"

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

#endif

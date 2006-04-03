// -*- C++ -*-
// $Id: gaugebcs.h,v 3.0 2006-04-03 04:58:54 edwards Exp $

/*! \file
 * \brief Gauge boundary conditions
 *
 * Support for gauge boundary conditions
 */

/*! \defgroup gaugebc Gauge boundary conditions
 * \ingroup gaugeact
 *
 * Support for gauge boundary conditions
 */

#ifndef __gaugebcsss_h__
#define __gaugebcsss_h__

#include "gaugebc_factory.h"
#include "gaugebc_aggregate.h"

#include "simple_gaugebc.h"
#include "periodic_gaugebc.h"
#include "schroedinger_gaugebc.h"
#include "schr_sf_gaugebc.h"
#include "schr_triv_gaugebc.h"
#include "schr_nonpert_gaugebc.h"
#include "schr_coupling_gaugebc.h"
#include "schr_chromomag_gaugebc.h"
#include "schr_dirich_gaugebc.h"

#endif

// -*- C++ -*-
// $Id: gauge.h,v 1.5 2005-01-12 04:44:53 edwards Exp $

/*! \file
 * \brief Gauge action support
 *
 * Central include file for all gauge action support
 */

/*! \defgroup gaugeact Gauge action support
 * \ingroup actions
 *
 * Support for construction of gauge actions, including 
 * improved gauge actions
 */

#ifndef __actgaugee_h__
#define __actgaugee_h__

#include "actions/gauge/gaugebc_factory.h"
#include "actions/gauge/gaugebc_periodic.h"
#include "actions/gauge/gaugebc_simple.h"
#include "actions/gauge/gaugebc_schroedinger.h"
#include "actions/gauge/gaugebcs.h"

#include "actions/gauge/gaugeact_factory.h"
#include "actions/gauge/wilson_gaugeact.h"

#include "actions/gauge/plaq_gaugeact.h"
#include "actions/gauge/rect_gaugeact.h"
#include "actions/gauge/pg_gaugeact.h"

#endif



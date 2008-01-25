// -*- C++ -*-
// $Id: gauge.h,v 3.3 2008-01-25 22:22:39 edwards Exp $

/*! \file
 * \brief Include all gauge utility routines
 *
 * Utility routines for simple manipulation of gauge fields.
 */

/*! \defgroup gauge Utility routines for manipulating color matrices
 * \ingroup util
 *
 * Utility routines for simple manipulation of gauge fields.
 */

#ifndef __gauge_h__
#define __gauge_h__

#include "gauge_init.h"
#include "gauge_init_factory.h"
#include "gauge_init_aggregate.h"

#include "szinqio_gauge_init.h"

#include "gauge_startup.h"   // deprecated

#include "stout_utils.h"

#include "conjgauge.h"
#include "constgauge.h"
#include "eesu3.h"
#include "expm12.h"
#include "expsu3.h"
#include "expmat.h"
#include "hotst.h"
#include "reunit.h"
#include "unit_check.h"
#include "rgauge.h"
#include "taproj.h"
#include "sun_proj.h"
#include "su3proj.h"
#include "su2extract.h"
#include "sunfill.h"

#include "gauge_s.h"

#endif



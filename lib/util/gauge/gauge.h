// -*- C++ -*-
// $Id: gauge.h,v 1.10 2004-01-07 13:50:09 bjoo Exp $

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

#include "expm12.h"
#include "hotst.h"
#include "reunit.h"
#include "unit_check.h"
#include "rgauge.h"
#include "taproj.h"
#include "sun_proj.h"
#include "su3proj.h"
#include "su2extract.h"
#include "sunfill.h"

#if CHROMA_BUILD_STAGGERED
#include "gauge_s.h"
#endif

#endif




// $Id: io.h,v 1.17 2004-09-08 02:48:26 edwards Exp $

#ifndef __io_h__
#define __io_h__

/*! \file
 * \brief IO routines
 *
 * Readers and writers of gauge fields and propagators
 */

/*! \defgroup io IO routines
 * \ingroup lib
 *
 * Readers and writers of gauge fields and propagators
 */

#include "readszin.h"
#include "qprop_io.h"
#include "szin_io.h"

#include "writeszin.h"

#include "gauge_io.h"
#include "milc_io.h"
#include "kyugauge_io.h"
#include "readmilc.h"
#include "writemilc.h"

#include "param_io.h"
#if defined CHROMA_BUILD_WILSON
#include "io_w.h"
#elif defined CHROMA_BUILD_STAGGERED
#include "io_s.h"
#endif

#endif

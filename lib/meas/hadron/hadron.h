// -*- C++ -*-
// $Id: hadron.h,v 1.11 2004-03-03 10:43:37 mcneile Exp $

/*! \file
 * \brief Hadronic observables
 *
 * Central include file for all measurements of hadronic observables
 */

/*! \defgroup hadron Hadronic observables
 * \ingroup meas
 *
 * Measure hadronic observables like spectroscopy, form-factors,
 * structure functions. Also source construction routines.
 */

#ifndef __hadron_h__
#define __hadron_h__

#include "srcfil.h"
#include "z2_src.h"
#include "stoch_var.h"

#ifdef CHROMA_BUILD_WILSON
#include "hadron_w.h"
#elif CHROMA_BUILD_STAGGERED
#include "hadron_s.h"
#endif

#endif

// -*- C++ -*-
// $Id: hadron.h,v 1.10 2004-01-07 13:50:08 bjoo Exp $

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

#ifdef CHROMA_BUILD_WILSON
#include "hadron_w.h"
#elif CHROMA_BUILD_STAGGERED
#include "hadron_s.h"
#endif

#endif

//
// $Id: qprop.h,v 1.5 2004-05-06 16:48:21 bjoo Exp $

/*! \file
 * \brief Quark propagator solution routines
 *
 * Routines for computing a quark propagator with various fermion actions
 */

/*! \defgroup qprop Quark propagator solution routines
 * \ingroup actions
 *
 * Routines for computing a quark propagator with various fermion actions
 */

#ifndef __qprop_h__
#define __qprop_h__

#ifdef CHROMA_BUILD_WILSON
#include "qprop_w.h"
#elif defined CHROMA_BUILD_STAGGERED
#endif

#endif

// -*- C++ -*-
// $Id: linop.h,v 1.15 2004-01-07 13:50:08 bjoo Exp $

/*! \file
 * \brief Linear operators
 *
 * Various fermion linear operators
 */

/*! \defgroup linop Linear operators
 * \ingroup actions
 *
 * Various fermion linear operators
 */

#ifndef __linop_h__
#define __linop_h__

#include "lmdagm.h"

#ifdef CHROMA_BUILD_WILSON
#include "linop_w.h"
#elif CHROMA_BUILD_STAGGERED
#include "linop_s.h"
#endif

#endif



// -*- C++ -*-
// $Id: linop.h,v 1.6 2003-09-11 11:12:21 bjoo Exp $

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

// #include "lmpsim_w.h"
#include "lwldslash_w.h"

#ifdef BUILD_SSE_WILSON_DSLASH
#include "lwldslash_w_sse.h" 
#endif

#include "unprec_wilson_linop_w.h"
#include "overlapbu_linop_w.h"

#endif



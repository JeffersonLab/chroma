// -*- C++ -*-
// $Id: fermacts.h,v 1.10 2004-01-07 13:50:07 bjoo Exp $

/*! \file
 * \brief Fermion actions
 *
 * Various fermion actions
 */

/*! \defgroup fermact Fermion actions
 * \ingroup actions
 *
 * Various fermion actions
 */

#ifndef __fermactss_h__
#define __fermactss_h__

#ifdef CHROMA_BUILD_WILSON
#include "fermacts_w.h"
#elif  CHROMA_BUILD_STAGGERED
#include "fermacts_s.h"
#endif

#endif



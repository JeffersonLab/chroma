// -*- C++ -*-
// $Id: ferm.h,v 1.4 2004-05-14 00:22:24 edwards Exp $

/*! \file
 * \brief Utilities for manipulating fermion fields
 *
 * Central include file for all fermion manipulation routines
 */

/*! \defgroup ferm Fermion manipulation routines
 * \ingroup util
 *
 * Central include file for all fermion manipulation routines
 */

#ifndef __ferm_h__
#define __ferm_h__

#include "transf.h"

#if defined CHROMA_BUILD_WILSON
#include "ferm_w.h"
#elif defined CHROMA_BUILD_STAGGERED
#include "ferm_s.h"
#endif

#endif



// -*- C++ -*-
// $Id: eig.h,v 1.6 2004-01-16 15:24:40 kostas Exp $

/*! \file
 * \brief Eigenvalue measurements
 *
 * Central include file for all measurements of eigenvalues
 */

/*! \defgroup eig Eigenvalue measurements
 * \ingroup meas
 *
 * Central include file for all measurements related to calculations
 * of eigenvalues of various linear operators.
 */

#ifndef __eig_h__
#define __eig_h__

//#include "gramschm.h"

#ifdef CHROMA_BUILD_WILSON
#include "eig_w.h"
#elif CHROMA_BUILD_STAGGERED
#include "eig_s.h"
#endif

#endif

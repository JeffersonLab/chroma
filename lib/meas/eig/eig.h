// -*- C++ -*-
// $Id: eig.h,v 2.0 2005-09-25 21:04:33 edwards Exp $

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

#include "gramschm.h"
#include "gramschm_array.h"
#include "sn_jacob.h"
#include "sn_jacob_array.h"
#include "ritz.h"
#include "ritz_array.h"
#include "eig_spec.h"
#include "eig_spec_array.h"


#ifdef CHROMA_BUILD_WILSON
#include "eig_w.h"
#elif defined CHROMA_BUILD_STAGGERED
#include "eig_s.h"
#endif

#endif

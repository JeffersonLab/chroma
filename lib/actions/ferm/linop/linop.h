// -*- C++ -*-
// $Id: linop.h,v 1.14 2003-11-23 06:02:57 edwards Exp $

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

// the file dslash_w.h switches in the right dslash op for
// WilsonDslash
#include "dslash_w.h"

#include "unprec_wilson_linop_w.h"
#include "prec_wilson_linop_w.h"
#include "overlapbu_linop_w.h"

#include "unprec_dwf_linop_w.h"
#include "unprec_dwf_linop_array_w.h"
#include "unprec_ovdwf_linop_array_w.h"
#include "dwffld_w.h"
#include "ldwfdslash_w.h"

#include "prec_dwf_linop_array_w.h"

#include "unprec_ovext_linop_array_w.h"

#include "lmdagm_w.h"

#endif



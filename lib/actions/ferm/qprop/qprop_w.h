//
// $Id: qprop_w.h,v 1.4 2005-01-02 05:21:10 edwards Exp $

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

#ifndef __qprop_w_h__
#define __qprop_w_h__

#include "dwf_quarkprop4_w.h"
#include "nef_quarkprop4_w.h"

#if defined(BUILD_SSE_DWF_CG)
#include "prec_dwf_qprop_array_sse_w.h"
#endif

#endif

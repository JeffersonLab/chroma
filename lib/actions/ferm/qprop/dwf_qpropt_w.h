// $Id: dwf_qpropt_w.h,v 1.1 2005-01-07 05:00:10 edwards Exp $
/*! \file
 * \brief Pick up possibly optimized DWF inverters.
 *
 * Typedefs to pick up possibly optimized DWF inverters
 */

#ifndef DWF_QPROPT_W_H
#define DWF_QPROPT_W_H


// The QDP naive class: PrecFermAct5DQprop
#include "prec_fermact_qprop_array.h"


// The following is an ifdef lis that switches in optimised
// Dslash-es. Currently only optimised inverter is the SSE One;

#if defined(BUILD_SSE_DWF_CG)
// The file defines the SSE Dslash class
// The following typedef switches it in.
#include "prec_dwf_qprop_array_sse_w.h"
typedef SSEDWFQpropT DWFQpropT;

#else

// Bottom line, if no optimised DWF qpropT-s exist then the naive QDP qpropT
// becomes the DWFQpropT
typedef PrecFermAct5DQprop<LatticeFermion, multi1d<LatticeColorMatrix> > DWFQpropT;
#endif


#endif

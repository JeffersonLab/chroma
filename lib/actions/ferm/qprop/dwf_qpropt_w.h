// $Id: dwf_qpropt_w.h,v 2.1 2005-12-30 20:28:44 kostas Exp $
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
// Dslash-es. Currently only optimised inverters are the SSE and ALTICEC 

#if defined(BUILD_SSE_DWF_CG)
// The file defines the SSE Dslash class
// The following typedef switches it in.
#include "prec_dwf_qprop_array_sse_w.h"
namespace Chroma {
typedef SSEDWFQpropT DWFQpropT;
}  // end namespace Chroma

#elif defined(BUILD_ALTIVEC_DWF_CG)
// The file defines the ALTIVEC Dslash class
// The following typedef switches it in.
#include "prec_dwf_qprop_array_altivec_w.h"
namespace Chroma {
typedef ALTIVECDWFQpropT DWFQpropT;
}  // end namespace Chroma

#else

// Bottom line, if no optimised DWF qpropT-s exist then the naive QDP qpropT
// becomes the DWFQpropT
namespace Chroma {
typedef PrecFermAct5DQprop<LatticeFermion, multi1d<LatticeColorMatrix> > DWFQpropT;
}  // end namespace Chroma
#endif


#endif

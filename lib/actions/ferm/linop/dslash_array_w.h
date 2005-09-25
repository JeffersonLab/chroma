// -*- C++ -*-
// $Id: dslash_array_w.h,v 2.0 2005-09-25 21:04:28 edwards Exp $
/*! \file
 *  \brief Include possibly optimized Wilson dslash
 */

#ifndef __dslash_array_w_h
#define __dslash_array_w_h


// The QDP naive dslash class: QDPWilsonDslashArray
#include "lwldslash_array_w.h"


// The following is an ifdef lis that switches in optimised
// Dslash-es. Currently only optimised dslash is the SSE One;

#ifdef BUILD_SSE_WILSON_DSLASH
// The lwldslash_w_sse.h defines the SSE Dslash class
// The following typedef switches it in.
# include "lwldslash_array_sse_w.h"
namespace Chroma {
typedef SSEWilsonDslashArray WilsonDslashArray;
}  // end namespace Chroma


// Many #elif clauses could come in here for other opotimised Dslash-es
#elif defined BUILD_PAB_WILSON_DSLASH
# include "lwldslash_array_pab_w.h"
namespace Chroma {
typedef PABWilsonDslashArray WilsonDslashArray;
}  // end namespace Chroma

#else

// Bottom line, if no optimised Dslash-s exist then the naive QDP Dslash
// becomes the WilsonDslash
namespace Chroma {
typedef QDPWilsonDslashArray WilsonDslashArray;
}  // end namespace Chroma
#endif


#endif

// -*- C++ -*-
// $Id: dslash_w.h,v 2.1 2006-02-16 21:02:58 edwards Exp $
/*! \file
 *  \brief Include possibly optimized Wilson dslash
 */

#ifndef DSLASH_W_H
#define DSLASH_W_H

#include "chroma_config.h"

// The QDP naive dslash class: QDPWilsonDslash
#include "lwldslash_w.h"


// The following is an ifdef lis that switches in optimised
// Dslash-es. Currently only optimised dslash is the SSE One;

#ifdef BUILD_SSE_WILSON_DSLASH
// The lwldslash_w_sse.h defines the SSE Dslash class
// The following typedef switches it in.
# include "lwldslash_w_sse.h"
namespace Chroma {
typedef SSEWilsonDslash WilsonDslash;
}  // end namespace Chroma


// Many #elif clauses could come in here for other opotimised Dslash-es
#elif defined BUILD_PAB_WILSON_DSLASH
# include "lwldslash_w_pab.h"
namespace Chroma {
typedef PABWilsonDslash WilsonDslash;
}  // end namespace Chroma

#else

// Bottom line, if no optimised Dslash-s exist then the naive QDP Dslash
// becomes the WilsonDslash
namespace Chroma {
typedef QDPWilsonDslash WilsonDslash;
}  // end namespace Chroma
#endif


#endif

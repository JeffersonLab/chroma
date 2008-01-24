// -*- C++ -*-
// $Id: dslash_w.h,v 3.4 2008-01-24 15:18:53 edwards Exp $
/*! \file
 *  \brief Include possibly optimized Wilson dslash
 */

#ifndef DSLASH_W_H
#define DSLASH_W_H

#include "chroma_config.h"

// QDP Completely naive Dslash class
#include "lwldslash_w.h"

// The QDP dslash class: QDPWilsonDslashOpt 
// This is 'optimised' with temporaries and fused ops.
#include "lwldslash_qdpopt_w.h"


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
typedef QDPWilsonDslashOpt WilsonDslash;
}  // end namespace Chroma
#endif


// 3D Dslashes
// These guards make sure 3D is only ever considered in the right situations


#include "qdp_config.h"
#if QDP_NS==4
#if QDP_NC==3
#if QDP_ND==4

#include "lwldslash_3d_qdp_w.h"
#ifdef BUILD_SSE_WILSON_DSLASH
#include "lwldslash_3d_sse_w.h"

// For now this is the naive Wilson Dslash but 
// I put in this clause because 
namespace Chroma {
typedef SSEWilsonDslash3D WilsonDslash3D;
}

#else

// Bottom line, if no optimised 3d Dslash-s exist then the naive QDP Dslash3D
// becomes the WilsonDslash
namespace Chroma {
typedef QDPWilsonDslash3D WilsonDslash3D;
}  // end namespace Chroma
#endif


#endif
#endif
#endif


#endif

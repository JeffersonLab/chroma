// -*- C++ -*-
// $Id: dslash_w.h,v 3.7 2009-05-01 15:03:02 bjoo Exp $
/*! \file
 *  \brief Include possibly optimized Wilson dslash
 */

#ifndef DSLASH_W_H
#define DSLASH_W_H

#include "qdp_config.h"
#include "chroma_config.h"

// QDP Completely naive Dslash class
#include "lwldslash_w.h"

// The QDP dslash class: QDPWilsonDslashOpt 
// This is 'optimised' with temporaries and fused ops.
#include "lwldslash_qdpopt_w.h"


// The following is an ifdef lis that switches in optimised
// Dslash-es. Currently only optimised dslash is the SSE One;
#ifdef BUILD_CPP_WILSON_DSLASH

#warning "Using New Dslashen"

#include "lwldslash_w_cppf.h"
#include "lwldslash_w_cppd.h"
namespace Chroma {

typedef CPPWilsonDslashF WilsonDslashF;
typedef CPPWilsonDslashD WilsonDslashD;

#if BASE_PRECISION == 32
typedef CPPWilsonDslashF WilsonDslash;
#else 
typedef CPPWilsonDslashD WilsonDslash;
#endif

} // End Namespace Chroma

#elif defined BUILD_SSE_WILSON_DSLASH
// The lwldslash_w_sse.h defines the SSE Dslash class
// The following typedef switches it in.

# include "lwldslash_w_sse.h"
namespace Chroma {
  typedef SSEWilsonDslash WilsonDslash;
#if BASE_PRECISION == 32
  typedef SSEWilsonDslash WilsonDslashF;

  // original code:
  // typedef QDPWilsonDslashOptD WilsonDslashD
  // disabled for now until spin optimizations restored in QDP++
  typedef QDPWilsonDslashD WilsonDslashD;
#else
  typedef SSEWilsonDslash WilsonDslashD;
  typedef QDPWilsonDslashF WilsonDslashF;
#endif

}  // end namespace Chroma

// Many #elif clauses could come in here for other opotimised Dslash-es
#elif defined BUILD_PAB_WILSON_DSLASH
# include "lwldslash_w_pab.h"
namespace Chroma {

  // Assume a DP build
  typedef PABWilsonDslash WilsonDslash;

  // I should set up both a single and a double prec PAB dslash?
  typedef QDPWilsonDslashF WilsonDslashF;

#ifndef CHROMA_USE_SLOPPY_BAGEL_DSLASH
  // IF we are NOT Sloppy the PABWilsonDslash is the DP guy
  typedef PABWilsonDslash WilsonDslashD;
#else
  // If the Dsalsh is SLoppy it is a single prec dslash but
  // with a double prec exterior... so we fall back to QDP
  // Dslash for the true double
  typedef QDPWilsonDslashD WilsonDslashD;
#endif


}  // end namespace Chroma

#else

// Bottom line, if no optimised Dslash-s exist then the naive QDP Dslash
// becomes the WilsonDslash
namespace Chroma {

  typedef QDPWilsonDslash WilsonDslash;
  typedef QDPWilsonDslashF WilsonDslashF;
  typedef QDPWilsonDslashD WilsonDslashD;

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

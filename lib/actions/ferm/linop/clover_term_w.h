// -*- C++ -*-
// $Id: clover_term_w.h,v 2.1 2005-12-18 23:53:26 edwards Exp $
/*! \file
 *  \brief Include possibly optimized Clover terms
 */

#ifndef __clover_term_w_h__
#define __clover_term_w_h__


// The QDP naive clover term
#include "clover_term_qdp_w.h"


// The following is an ifdef lis that switches in optimised
// terms. Currently only optimised dslash is the SSE One;

#ifdef BUILD_SSE_CLOVER_TERM
// The clover_term_sse_w.h defines the SSE Dslash class
// The following typedef switches it in.
# include "clover_term_sse_w.h"
namespace Chroma {
typedef SSECloverTerm CloverTerm;
}  // end namespace Chroma

#else

// Bottom line, if no optimised Dslash-s exist then the naive QDP Dslash
// becomes the WilsonDslash
namespace Chroma {
typedef QDPCloverTerm CloverTerm;
}  // end namespace Chroma
#endif


#endif

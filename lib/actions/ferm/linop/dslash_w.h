#ifndef DSLASH_W_H
#define DSLASH_W_H


// The QDP naive dslash class: QDPWilsonDslash
#include "lwldslash_w.h"


// The following is an ifdef lis that switches in optimised
// Dslash-es. Currently only optimised dslash is the SSE One;

#ifdef BUILD_SSE_WILSON_DSLASH
// The lwldslash_w_sse.h defines the SSE Dslash class
// The following typedef switches it in.
# include "lwldslash_w_sse.h"
typedef SSEWilsonDslash WilsonDslash;

// Many #elif clauses could come in here for other opotimised Dslash-es

#else

// Bottom line, if no optimised Dslash-s exist then the naive QDP Dslash
// becomes the WilsonDslash
typedef QDPWilsonDslash WilsonDslash;
#endif


#endif

// -*- C++ -*-
// $Id: clover_term_w.h,v 3.1 2007-06-18 19:24:12 bjoo Exp $
/*! \file
 *  \brief Include possibly optimized Clover terms
 */

#ifndef __clover_term_w_h__
#define __clover_term_w_h__

#include "chroma_config.h"

// The QDP naive clover term
#include "clover_term_qdp_w.h"


// The following is an ifdef lis that switches in optimised
// terms. Currently only optimised dslash is the SSE One;

#ifdef BUILD_BAGEL_CLOVER_TERM
# include "clover_term_bagel_clover.h"
namespace Chroma {
typedef BAGELCloverTerm CloverTerm;
}  // end namespace Chroma

#else

// Bottom line, if no optimised Dslash-s exist then the naive QDP Dslash
// becomes the WilsonDslash
namespace Chroma {
typedef QDPCloverTerm CloverTerm;
}  // end namespace Chroma
#endif


#endif

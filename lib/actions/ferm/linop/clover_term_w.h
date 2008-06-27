// -*- C++ -*-
// $Id: clover_term_w.h,v 3.2 2008-06-27 20:24:09 bjoo Exp $
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

#if defined(BUILD_BAGEL_CLOVER_TERM)
# include "clover_term_bagel_clover.h"
namespace Chroma {
typedef BAGELCloverTerm CloverTerm;
}  // end namespace Chroma

#elif defined(BUILD_SSED_CLOVER_TERM)

# include "clover_term_ssed.h"
namespace Chroma {
typedef SSEDCloverTerm CloverTerm;
} 
#else 

// Bottom line, if no optimised Dslash-s exist then the naive QDP Dslash
// becomes the WilsonDslash
namespace Chroma {
typedef QDPCloverTerm CloverTerm;
}  // end namespace Chroma
#endif


#endif

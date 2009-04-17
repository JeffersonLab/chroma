// -*- C++ -*-
// $Id: clover_term_w.h,v 3.3 2009-04-17 02:05:33 bjoo Exp $
/*! \file
 *  \brief Include possibly optimized Clover terms
 */

#ifndef __clover_term_w_h__
#define __clover_term_w_h__

#include "chroma_config.h"
#include "qdp_config.h"

// The QDP naive clover term
#include "clover_term_qdp_w.h"


// The following is an ifdef lis that switches in optimised
// terms. Currently only optimised dslash is the SSE One;

#if defined(BUILD_BAGEL_CLOVER_TERM)
# include "clover_term_bagel_clover.h"
namespace Chroma {
  // Assume DP Build
  typedef BAGELCloverTerm CloverTerm;
  typedef BAGELCloverTerm CloverTermD;
  typedef QDPCloverTermF CloverTermF;

}  // end namespace Chroma

#elif defined(BUILD_SSED_CLOVER_TERM)

# include "clover_term_ssed.h"
namespace Chroma {
 

#if BASE_PRECISION==32
  typedef QDPCloverTerm CloverTerm;
#else
  typedef SSEDCloverTerm CloverTerm;
#endif

  typedef QDPCloverTermF CloverTermF;
  typedef SSEDCloverTerm CloverTermD;

} 
#else 

// Bottom line, if no optimised Dslash-s exist then the naive QDP Dslash
// becomes the WilsonDslash
namespace Chroma {
  typedef QDPCloverTerm CloverTerm;
  typedef QDPCloverTermF CloverTermF;
  typedef QDPCloverTermD CloverTermD;
}  // end namespace Chroma
#endif


#endif

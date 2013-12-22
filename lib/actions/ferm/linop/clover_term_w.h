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
//


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
#include "clover_term_qdp_w.h"
namespace Chroma {
 

#if BASE_PRECISION==32
  typedef QDPCloverTerm CloverTerm;
#else
  typedef SSEDCloverTerm CloverTerm;
  template<typename T,typename U>
  struct CloverTermT {
    typedef QDPCloverTermT<T,U> Type_t;
  };
#endif

  typedef QDPCloverTermF CloverTermF;
  typedef SSEDCloverTerm CloverTermD;

}
#elif defined(BUILD_JIT_CLOVER_TERM)

# include "clover_term_ptx_w.h"
namespace Chroma {
 
  typedef PTXCloverTerm CloverTerm;
  typedef PTXCloverTermF CloverTermF;
  typedef PTXCloverTermD CloverTermD;

  template<typename T,typename U>
  struct CloverTermT {
    typedef PTXCloverTermT<T,U> Type_t;
  };


}
#else 

// Bottom line, if no optimised Dslash-s exist then the naive QDP Dslash
// becomes the WilsonDslash

#include "clover_term_qdp_w.h"
namespace Chroma {
  typedef QDPCloverTerm CloverTerm;
  typedef QDPCloverTermF CloverTermF;
  typedef QDPCloverTermD CloverTermD;

  template<typename T,typename U>
  struct CloverTermT {
    typedef QDPCloverTermT<T,U> Type_t;
  };
}  // end namespace Chroma
#endif


#endif

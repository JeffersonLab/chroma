// -*- C++ -*-
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

#if defined(QDPJIT_IS_QDPJITPTX)
#include "clover_term_ptx_w.h"
namespace Chroma {
  typedef PTXCloverTerm CloverTerm;
  typedef PTXCloverTermF CloverTermF;
  typedef PTXCloverTermD CloverTermD;
  template<typename T,typename U>
  struct CloverTermT {
    typedef PTXCloverTermT<T,U> Type_t;
  };
}
#elif defined(QDPJIT_IS_QDPJITNVVM)
#include "clover_term_nvvm_w.h"
namespace Chroma {
  typedef NVVMCloverTerm CloverTerm;
  typedef NVVMCloverTermF CloverTermF;
  typedef NVVMCloverTermD CloverTermD;
  template<typename T,typename U>
  struct CloverTermT {
    typedef NVVMCloverTermT<T,U> Type_t;
  };
}
#else
#include "clover_term_llvm_w.h"
namespace Chroma {
  typedef LLVMCloverTerm CloverTerm;
  typedef LLVMCloverTermF CloverTermF;
  typedef LLVMCloverTermD CloverTermD;

  template<typename T,typename U>
  using CloverTermT = LLVMCloverTermT<T,U>;

}
#endif

#else 

// Bottom line, if no optimised Dslash-s exist then the naive QDP Dslash
// becomes the WilsonDslash

#include "clover_term_qdp_w.h"
namespace Chroma {
  typedef QDPCloverTerm CloverTerm;
  typedef QDPCloverTermF CloverTermF;
  typedef QDPCloverTermD CloverTermD;

  template<typename T,typename U>
  using CloverTermT = QDPCloverTermT<T,U>;
}  // end namespace Chroma
#endif


#endif

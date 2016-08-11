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
   using CloverTerm = QDPCloverTerm;
#else
   using CloverTerm = SSEDCloverTerm;
#endif

  using CloverTermF = QDPCloverTermF;
  using CloverTermD = SSEDCloverTerm;

  template<typename T, typename U>
  using CloverTermT = QDPCloverTermT<T,U>;
}
#elif defined(BUILD_JIT_CLOVER_TERM)

#if defined(QDPJIT_IS_QDPJITPTX)
#include "clover_term_ptx_w.h"
namespace Chroma {
  using CloverTerm = PTXCloverTerm;
  using CloverTermF = PTXCloverTermF;
  using CloverTermD = PTXCloverTermD;

  template<typename T, typename U>
  using CloverTermT = PTXCloverTermT<T,U>;

}
#elif defined(QDPJIT_IS_QDPJITNVVM)
#include "clover_term_nvvm_w.h"
namespace Chroma {
  using CloverTerm  = NVVMCloverTerm;
  using CloverTermF = NVVMCloverTermF;
  using CloverTermD = NVVMCloverTermD;

  template<typename T, typename U>
  using CloverTermT = NVVMCloverTermT<T,U>;

}
#else
#include "clover_term_llvm_w.h"
namespace Chroma {
  using CloverTerm  = LLVMCloverTerm;
  using CloverTermF = LLVMCloverTermF;
  using CloverTermD = LLVMCloverTermD;

  template<typename T,typename U>
  using CloverTermT = LLVMCloverTermT<T,U>;

}
#endif

#else 

// Bottom line, if no optimised Dslash-s exist then the naive QDP Dslash
// becomes the WilsonDslash

#include "clover_term_qdp_w.h"
namespace Chroma {
  using CloverTerm  = QDPCloverTerm;
  using CloverTermF = QDPCloverTermF;
  using CloverTermD = QDPCloverTermD;

  template<typename T,typename U>
  using CloverTermT = QDPCloverTermT<T,U>;
}  // end namespace Chroma
#endif


#endif

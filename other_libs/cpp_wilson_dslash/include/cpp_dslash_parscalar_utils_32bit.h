#ifndef CPP_DSLASH_PARSCALAR_UTILS_32BIT_H
#define CPP_DSLASH_PARSCALAR_UTILS_32BIT_H

#include "dslash_config.h"
#include "cpp_dslash_types.h"

#ifdef DSLASH_USE_SSE2
#include <xmmintrin.h>
/* SSE DECOMP/RECONS functions */
#include "cpp_dslash_parscalar_decomp_32bit_sse2.h"
#include "cpp_dslash_parscalar_decomp_hvv_32bit_sse2.h"
#include "cpp_dslash_parscalar_mvv_recons_32bit_sse2.h"
#include "cpp_dslash_parscalar_recons_32bit_sse2.h"
#else 
/* C ones */
#include "cpp_dslash_parscalar_decomp_32bit_c.h"
#include "cpp_dslash_parscalar_decomp_hvv_32bit_c.h"
#include "cpp_dslash_parscalar_mvv_recons_32bit_c.h"
#include "cpp_dslash_parscalar_recons_32bit_c.h"

#define ALIGN __attribute__ ((aligned(16)))
#endif

#ifdef DSLASH_PREFETCH

#define PREFETCH(addr,hint)  _mm_prefetch((addr),(hint))

#else

#define PREFETCH(addr,hint)  

#endif



#include <cstring>

namespace CPlusPlusWilsonDslash { 
  using namespace Dslash32BitTypes;

  namespace DslashParscalar32Bit { 

    void decomp_plus(size_t lo,size_t hi, int id, const void *ptr);

    void decomp_hvv_plus(size_t lo,size_t hi, int id, const void *ptr);

    void mvv_recons_plus(size_t lo,size_t hi, int id, const void *ptr);

    void recons_plus(size_t lo,size_t hi, int id, const void *ptr );

    void decomp_minus(size_t lo,size_t hi, int id, const void *ptr );

    void decomp_hvv_minus(size_t lo,size_t hi, int id, const void *ptr );

    void mvv_recons_minus(size_t lo,size_t hi, int id, const void *ptr);
    void recons_minus(size_t lo,size_t hi, int id, const void *ptr );

  }

}

#endif

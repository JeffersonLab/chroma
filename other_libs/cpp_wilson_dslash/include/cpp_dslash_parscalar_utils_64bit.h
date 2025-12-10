#ifndef CPP_DSLASH_PARSCALAR_UTILS_64BIT_H
#define CPP_DSLASH_PARSCALAR_UTILS_64BIT_H

#include "dslash_config.h"
#include "cpp_dslash_types.h"
#include "shift_table_parscalar.h"
#include "dispatch_parscalar.h"

#ifdef DSLASH_USE_SSE2
/* SSE DECOMP/RECONS functions */
#include "cpp_dslash_parscalar_decomp_64bit_sse2.h"
#include "cpp_dslash_parscalar_decomp_hvv_64bit_sse2.h"
#include "cpp_dslash_parscalar_mvv_recons_64bit_sse2.h"
#include "cpp_dslash_parscalar_recons_64bit_sse2.h"
#else
/* C Equiv functions */
#include "cpp_dslash_parscalar_decomp_64bit_c.h"
#include "cpp_dslash_parscalar_decomp_hvv_64bit_c.h"
#include "cpp_dslash_parscalar_mvv_recons_64bit_c.h"
#include  "cpp_dslash_parscalar_recons_64bit_c.h"
#define ALIGN __attribute__ ((aligned(16)))
#endif



namespace CPlusPlusWilsonDslash {
  using namespace Dslash64BitTypes;

  namespace DslashParscalar64Bit {

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

#ifndef SSE_SIGN_64BIT_H
#define SSE_SIGN_64BIT_H

#include <xmmintrin.h>
namespace CPlusPlusWilsonDslash { 

  namespace  DslashParscalar64Bit { 
    union SSEMask { 
      unsigned int a[4];
      __m128d vector;
    };

    union SSEMask2 {
      double a[2];
      __m128d vector;
    };

  }
}
#ifndef ALIGN
#define ALIGN __attribute__ ((aligned(16)))
#endif

#endif

#ifndef SSE_SIGN_32BIT_H
#define SSE_SIGN_32BIT_H


namespace CPlusPlusWilsonDslash { 

  namespace  DslashParscalar32Bit { 

    union SSESign { 
      unsigned int a[4];
      __m128 vector;
    } ;
  }
}
#ifndef ALIGN
#define ALIGN __attribute__ ((aligned(16)))
#endif

#endif

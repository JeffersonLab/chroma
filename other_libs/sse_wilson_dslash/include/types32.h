#ifndef TYPES_32_H
#define TYPES_32_H

#ifdef __cplusplus
extern "C" { 
#endif

#ifndef  __INCLUDED_SSE_ALIGN_H__
#include <sse_align.h>
#endif

  /* now overlays for spinors as arrays or structs. Most important are: 
     u_mat_array - is a single link matrix 
     my_mat_array - is a 4-std::vector of u_mat_array pointers
     half_spinor_array - is a 2 component std::vector of color vectors 
     spinor_array - which is a 4 component std::vector of color vectors */
  
  typedef float u_mat_array[3][3][2]  ALIGN;       /* color color re/im */ 
  typedef float spinor_array[4][3][2] ALIGN;       /* Nspin4 color re/im */
  typedef float halfspinor_array[3][2][2]    ALIGN;    /* Half Spinor re/im,spin2,color */
  typedef float colorvec[3][2];

  typedef u_mat_array (*my_mat_array)[4] ALIGN;  

  #include <xmmintrin.h>
  
#define prefetch_spinor(addr) \
  _mm_prefetch( (((char*)(((unsigned long)(addr))&~0x7f))), _MM_HINT_T0 ); \
  _mm_prefetch( (((char*)(((unsigned long)(addr))&~0x7f))+128), _MM_HINT_T0 )


#define prefetch_spinor_nta(addr) \
  _mm_prefetch( (((char*)(((unsigned long)(addr))&~0x7f))), _MM_HINT_NTA ); \
  _mm_prefetch( (((char*)(((unsigned long)(addr))&~0x7f))+128), _MM_HINT_NTA) 

#define prefetch_single(addr) \
  _mm_prefetch( (((char*)(((unsigned long)(addr))&~0x7f))), _MM_HINT_T0 )

#define prefetch_single_nta(addr) \
  _mm_prefetch( (((char*)(((unsigned long)(addr))&~0x7f))), _MM_HINT_NTA )

#define prefetch_su3(addr) \
  _mm_prefetch( (((char*)(((unsigned long)(addr))&~0x7f))), _MM_HINT_T0 ); \
  _mm_prefetch( (((char*)(((unsigned long)(addr))&~0x7f))+128), _MM_HINT_T0) 

#ifdef __cplusplus
}
#endif

#endif

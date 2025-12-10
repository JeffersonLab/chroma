#ifndef CPP_DSLASH_PARSCALAR_DECOMP_32BIT_SSE2_H
#define CPP_DSLASH_PARSCALAR_DECOMP_32BIT_SSE2_H

#include <xmmintrin.h>
#include <cpp_dslash_types.h>
#include <sse_sign_32bit.h>

using namespace CPlusPlusWilsonDslash::Dslash32BitTypes;

namespace CPlusPlusWilsonDslash { 

  namespace  DslashParscalar32Bit { 


inline
void decomp_gamma0_minus( FourSpinor src, HalfSpinor dst) 
{
  SSESign signs24 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x00000000, 0x80000000 }};

  /* c <-> color, s <-> spin */

  /* Space for upper components */
  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2;

  /* Space for lower components */
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;

  /* Swap upper and lower components */
  /* Compiler should spill, or use 64 bit extras */
  __m128 xmm6;
  __m128 xmm7;
  __m128 xmm8;

  /* Swap upper and lower components */
  /* Compiler should spill, or use 64 bit extras */
  __m128 xmm9;
  __m128 xmm10;
  __m128 xmm11;



  xmm0 = _mm_load_ps(&src[0][0][0]);
  xmm2 = _mm_load_ps(&src[0][2][0]);
  xmm6 = _mm_load_ps(&src[1][1][0]);
  
  xmm3 = _mm_load_ps(&src[2][0][0]);
  xmm5 = _mm_load_ps(&src[2][2][0]);
  xmm7 = _mm_load_ps(&src[3][1][0]);

  xmm1 = _mm_xor_ps(xmm1,xmm1); // This should zero 
  xmm4 = _mm_xor_ps(xmm4,xmm4);

  xmm1 = _mm_movelh_ps(xmm1,xmm6);
  xmm4 = _mm_movelh_ps(xmm4,xmm7);

  xmm1 = _mm_movehl_ps(xmm1, xmm0);
  xmm4 = _mm_movehl_ps(xmm4, xmm3);


  xmm0 = _mm_shuffle_ps(xmm0, xmm2, 0xe4);
  xmm3 = _mm_shuffle_ps(xmm3, xmm5, 0xe4);

  xmm2 = _mm_shuffle_ps(xmm2, xmm6, 0xe4);
  xmm5 = _mm_shuffle_ps(xmm5, xmm7, 0xe4);

 
  /* Swap the lower components  and multiply by -i*/
  xmm6 = _mm_shuffle_ps(xmm3, xmm3, 0x1b);
  xmm7 = _mm_shuffle_ps(xmm4, xmm4, 0x1b);
  xmm8 = _mm_shuffle_ps(xmm5, xmm5, 0x1b);


  xmm9 = _mm_xor_ps(xmm6, signs24.vector);
  xmm10 = _mm_xor_ps(xmm7, signs24.vector);
  xmm11 = _mm_xor_ps(xmm8, signs24.vector);


  /* Add */
  xmm0 = _mm_add_ps(xmm0, xmm9);
  xmm1 = _mm_add_ps(xmm1, xmm10);
  xmm2 = _mm_add_ps(xmm2, xmm11);

  /* Store */
  _mm_store_ps(&dst[0][0][0],xmm0);
  _mm_store_ps(&dst[1][0][0],xmm1);
  _mm_store_ps(&dst[2][0][0],xmm2);
  
  
}

inline
void decomp_gamma1_minus( FourSpinor src, HalfSpinor dst)
{


  SSESign signs34 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x00000000, 0x80000000, 0x80000000 }};

  /* Space for upper components */
  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2;

  /* Space for lower components */
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;

  /* Swap upper and lower components */
  /* Compiler should spill, or use 64 bit extras */
  __m128 xmm6;
  __m128 xmm7;
  __m128 xmm8;

  /* Swap upper and lower components */
  /* Compiler should spill, or use 64 bit extras */
  __m128 xmm9;
  __m128 xmm10;
  __m128 xmm11;

  xmm0 = _mm_load_ps(&src[0][0][0]);
  xmm2 = _mm_load_ps(&src[0][2][0]);
  xmm6 = _mm_load_ps(&src[1][1][0]);
  
  xmm3 = _mm_load_ps(&src[2][0][0]);
  xmm5 = _mm_load_ps(&src[2][2][0]);
  xmm7 = _mm_load_ps(&src[3][1][0]);

  xmm1 = _mm_xor_ps(xmm1,xmm1); // This should zero 
  xmm4 = _mm_xor_ps(xmm4,xmm4);

  xmm1 = _mm_movelh_ps(xmm1,xmm6);
  xmm4 = _mm_movelh_ps(xmm4,xmm7);

  xmm1 = _mm_movehl_ps(xmm1, xmm0);
  xmm4 = _mm_movehl_ps(xmm4, xmm3);


  xmm0 = _mm_shuffle_ps(xmm0, xmm2, 0xe4);
  xmm3 = _mm_shuffle_ps(xmm3, xmm5, 0xe4);

  xmm2 = _mm_shuffle_ps(xmm2, xmm6, 0xe4);
  xmm5 = _mm_shuffle_ps(xmm5, xmm7, 0xe4);

 
  /* Swap the lower components */
  xmm6 = _mm_shuffle_ps(xmm3, xmm3, 0x4e);
  xmm7 = _mm_shuffle_ps(xmm4, xmm4, 0x4e);
  xmm8 = _mm_shuffle_ps(xmm5, xmm5, 0x4e);

  xmm9 = _mm_xor_ps(xmm6, signs34.vector);
  xmm10 = _mm_xor_ps(xmm7, signs34.vector);
  xmm11 = _mm_xor_ps(xmm8, signs34.vector);

  /* Add */
  xmm0 = _mm_add_ps(xmm0, xmm9);
  xmm1 = _mm_add_ps(xmm1, xmm10);
  xmm2 = _mm_add_ps(xmm2, xmm11);

  /* Store */
  _mm_store_ps(&dst[0][0][0],xmm0);
  _mm_store_ps(&dst[1][0][0],xmm1);
  _mm_store_ps(&dst[2][0][0],xmm2);

  

}

inline
void decomp_gamma2_minus( FourSpinor src, HalfSpinor dst) 
{


  SSESign signs23 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x80000000, 0x00000000 }};

  /* Space for upper components */
  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2;

  /* Space for lower components */
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;

  /* Swap upper and lower components */
  /* Compiler should spill, or use 64 bit extras */
  __m128 xmm6;
  __m128 xmm7;
  __m128 xmm8;

  /* Swap upper and lower components */
  /* Compiler should spill, or use 64 bit extras */
  __m128 xmm9;
  __m128 xmm10;
  __m128 xmm11;

  xmm0 = _mm_load_ps(&src[0][0][0]);
  xmm2 = _mm_load_ps(&src[0][2][0]);
  xmm6 = _mm_load_ps(&src[1][1][0]);
  
  xmm3 = _mm_load_ps(&src[2][0][0]);
  xmm5 = _mm_load_ps(&src[2][2][0]);
  xmm7 = _mm_load_ps(&src[3][1][0]);

  xmm1 = _mm_xor_ps(xmm1,xmm1); // This should zero 
  xmm4 = _mm_xor_ps(xmm4,xmm4);

  xmm1 = _mm_movelh_ps(xmm1,xmm6);
  xmm4 = _mm_movelh_ps(xmm4,xmm7);

  xmm1 = _mm_movehl_ps(xmm1, xmm0);
  xmm4 = _mm_movehl_ps(xmm4, xmm3);


  xmm0 = _mm_shuffle_ps(xmm0, xmm2, 0xe4);
  xmm3 = _mm_shuffle_ps(xmm3, xmm5, 0xe4);

  xmm2 = _mm_shuffle_ps(xmm2, xmm6, 0xe4);
  xmm5 = _mm_shuffle_ps(xmm5, xmm7, 0xe4);
 
  /* Swap the lower components */
  xmm6 = _mm_shuffle_ps(xmm3, xmm3, 0xb1);
  xmm7 = _mm_shuffle_ps(xmm4, xmm4, 0xb1);
  xmm8 = _mm_shuffle_ps(xmm5, xmm5, 0xb1);

  xmm9 = _mm_xor_ps(xmm6, signs23.vector);
  xmm10 = _mm_xor_ps(xmm7, signs23.vector);
  xmm11 = _mm_xor_ps(xmm8, signs23.vector);

  /* Add */
  xmm0 = _mm_add_ps(xmm0, xmm9);
  xmm1 = _mm_add_ps(xmm1, xmm10);
  xmm2 = _mm_add_ps(xmm2, xmm11);

  /* Store */
  _mm_store_ps(&dst[0][0][0],xmm0);
  _mm_store_ps(&dst[1][0][0],xmm1);
  _mm_store_ps(&dst[2][0][0],xmm2);

  

}

inline
void decomp_gamma3_minus( FourSpinor src, HalfSpinor dst) 
{


  /* Space for upper components */
  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2;

  /* Space for lower components */
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;

  __m128 xmm6;
  __m128 xmm7;

  xmm0 = _mm_load_ps(&src[0][0][0]);
  xmm2 = _mm_load_ps(&src[0][2][0]);
  xmm6 = _mm_load_ps(&src[1][1][0]);
  
  xmm3 = _mm_load_ps(&src[2][0][0]);
  xmm5 = _mm_load_ps(&src[2][2][0]);
  xmm7 = _mm_load_ps(&src[3][1][0]);

  xmm1 = _mm_xor_ps(xmm1,xmm1); // This should zero 
  xmm4 = _mm_xor_ps(xmm4,xmm4);

  xmm1 = _mm_movelh_ps(xmm1,xmm6);
  xmm4 = _mm_movelh_ps(xmm4,xmm7);

  xmm1 = _mm_movehl_ps(xmm1, xmm0);
  xmm4 = _mm_movehl_ps(xmm4, xmm3);


  xmm0 = _mm_shuffle_ps(xmm0, xmm2, 0xe4);
  xmm3 = _mm_shuffle_ps(xmm3, xmm5, 0xe4);

  xmm2 = _mm_shuffle_ps(xmm2, xmm6, 0xe4);
  xmm5 = _mm_shuffle_ps(xmm5, xmm7, 0xe4);

 
  /* sub */
  xmm0 = _mm_sub_ps(xmm0, xmm3);
  xmm1 = _mm_sub_ps(xmm1, xmm4);
  xmm2 = _mm_sub_ps(xmm2, xmm5);

  /* Store */
  _mm_store_ps(&dst[0][0][0],xmm0);
  _mm_store_ps(&dst[1][0][0],xmm1);
  _mm_store_ps(&dst[2][0][0],xmm2);

}


inline
void decomp_gamma0_plus( FourSpinor src, HalfSpinor dst) 
{
  
  SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};


  /* c <-> color, s <-> spin */

  /* Space for upper components */
  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2;

  /* Space for lower components */
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;

  /* Swap upper and lower components */
  /* Compiler should spill, or use 64 bit extras */
  __m128 xmm6;
  __m128 xmm7;
  __m128 xmm8;

  /* Swap upper and lower components */
  /* Compiler should spill, or use 64 bit extras */
  __m128 xmm9;
  __m128 xmm10;
  __m128 xmm11;



  xmm0 = _mm_load_ps(&src[0][0][0]);
  xmm2 = _mm_load_ps(&src[0][2][0]);
  xmm6 = _mm_load_ps(&src[1][1][0]);
  
  xmm3 = _mm_load_ps(&src[2][0][0]);
  xmm5 = _mm_load_ps(&src[2][2][0]);
  xmm7 = _mm_load_ps(&src[3][1][0]);

  xmm1 = _mm_xor_ps(xmm1,xmm1); // This should zero 
  xmm4 = _mm_xor_ps(xmm4,xmm4);

  xmm1 = _mm_movelh_ps(xmm1,xmm6);
  xmm4 = _mm_movelh_ps(xmm4,xmm7);

  xmm1 = _mm_movehl_ps(xmm1, xmm0);
  xmm4 = _mm_movehl_ps(xmm4, xmm3);


  xmm0 = _mm_shuffle_ps(xmm0, xmm2, 0xe4);
  xmm3 = _mm_shuffle_ps(xmm3, xmm5, 0xe4);

  xmm2 = _mm_shuffle_ps(xmm2, xmm6, 0xe4);
  xmm5 = _mm_shuffle_ps(xmm5, xmm7, 0xe4);

 
  /* Swap the lower components  and multiply by +i*/
  xmm6 = _mm_shuffle_ps(xmm3, xmm3, 0x1b);
  xmm7 = _mm_shuffle_ps(xmm4, xmm4, 0x1b);
  xmm8 = _mm_shuffle_ps(xmm5, xmm5, 0x1b);

  xmm9 = _mm_xor_ps(xmm6, signs13.vector);
  xmm10 = _mm_xor_ps(xmm7, signs13.vector);
  xmm11 = _mm_xor_ps(xmm8, signs13.vector);

  /* Add */
  xmm0 = _mm_add_ps(xmm0, xmm9);
  xmm1 = _mm_add_ps(xmm1, xmm10);
  xmm2 = _mm_add_ps(xmm2, xmm11);

  /* Store */
  _mm_store_ps(&dst[0][0][0],xmm0);
  _mm_store_ps(&dst[1][0][0],xmm1);
  _mm_store_ps(&dst[2][0][0],xmm2);
   

}

inline
void decomp_gamma1_plus( FourSpinor src, HalfSpinor dst) 
{

  SSESign signs12 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x80000000, 0x00000000, 0x00000000 }};


  /* Space for upper components */
  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2;

  /* Space for lower components */
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;

  /* Swap upper and lower components */
  /* Compiler should spill, or use 64 bit extras */
  __m128 xmm6;
  __m128 xmm7;
  __m128 xmm8;

  /* Swap upper and lower components */
  /* Compiler should spill, or use 64 bit extras */
  __m128 xmm9;
  __m128 xmm10;
  __m128 xmm11;


  xmm0 = _mm_load_ps(&src[0][0][0]);
  xmm2 = _mm_load_ps(&src[0][2][0]);
  xmm6 = _mm_load_ps(&src[1][1][0]);
  
  xmm3 = _mm_load_ps(&src[2][0][0]);
  xmm5 = _mm_load_ps(&src[2][2][0]);
  xmm7 = _mm_load_ps(&src[3][1][0]);

  xmm1 = _mm_xor_ps(xmm1,xmm1); // This should zero 
  xmm4 = _mm_xor_ps(xmm4,xmm4);

  xmm1 = _mm_movelh_ps(xmm1,xmm6);
  xmm4 = _mm_movelh_ps(xmm4,xmm7);

  xmm1 = _mm_movehl_ps(xmm1, xmm0);
  xmm4 = _mm_movehl_ps(xmm4, xmm3);


  xmm0 = _mm_shuffle_ps(xmm0, xmm2, 0xe4);
  xmm3 = _mm_shuffle_ps(xmm3, xmm5, 0xe4);

  xmm2 = _mm_shuffle_ps(xmm2, xmm6, 0xe4);
  xmm5 = _mm_shuffle_ps(xmm5, xmm7, 0xe4);

 
  /* Swap the lower components */
  xmm6 = _mm_shuffle_ps(xmm3, xmm3, 0x4e);
  xmm7 = _mm_shuffle_ps(xmm4, xmm4, 0x4e);
  xmm8 = _mm_shuffle_ps(xmm5, xmm5, 0x4e);

  xmm9 = _mm_xor_ps(xmm6, signs12.vector);
  xmm10 = _mm_xor_ps(xmm7, signs12.vector);
  xmm11 = _mm_xor_ps(xmm8, signs12.vector);

  /* Add */
  xmm0 = _mm_add_ps(xmm0, xmm9);
  xmm1 = _mm_add_ps(xmm1, xmm10);
  xmm2 = _mm_add_ps(xmm2, xmm11);

  /* Store */
  _mm_store_ps(&dst[0][0][0],xmm0);
  _mm_store_ps(&dst[1][0][0],xmm1);
  _mm_store_ps(&dst[2][0][0],xmm2);


}

inline
void decomp_gamma2_plus( FourSpinor src, HalfSpinor dst) 
{

  SSESign signs14 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x00000000, 0x80000000 }};



  /* Space for upper components */
  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2;

  /* Space for lower components */
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;

  /* Swap upper and lower components */
  /* Compiler should spill, or use 64 bit extras */
  __m128 xmm6;
  __m128 xmm7;
  __m128 xmm8;

  /* Swap upper and lower components */
  /* Compiler should spill, or use 64 bit extras */
  __m128 xmm9;
  __m128 xmm10;
  __m128 xmm11;


  xmm0 = _mm_load_ps(&src[0][0][0]);
  xmm2 = _mm_load_ps(&src[0][2][0]);
  xmm6 = _mm_load_ps(&src[1][1][0]);
  
  xmm3 = _mm_load_ps(&src[2][0][0]);
  xmm5 = _mm_load_ps(&src[2][2][0]);
  xmm7 = _mm_load_ps(&src[3][1][0]);

  xmm1 = _mm_xor_ps(xmm1,xmm1); // This should zero 
  xmm4 = _mm_xor_ps(xmm4,xmm4);

  xmm1 = _mm_movelh_ps(xmm1,xmm6);
  xmm4 = _mm_movelh_ps(xmm4,xmm7);

  xmm1 = _mm_movehl_ps(xmm1, xmm0);
  xmm4 = _mm_movehl_ps(xmm4, xmm3);


  xmm0 = _mm_shuffle_ps(xmm0, xmm2, 0xe4);
  xmm3 = _mm_shuffle_ps(xmm3, xmm5, 0xe4);

  xmm2 = _mm_shuffle_ps(xmm2, xmm6, 0xe4);
  xmm5 = _mm_shuffle_ps(xmm5, xmm7, 0xe4);
 
  /* Swap the lower components */
  xmm6 = _mm_shuffle_ps(xmm3, xmm3, 0xb1);
  xmm7 = _mm_shuffle_ps(xmm4, xmm4, 0xb1);
  xmm8 = _mm_shuffle_ps(xmm5, xmm5, 0xb1);

  xmm9 = _mm_xor_ps(xmm6, signs14.vector);
  xmm10 = _mm_xor_ps(xmm7, signs14.vector);
  xmm11 = _mm_xor_ps(xmm8, signs14.vector);

  /* Add */
  xmm0 = _mm_add_ps(xmm0, xmm9);
  xmm1 = _mm_add_ps(xmm1, xmm10);
  xmm2 = _mm_add_ps(xmm2, xmm11);

  /* Store */
  _mm_store_ps(&dst[0][0][0],xmm0);
  _mm_store_ps(&dst[1][0][0],xmm1);
  _mm_store_ps(&dst[2][0][0],xmm2);


}

inline
void decomp_gamma3_plus( FourSpinor src, HalfSpinor dst) 
{

  /* Space for upper components */
  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2;

  /* Space for lower components */
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;

  __m128 xmm6;
  __m128 xmm7;


  xmm0 = _mm_load_ps(&src[0][0][0]);
  xmm2 = _mm_load_ps(&src[0][2][0]);
  xmm6 = _mm_load_ps(&src[1][1][0]);
  
  xmm3 = _mm_load_ps(&src[2][0][0]);
  xmm5 = _mm_load_ps(&src[2][2][0]);
  xmm7 = _mm_load_ps(&src[3][1][0]);

  xmm1 = _mm_xor_ps(xmm1,xmm1); // This should zero 
  xmm4 = _mm_xor_ps(xmm4,xmm4);

  xmm1 = _mm_movelh_ps(xmm1,xmm6);
  xmm4 = _mm_movelh_ps(xmm4,xmm7);

  xmm1 = _mm_movehl_ps(xmm1, xmm0);
  xmm4 = _mm_movehl_ps(xmm4, xmm3);


  xmm0 = _mm_shuffle_ps(xmm0, xmm2, 0xe4);
  xmm3 = _mm_shuffle_ps(xmm3, xmm5, 0xe4);

  xmm2 = _mm_shuffle_ps(xmm2, xmm6, 0xe4);
  xmm5 = _mm_shuffle_ps(xmm5, xmm7, 0xe4);

 
  /* sub */
  xmm0 = _mm_add_ps(xmm0, xmm3);
  xmm1 = _mm_add_ps(xmm1, xmm4);
  xmm2 = _mm_add_ps(xmm2, xmm5);

  /* Store */
  _mm_store_ps(&dst[0][0][0],xmm0);
  _mm_store_ps(&dst[1][0][0],xmm1);
  _mm_store_ps(&dst[2][0][0],xmm2);


}

  }
}

#endif

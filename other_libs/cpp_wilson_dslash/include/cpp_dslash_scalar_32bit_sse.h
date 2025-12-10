#ifndef CPP_DSLASH_SCALAR_32BIT_SSE_H
#define CPP_DSLASH_SCALAR_32BIT_SSE_H

#include <iostream>
#include <xmmintrin.h>

#ifndef ALIGN
#define ALIGN __attribute__ ((aligned(16)))
#endif

#include <cpp_dslash_types.h>

namespace CPlusPlusWilsonDslash {

  namespace DslashScalar32Bit {

       /* Thread dispatch function for D */
    void DPsiPlus(size_t lo, size_t hi, int id, const void *ptr);
    
    /* Thread dispatch function for D^\dagger */
    void DPsiMinus(size_t lo, size_t hi, int id, const void *ptr);

       /* Thread dispatch function for D */
    void DPsiPlus3D(size_t lo, size_t hi, int id, const void *ptr);
    
    /* Thread dispatch function for D^\dagger */
    void DPsiMinus3D(size_t lo, size_t hi, int id, const void *ptr);


    typedef union { 
      unsigned int a[4];
      __m128 vector;
    } SSESign;

    using namespace Dslash32BitTypes;
    

    inline
    void dslash_plus_dir0_forward(FourSpinor spinor_in,
				  GaugeMatrix u,
				  HalfSpinor  upper_sum,
				  HalfSpinor  lower_sum   )
      
    {
      SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
      SSESign signs12 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x80000000, 0x00000000, 0x00000000 }};
      SSESign signs14 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x00000000, 0x80000000 }};
      SSESign signs24 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x00000000, 0x80000000 }};
      SSESign signs34 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x00000000, 0x80000000, 0x80000000 }};
      SSESign signs23 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x80000000, 0x00000000 }};
      
      __m128 xmm0 ALIGN;
      __m128 xmm1 ALIGN;
      __m128 xmm2 ALIGN;
      __m128 xmm3 ALIGN;
      __m128 xmm4 ALIGN;
      __m128 xmm5 ALIGN;
      __m128 xmm6 ALIGN;
      __m128 xmm7 ALIGN;
      
      

      /* Component 0 into the low 2 floats */
      xmm0 = _mm_loadl_pi(xmm0, (__m64 *)&spinor_in[0][0][0]);
      xmm1 = _mm_loadl_pi(xmm1, (__m64 *)&spinor_in[0][1][0]);
      xmm2 = _mm_loadl_pi(xmm2, (__m64 *)&spinor_in[0][2][0]);
      
      /* Component 1 into the high 2 floats */
      xmm0 = _mm_loadh_pi(xmm0, (__m64 *)&spinor_in[1][0][0]);
      xmm1 = _mm_loadh_pi(xmm1, (__m64 *)&spinor_in[1][1][0]);
      xmm2 = _mm_loadh_pi(xmm2, (__m64 *)&spinor_in[1][2][0]);

      /* Component 2 into low 2 floats */
      xmm3 = _mm_loadl_pi(xmm3, (__m64 *)&spinor_in[2][0][0]);
      xmm4 = _mm_loadl_pi(xmm4, (__m64 *)&spinor_in[2][1][0]);
      xmm5 = _mm_loadl_pi(xmm5, (__m64 *)&spinor_in[2][2][0]);

      /* Component 3 into the high 2 floats */
      xmm3 = _mm_loadh_pi(xmm3, (__m64 *)&spinor_in[3][0][0]);
      xmm4 = _mm_loadh_pi(xmm4, (__m64 *)&spinor_in[3][1][0]);
      xmm5 = _mm_loadh_pi(xmm5, (__m64 *)&spinor_in[3][2][0]);
      

      /* Spin Projection. Results into xmm0-xmm2 */
      /* gamma0 minus projection */
      xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x1b);
      xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x1b);
      xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x1b);
      
      xmm3 = _mm_xor_ps(signs24.vector, xmm3);
      xmm4 = _mm_xor_ps(signs24.vector, xmm4);
      xmm5 = _mm_xor_ps(signs24.vector, xmm5);
      
      xmm0 = _mm_add_ps(xmm0, xmm3);
      xmm1 = _mm_add_ps(xmm1, xmm4);
      xmm2 = _mm_add_ps(xmm2, xmm5);
      

      /* SU3 Multiply */
      xmm3 = _mm_load_ss(&u[0][0][0]);
      xmm6 = _mm_load_ss(&u[1][0][0]);
      xmm4 = _mm_load_ss(&u[0][1][0]);
      xmm7 = _mm_load_ss(&u[2][1][0]);
      xmm5 = _mm_load_ss(&u[0][2][0]);
      xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
      xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
      xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
      xmm3 = _mm_mul_ps(xmm0,xmm3);
      xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
      xmm6 = _mm_mul_ps(xmm1,xmm6);
      xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
      xmm4 = _mm_mul_ps(xmm0, xmm4);
      xmm3 = _mm_add_ps(xmm6, xmm3);
      xmm7 = _mm_mul_ps(xmm2, xmm7);
      xmm5 = _mm_mul_ps(xmm0, xmm5);
      xmm4 = _mm_add_ps(xmm7, xmm4);
      xmm6 = _mm_load_ss(&u[1][2][0]);
      xmm7 = _mm_load_ss(&u[2][0][0]);
      xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
      xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
      xmm6 = _mm_mul_ps(xmm1, xmm6);
      xmm7 = _mm_mul_ps(xmm2, xmm7);
      xmm5 = _mm_add_ps(xmm6, xmm5);
      xmm3 = _mm_add_ps(xmm7, xmm3);
      xmm6 = _mm_load_ss(&u[1][1][0]);
      xmm7 = _mm_load_ss(&u[2][2][0]);
      xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
      xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
      xmm6 = _mm_mul_ps(xmm1, xmm6);
      xmm7 = _mm_mul_ps(xmm2, xmm7);
      xmm4 = _mm_add_ps(xmm6, xmm4);
      xmm5 = _mm_add_ps(xmm7, xmm5);
      xmm6 = _mm_load_ss( &u[0][0][1] );
      xmm7 = _mm_load_ss( &u[1][1][1] );
      xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
      xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
      xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
      xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
      xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
      xmm0 = _mm_xor_ps(signs13.vector, xmm0);
      xmm1 = _mm_xor_ps(signs13.vector, xmm1);
      xmm2 = _mm_xor_ps(signs13.vector, xmm2);
      xmm6 = _mm_mul_ps(xmm0,xmm6);
      xmm7 = _mm_mul_ps(xmm1,xmm7);
      xmm3 = _mm_add_ps(xmm6,xmm3);
      xmm4 = _mm_add_ps(xmm7,xmm4);
      xmm6 = _mm_load_ss( &u[2][2][1] );
      xmm7 = _mm_load_ss( &u[0][1][1] );
      xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
      xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
      xmm6 = _mm_mul_ps(xmm2, xmm6);
      xmm7 = _mm_mul_ps(xmm0, xmm7);
      xmm5 = _mm_add_ps(xmm6, xmm5);
    xmm4 = _mm_add_ps(xmm7, xmm4);
    xmm6 = _mm_load_ss(&u[1][0][1] );
    xmm7 = _mm_load_ss(&u[0][2][1] );
    xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
    xmm6 = _mm_mul_ps(xmm1, xmm6);
    xmm7 = _mm_mul_ps(xmm0, xmm7);
    xmm3 = _mm_add_ps(xmm6, xmm3);
    xmm5 = _mm_add_ps(xmm7, xmm5);
    xmm0 = _mm_load_ss( &u[2][0][1] );
    xmm6 = _mm_load_ss( &u[1][2][1] );
    xmm7 = _mm_load_ss( &u[2][1][1] );
    xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0x0);
    xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
    xmm0 = _mm_mul_ps(xmm2, xmm0);
    xmm6 = _mm_mul_ps(xmm1, xmm6);
    xmm7 = _mm_mul_ps(xmm2, xmm7);
    xmm3 = _mm_add_ps(xmm0, xmm3);
    xmm5 = _mm_add_ps(xmm6, xmm5);
    xmm4 = _mm_add_ps(xmm7, xmm4);


    /* Reconstruction: Upper components just go
       Reshuffle spin and color indices so we can use movaps
       to store - it is aligned and faster... */

    _mm_store_ps(&upper_sum[0][0][0], xmm3);
    _mm_store_ps(&upper_sum[1][0][0], xmm4);
    _mm_store_ps(&upper_sum[2][0][0], xmm5);
    
    /* Lower components - do projection */
    xmm3 = _mm_shuffle_ps( xmm3, xmm3, 0x1b);
    xmm4 = _mm_shuffle_ps( xmm4, xmm4, 0x1b);
    xmm5 = _mm_shuffle_ps( xmm5, xmm5, 0x1b);


    xmm3 = _mm_xor_ps(signs13.vector, xmm3);
    xmm4 = _mm_xor_ps(signs13.vector, xmm4);
    xmm5 = _mm_xor_ps(signs13.vector, xmm5);

    /* Store */
    _mm_store_ps(&lower_sum[0][0][0], xmm3);
    _mm_store_ps(&lower_sum[1][0][0], xmm4);
    _mm_store_ps(&lower_sum[2][0][0], xmm5);


  }
    
  inline
  void dslash_plus_dir0_backward_add( FourSpinor spinor_in,
				      GaugeMatrix u,
				      HalfSpinor  upper_sum,
				      HalfSpinor  lower_sum   )

  {

    SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
    SSESign signs12 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x80000000, 0x00000000, 0x00000000 }};
    SSESign signs14 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x00000000, 0x80000000 }};
    SSESign signs24 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x00000000, 0x80000000 }};
    SSESign signs34 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x00000000, 0x80000000, 0x80000000 }};
    SSESign signs23 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x80000000, 0x00000000 }};


    __m128 xmm0 ALIGN;
    __m128 xmm1 ALIGN;
    __m128 xmm2 ALIGN;
    __m128 xmm3 ALIGN;
    __m128 xmm4 ALIGN;
    __m128 xmm5 ALIGN;
    __m128 xmm6 ALIGN;
    __m128 xmm7 ALIGN;



    /* Component 0 into the low 2 floats */
    xmm0 = _mm_loadl_pi(xmm0, (__m64 *)&spinor_in[0][0][0]);
    xmm1 = _mm_loadl_pi(xmm1, (__m64 *)&spinor_in[0][1][0]);
    xmm2 = _mm_loadl_pi(xmm2, (__m64 *)&spinor_in[0][2][0]);

    /* Component 1 into the high 2 floats */
    xmm0 = _mm_loadh_pi(xmm0, (__m64 *)&spinor_in[1][0][0]);
    xmm1 = _mm_loadh_pi(xmm1, (__m64 *)&spinor_in[1][1][0]);
    xmm2 = _mm_loadh_pi(xmm2, (__m64 *)&spinor_in[1][2][0]);
    
    /* Component 2 into low 2 floats */
    xmm3 = _mm_loadl_pi(xmm3, (__m64 *)&spinor_in[2][0][0]);
    xmm4 = _mm_loadl_pi(xmm4, (__m64 *)&spinor_in[2][1][0]);
    xmm5 = _mm_loadl_pi(xmm5, (__m64 *)&spinor_in[2][2][0]);

    /* Component 3 into the high 2 floats */
    xmm3 = _mm_loadh_pi(xmm3, (__m64 *)&spinor_in[3][0][0]);
    xmm4 = _mm_loadh_pi(xmm4, (__m64 *)&spinor_in[3][1][0]);
    xmm5 = _mm_loadh_pi(xmm5, (__m64 *)&spinor_in[3][2][0]);

    /* Spin Projection. Results into xmm0-xmm2 */
    /* gamma0 plus projection */
    xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x1b);
    xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x1b);
    xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x1b);

    xmm3 = _mm_xor_ps(signs13.vector, xmm3);
    xmm4 = _mm_xor_ps(signs13.vector, xmm4);
    xmm5 = _mm_xor_ps(signs13.vector, xmm5);

    xmm0 = _mm_add_ps(xmm0, xmm3);
    xmm1 = _mm_add_ps(xmm1, xmm4);
    xmm2 = _mm_add_ps(xmm2, xmm5);


    /* Adj SU(3) multiply */

  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm6 = _mm_load_ss(&u[0][1][0]);
  xmm4 = _mm_load_ss(&u[1][0][0]);
  xmm7 = _mm_load_ss(&u[1][2][0]);
  xmm5 = _mm_load_ss(&u[2][0][0]);
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
  xmm3 = _mm_mul_ps(xmm0,xmm3);
  xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
  xmm6 = _mm_mul_ps(xmm1,xmm6);
  xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
  xmm4 = _mm_mul_ps(xmm0, xmm4);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_mul_ps(xmm0, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[2][1][0]);
  xmm7 = _mm_load_ss(&u[0][2][0]);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm3 = _mm_add_ps(xmm7, xmm3);
  xmm6 = _mm_load_ss(&u[1][1][0]);
  xmm7 = _mm_load_ss(&u[2][2][0]);
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm4 = _mm_add_ps(xmm6, xmm4);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_xor_ps(signs24.vector, xmm0);
  xmm1 = _mm_xor_ps(signs24.vector, xmm1);
  xmm2 = _mm_xor_ps(signs24.vector, xmm2);
  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);
  xmm4 = _mm_add_ps(xmm7,xmm4);
  xmm6 = _mm_load_ss( &u[2][2][1] );
  xmm7 = _mm_load_ss( &u[1][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm2, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[0][1][1] );
  xmm7 = _mm_load_ss(&u[2][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm0 = _mm_load_ss( &u[0][2][1] );
  xmm6 = _mm_load_ss( &u[2][1][1] );
  xmm7 = _mm_load_ss( &u[1][2][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm0 = _mm_mul_ps(xmm2, xmm0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm3 = _mm_add_ps(xmm0, xmm3);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);

  /* Result in      xmm3,4,5 */
  /* END MVV */

  /* Reconstruction */

  /* Load up upper partial sum */
  xmm0 = _mm_load_ps(&upper_sum[0][0][0]);
  xmm1 = _mm_load_ps(&upper_sum[1][0][0]);
  xmm2 = _mm_load_ps(&upper_sum[2][0][0]);

  /* Add upper component */
  xmm0 = _mm_add_ps(xmm0, xmm3);
  xmm1 = _mm_add_ps(xmm1, xmm4);
  xmm2 = _mm_add_ps(xmm2, xmm5);
  
  _mm_store_ps(&upper_sum[0][0][0], xmm0);
  _mm_store_ps(&upper_sum[1][0][0], xmm1);
  _mm_store_ps(&upper_sum[2][0][0], xmm2);
    
  

  /* Lower components - do projection */
  xmm0 = _mm_load_ps(&lower_sum[0][0][0]);
  xmm1 = _mm_load_ps(&lower_sum[1][0][0]);
  xmm2 = _mm_load_ps(&lower_sum[2][0][0]);
  

  xmm3 = _mm_shuffle_ps( xmm3, xmm3, 0x1b);
  xmm4 = _mm_shuffle_ps( xmm4, xmm4, 0x1b);
  xmm5 = _mm_shuffle_ps( xmm5, xmm5, 0x1b);

  xmm3 = _mm_xor_ps(signs24.vector, xmm3);
  xmm4 = _mm_xor_ps(signs24.vector, xmm4);
  xmm5 = _mm_xor_ps(signs24.vector, xmm5);

  xmm0 = _mm_add_ps(xmm0, xmm3);
  xmm1 = _mm_add_ps(xmm1, xmm4);
  xmm2 = _mm_add_ps(xmm2, xmm5);
  
  /* Store */
  _mm_store_ps(&lower_sum[0][0][0], xmm0);
  _mm_store_ps(&lower_sum[1][0][0], xmm1);
  _mm_store_ps(&lower_sum[2][0][0], xmm2);


  }

  inline
  void dslash_plus_dir1_forward_add( FourSpinor  spinor_in,
				     GaugeMatrix  u,
				     HalfSpinor  upper_sum,
				     HalfSpinor  lower_sum   )
  {

    SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
    SSESign signs12 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x80000000, 0x00000000, 0x00000000 }};
    SSESign signs14 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x00000000, 0x80000000 }};
    SSESign signs24 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x00000000, 0x80000000 }};
    SSESign signs34 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x00000000, 0x80000000, 0x80000000 }};
    SSESign signs23 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x80000000, 0x00000000 }};


    __m128 xmm0 ALIGN;
    __m128 xmm1 ALIGN;
    __m128 xmm2 ALIGN;
    __m128 xmm3 ALIGN;
    __m128 xmm4 ALIGN;
    __m128 xmm5 ALIGN;
    __m128 xmm6 ALIGN;
    __m128 xmm7 ALIGN;

  
    /* Component 0 into the low 2 floats */
    xmm0 = _mm_loadl_pi(xmm0, (__m64 *)&spinor_in[0][0][0]);
    xmm1 = _mm_loadl_pi(xmm1, (__m64 *)&spinor_in[0][1][0]);
    xmm2 = _mm_loadl_pi(xmm2, (__m64 *)&spinor_in[0][2][0]);

    /* Component 1 into the high 2 floats */
    xmm0 = _mm_loadh_pi(xmm0, (__m64 *)&spinor_in[1][0][0]);
    xmm1 = _mm_loadh_pi(xmm1, (__m64 *)&spinor_in[1][1][0]);
    xmm2 = _mm_loadh_pi(xmm2, (__m64 *)&spinor_in[1][2][0]);
    
    /* Component 2 into low 2 floats */
    xmm3 = _mm_loadl_pi(xmm3, (__m64 *)&spinor_in[2][0][0]);
    xmm4 = _mm_loadl_pi(xmm4, (__m64 *)&spinor_in[2][1][0]);
    xmm5 = _mm_loadl_pi(xmm5, (__m64 *)&spinor_in[2][2][0]);

    /* Component 3 into the high 2 floats */
    xmm3 = _mm_loadh_pi(xmm3, (__m64 *)&spinor_in[3][0][0]);
    xmm4 = _mm_loadh_pi(xmm4, (__m64 *)&spinor_in[3][1][0]);
    xmm5 = _mm_loadh_pi(xmm5, (__m64 *)&spinor_in[3][2][0]);

    /* Spin Projection. Results into xmm0-xmm2 */
    /* gamma1 minus projection */
    xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x4e);
    xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x4e);
    xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x4e);

    xmm3 = _mm_xor_ps(signs34.vector, xmm3);
    xmm4 = _mm_xor_ps(signs34.vector, xmm4);
    xmm5 = _mm_xor_ps(signs34.vector, xmm5);

    xmm0 = _mm_add_ps(xmm0, xmm3);
    xmm1 = _mm_add_ps(xmm1, xmm4);
    xmm2 = _mm_add_ps(xmm2, xmm5);


    /* SU3 Multiply */
    xmm3 = _mm_load_ss(&u[0][0][0]);
    xmm6 = _mm_load_ss(&u[1][0][0]);
    xmm4 = _mm_load_ss(&u[0][1][0]);
    xmm7 = _mm_load_ss(&u[2][1][0]);
    xmm5 = _mm_load_ss(&u[0][2][0]);
    xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
    xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
    xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
    xmm3 = _mm_mul_ps(xmm0,xmm3);
    xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
    xmm6 = _mm_mul_ps(xmm1,xmm6);
    xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
    xmm4 = _mm_mul_ps(xmm0, xmm4);
    xmm3 = _mm_add_ps(xmm6, xmm3);
    xmm7 = _mm_mul_ps(xmm2, xmm7);
    xmm5 = _mm_mul_ps(xmm0, xmm5);
    xmm4 = _mm_add_ps(xmm7, xmm4);
    xmm6 = _mm_load_ss(&u[1][2][0]);
    xmm7 = _mm_load_ss(&u[2][0][0]);
    xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
    xmm6 = _mm_mul_ps(xmm1, xmm6);
    xmm7 = _mm_mul_ps(xmm2, xmm7);
    xmm5 = _mm_add_ps(xmm6, xmm5);
    xmm3 = _mm_add_ps(xmm7, xmm3);
    xmm6 = _mm_load_ss(&u[1][1][0]);
    xmm7 = _mm_load_ss(&u[2][2][0]);
    xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
    xmm6 = _mm_mul_ps(xmm1, xmm6);
    xmm7 = _mm_mul_ps(xmm2, xmm7);
    xmm4 = _mm_add_ps(xmm6, xmm4);
    xmm5 = _mm_add_ps(xmm7, xmm5);
    xmm6 = _mm_load_ss( &u[0][0][1] );
    xmm7 = _mm_load_ss( &u[1][1][1] );
    xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
    xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
    xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
    xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
    xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
    xmm0 = _mm_xor_ps(signs13.vector, xmm0);
    xmm1 = _mm_xor_ps(signs13.vector, xmm1);
    xmm2 = _mm_xor_ps(signs13.vector, xmm2);
    xmm6 = _mm_mul_ps(xmm0,xmm6);
    xmm7 = _mm_mul_ps(xmm1,xmm7);
    xmm3 = _mm_add_ps(xmm6,xmm3);
    xmm4 = _mm_add_ps(xmm7,xmm4);
    xmm6 = _mm_load_ss( &u[2][2][1] );
    xmm7 = _mm_load_ss( &u[0][1][1] );
    xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
    xmm6 = _mm_mul_ps(xmm2, xmm6);
    xmm7 = _mm_mul_ps(xmm0, xmm7);
    xmm5 = _mm_add_ps(xmm6, xmm5);
    xmm4 = _mm_add_ps(xmm7, xmm4);
    xmm6 = _mm_load_ss(&u[1][0][1] );
    xmm7 = _mm_load_ss(&u[0][2][1] );
    xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
    xmm6 = _mm_mul_ps(xmm1, xmm6);
    xmm7 = _mm_mul_ps(xmm0, xmm7);
    xmm3 = _mm_add_ps(xmm6, xmm3);
    xmm5 = _mm_add_ps(xmm7, xmm5);
    xmm0 = _mm_load_ss( &u[2][0][1] );
    xmm6 = _mm_load_ss( &u[1][2][1] );
    xmm7 = _mm_load_ss( &u[2][1][1] );
    xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0x0);
    xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
    xmm0 = _mm_mul_ps(xmm2, xmm0);
    xmm6 = _mm_mul_ps(xmm1, xmm6);
    xmm7 = _mm_mul_ps(xmm2, xmm7);
    xmm3 = _mm_add_ps(xmm0, xmm3);
    xmm5 = _mm_add_ps(xmm6, xmm5);
    xmm4 = _mm_add_ps(xmm7, xmm4);

    /* Reconstruction: Upper components just go
       Reshuffle spin and color indices so we can use movaps
       to store - it is aligned and faster... */
    xmm0 = _mm_load_ps(&upper_sum[0][0][0]);
    xmm1 = _mm_load_ps(&upper_sum[1][0][0]);
    xmm2 = _mm_load_ps(&upper_sum[2][0][0]);

    xmm0 = _mm_add_ps(xmm0, xmm3);
    xmm1 = _mm_add_ps(xmm1, xmm4);
    xmm2 = _mm_add_ps(xmm2, xmm5);

    _mm_store_ps(&upper_sum[0][0][0], xmm0);
    _mm_store_ps(&upper_sum[1][0][0], xmm1);
    _mm_store_ps(&upper_sum[2][0][0], xmm2);
    
    /* Lower components - do projection */
    xmm0 = _mm_load_ps(&lower_sum[0][0][0]);
    xmm1 = _mm_load_ps(&lower_sum[1][0][0]);
    xmm2 = _mm_load_ps(&lower_sum[2][0][0]);

    /* Gamma_minus 1 reconstruction */

    xmm3 = _mm_shuffle_ps( xmm3, xmm3, 0x4e);
    xmm4 = _mm_shuffle_ps( xmm4, xmm4, 0x4e);
    xmm5 = _mm_shuffle_ps( xmm5, xmm5, 0x4e);

    xmm3 = _mm_xor_ps(signs12.vector, xmm3);
    xmm4 = _mm_xor_ps(signs12.vector, xmm4);
    xmm5 = _mm_xor_ps(signs12.vector, xmm5);

    /* Accumulate */
    xmm0 = _mm_add_ps(xmm0, xmm3);
    xmm1 = _mm_add_ps(xmm1, xmm4);
    xmm2 = _mm_add_ps(xmm2, xmm5);

    /* Store */
    _mm_store_ps(&lower_sum[0][0][0], xmm0);
    _mm_store_ps(&lower_sum[1][0][0], xmm1);
    _mm_store_ps(&lower_sum[2][0][0], xmm2);

  }

  inline
  void dslash_plus_dir1_backward_add( FourSpinor  spinor_in,
				      GaugeMatrix  u,
				      HalfSpinor  upper_sum,
				      HalfSpinor  lower_sum   )
  {

    SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
    SSESign signs12 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x80000000, 0x00000000, 0x00000000 }};
    SSESign signs14 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x00000000, 0x80000000 }};
    SSESign signs24 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x00000000, 0x80000000 }};
    SSESign signs34 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x00000000, 0x80000000, 0x80000000 }};
    SSESign signs23 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x80000000, 0x00000000 }};

    __m128 xmm0 ALIGN;
    __m128 xmm1 ALIGN;
    __m128 xmm2 ALIGN;
    __m128 xmm3 ALIGN;
    __m128 xmm4 ALIGN;
    __m128 xmm5 ALIGN;
    __m128 xmm6 ALIGN;
    __m128 xmm7 ALIGN;



    /* Component 0 into the low 2 floats */
    xmm0 = _mm_loadl_pi(xmm0, (__m64 *)&spinor_in[0][0][0]);
    xmm1 = _mm_loadl_pi(xmm1, (__m64 *)&spinor_in[0][1][0]);
    xmm2 = _mm_loadl_pi(xmm2, (__m64 *)&spinor_in[0][2][0]);

    /* Component 1 into the high 2 floats */
    xmm0 = _mm_loadh_pi(xmm0, (__m64 *)&spinor_in[1][0][0]);
    xmm1 = _mm_loadh_pi(xmm1, (__m64 *)&spinor_in[1][1][0]);
    xmm2 = _mm_loadh_pi(xmm2, (__m64 *)&spinor_in[1][2][0]);
    
    /* Component 2 into low 2 floats */
    xmm3 = _mm_loadl_pi(xmm3, (__m64 *)&spinor_in[2][0][0]);
    xmm4 = _mm_loadl_pi(xmm4, (__m64 *)&spinor_in[2][1][0]);
    xmm5 = _mm_loadl_pi(xmm5, (__m64 *)&spinor_in[2][2][0]);

    /* Component 3 into the high 2 floats */
    xmm3 = _mm_loadh_pi(xmm3, (__m64 *)&spinor_in[3][0][0]);
    xmm4 = _mm_loadh_pi(xmm4, (__m64 *)&spinor_in[3][1][0]);
    xmm5 = _mm_loadh_pi(xmm5, (__m64 *)&spinor_in[3][2][0]);

    /* Spin Projection. Results into xmm0-xmm2 */
    /* gamma1 plus projection */
    xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x4e);
    xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x4e);
    xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x4e);

    xmm3 = _mm_xor_ps(signs12.vector, xmm3);
    xmm4 = _mm_xor_ps(signs12.vector, xmm4);
    xmm5 = _mm_xor_ps(signs12.vector, xmm5);

    xmm0 = _mm_add_ps(xmm0, xmm3);
    xmm1 = _mm_add_ps(xmm1, xmm4);
    xmm2 = _mm_add_ps(xmm2, xmm5);


    /* Adj SU(3) multiply */

  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm6 = _mm_load_ss(&u[0][1][0]);
  xmm4 = _mm_load_ss(&u[1][0][0]);
  xmm7 = _mm_load_ss(&u[1][2][0]);
  xmm5 = _mm_load_ss(&u[2][0][0]);
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
  xmm3 = _mm_mul_ps(xmm0,xmm3);
  xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
  xmm6 = _mm_mul_ps(xmm1,xmm6);
  xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
  xmm4 = _mm_mul_ps(xmm0, xmm4);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_mul_ps(xmm0, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[2][1][0]);
  xmm7 = _mm_load_ss(&u[0][2][0]);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm3 = _mm_add_ps(xmm7, xmm3);
  xmm6 = _mm_load_ss(&u[1][1][0]);
  xmm7 = _mm_load_ss(&u[2][2][0]);
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm4 = _mm_add_ps(xmm6, xmm4);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_xor_ps(signs24.vector, xmm0);
  xmm1 = _mm_xor_ps(signs24.vector, xmm1);
  xmm2 = _mm_xor_ps(signs24.vector, xmm2);
  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);
  xmm4 = _mm_add_ps(xmm7,xmm4);
  xmm6 = _mm_load_ss( &u[2][2][1] );
  xmm7 = _mm_load_ss( &u[1][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm2, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[0][1][1] );
  xmm7 = _mm_load_ss(&u[2][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm0 = _mm_load_ss( &u[0][2][1] );
  xmm6 = _mm_load_ss( &u[2][1][1] );
  xmm7 = _mm_load_ss( &u[1][2][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm0 = _mm_mul_ps(xmm2, xmm0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm3 = _mm_add_ps(xmm0, xmm3);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);

  /* Result in      xmm3,4,5 */
  /* END MVV */

  /* Reconstruction */

  /* Load up upper partial sum */
  xmm0 = _mm_load_ps(&upper_sum[0][0][0]);
  xmm1 = _mm_load_ps(&upper_sum[1][0][0]);
  xmm2 = _mm_load_ps(&upper_sum[2][0][0]);

  /* Add upper component */
  xmm0 = _mm_add_ps(xmm0, xmm3);
  xmm1 = _mm_add_ps(xmm1, xmm4);
  xmm2 = _mm_add_ps(xmm2, xmm5);
  
  _mm_store_ps(&upper_sum[0][0][0], xmm0);
  _mm_store_ps(&upper_sum[1][0][0], xmm1);
  _mm_store_ps(&upper_sum[2][0][0], xmm2);
    
  

  /* Lower components - do projection */
  xmm0 = _mm_load_ps(&lower_sum[0][0][0]);
  xmm1 = _mm_load_ps(&lower_sum[1][0][0]);
  xmm2 = _mm_load_ps(&lower_sum[2][0][0]);
  

  xmm3 = _mm_shuffle_ps( xmm3, xmm3, 0x4e);
  xmm4 = _mm_shuffle_ps( xmm4, xmm4, 0x4e);
  xmm5 = _mm_shuffle_ps( xmm5, xmm5, 0x4e);

  xmm3 = _mm_xor_ps(signs34.vector, xmm3);
  xmm4 = _mm_xor_ps(signs34.vector, xmm4);
  xmm5 = _mm_xor_ps(signs34.vector, xmm5);

  xmm0 = _mm_add_ps(xmm0, xmm3);
  xmm1 = _mm_add_ps(xmm1, xmm4);
  xmm2 = _mm_add_ps(xmm2, xmm5);
  
  /* Store */
  _mm_store_ps(&lower_sum[0][0][0], xmm0);
  _mm_store_ps(&lower_sum[1][0][0], xmm1);
  _mm_store_ps(&lower_sum[2][0][0], xmm2);


  }

  inline
  void dslash_plus_dir2_forward_add( FourSpinor  spinor_in,
				     GaugeMatrix  u,
				     HalfSpinor  upper_sum,
				     HalfSpinor  lower_sum   )
  {

    SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
    SSESign signs12 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x80000000, 0x00000000, 0x00000000 }};
    SSESign signs14 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x00000000, 0x80000000 }};
    SSESign signs24 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x00000000, 0x80000000 }};
    SSESign signs34 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x00000000, 0x80000000, 0x80000000 }};
    SSESign signs23 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x80000000, 0x00000000 }};

    __m128 xmm0 ALIGN;
    __m128 xmm1 ALIGN;
    __m128 xmm2 ALIGN;
    __m128 xmm3 ALIGN;
    __m128 xmm4 ALIGN;
    __m128 xmm5 ALIGN;
    __m128 xmm6 ALIGN;
    __m128 xmm7 ALIGN;

  
    /* Component 0 into the low 2 floats */
    xmm0 = _mm_loadl_pi(xmm0, (__m64 *)&spinor_in[0][0][0]);
    xmm1 = _mm_loadl_pi(xmm1, (__m64 *)&spinor_in[0][1][0]);
    xmm2 = _mm_loadl_pi(xmm2, (__m64 *)&spinor_in[0][2][0]);

    /* Component 1 into the high 2 floats */
    xmm0 = _mm_loadh_pi(xmm0, (__m64 *)&spinor_in[1][0][0]);
    xmm1 = _mm_loadh_pi(xmm1, (__m64 *)&spinor_in[1][1][0]);
    xmm2 = _mm_loadh_pi(xmm2, (__m64 *)&spinor_in[1][2][0]);
    
    /* Component 2 into low 2 floats */
    xmm3 = _mm_loadl_pi(xmm3, (__m64 *)&spinor_in[2][0][0]);
    xmm4 = _mm_loadl_pi(xmm4, (__m64 *)&spinor_in[2][1][0]);
    xmm5 = _mm_loadl_pi(xmm5, (__m64 *)&spinor_in[2][2][0]);

    /* Component 3 into the high 2 floats */
    xmm3 = _mm_loadh_pi(xmm3, (__m64 *)&spinor_in[3][0][0]);
    xmm4 = _mm_loadh_pi(xmm4, (__m64 *)&spinor_in[3][1][0]);
    xmm5 = _mm_loadh_pi(xmm5, (__m64 *)&spinor_in[3][2][0]);

    /* Spin Projection. Results into xmm0-xmm2 */
    /* gamma1 minus projection */
    xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0xb1);
    xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0xb1);
    xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0xb1);

    xmm3 = _mm_xor_ps(signs23.vector, xmm3);
    xmm4 = _mm_xor_ps(signs23.vector, xmm4);
    xmm5 = _mm_xor_ps(signs23.vector, xmm5);

    xmm0 = _mm_add_ps(xmm0, xmm3);
    xmm1 = _mm_add_ps(xmm1, xmm4);
    xmm2 = _mm_add_ps(xmm2, xmm5);


    /* SU3 Multiply */
    xmm3 = _mm_load_ss(&u[0][0][0]);
    xmm6 = _mm_load_ss(&u[1][0][0]);
    xmm4 = _mm_load_ss(&u[0][1][0]);
    xmm7 = _mm_load_ss(&u[2][1][0]);
    xmm5 = _mm_load_ss(&u[0][2][0]);
    xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
    xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
    xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
    xmm3 = _mm_mul_ps(xmm0,xmm3);
    xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
    xmm6 = _mm_mul_ps(xmm1,xmm6);
    xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
    xmm4 = _mm_mul_ps(xmm0, xmm4);
    xmm3 = _mm_add_ps(xmm6, xmm3);
    xmm7 = _mm_mul_ps(xmm2, xmm7);
    xmm5 = _mm_mul_ps(xmm0, xmm5);
    xmm4 = _mm_add_ps(xmm7, xmm4);
    xmm6 = _mm_load_ss(&u[1][2][0]);
    xmm7 = _mm_load_ss(&u[2][0][0]);
    xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
    xmm6 = _mm_mul_ps(xmm1, xmm6);
    xmm7 = _mm_mul_ps(xmm2, xmm7);
    xmm5 = _mm_add_ps(xmm6, xmm5);
    xmm3 = _mm_add_ps(xmm7, xmm3);
    xmm6 = _mm_load_ss(&u[1][1][0]);
    xmm7 = _mm_load_ss(&u[2][2][0]);
    xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
    xmm6 = _mm_mul_ps(xmm1, xmm6);
    xmm7 = _mm_mul_ps(xmm2, xmm7);
    xmm4 = _mm_add_ps(xmm6, xmm4);
    xmm5 = _mm_add_ps(xmm7, xmm5);
    xmm6 = _mm_load_ss( &u[0][0][1] );
    xmm7 = _mm_load_ss( &u[1][1][1] );
    xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
    xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
    xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
    xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
    xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
    xmm0 = _mm_xor_ps(signs13.vector, xmm0);
    xmm1 = _mm_xor_ps(signs13.vector, xmm1);
    xmm2 = _mm_xor_ps(signs13.vector, xmm2);
    xmm6 = _mm_mul_ps(xmm0,xmm6);
    xmm7 = _mm_mul_ps(xmm1,xmm7);
    xmm3 = _mm_add_ps(xmm6,xmm3);
    xmm4 = _mm_add_ps(xmm7,xmm4);
    xmm6 = _mm_load_ss( &u[2][2][1] );
    xmm7 = _mm_load_ss( &u[0][1][1] );
    xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
    xmm6 = _mm_mul_ps(xmm2, xmm6);
    xmm7 = _mm_mul_ps(xmm0, xmm7);
    xmm5 = _mm_add_ps(xmm6, xmm5);
    xmm4 = _mm_add_ps(xmm7, xmm4);
    xmm6 = _mm_load_ss(&u[1][0][1] );
    xmm7 = _mm_load_ss(&u[0][2][1] );
    xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
    xmm6 = _mm_mul_ps(xmm1, xmm6);
    xmm7 = _mm_mul_ps(xmm0, xmm7);
    xmm3 = _mm_add_ps(xmm6, xmm3);
    xmm5 = _mm_add_ps(xmm7, xmm5);
    xmm0 = _mm_load_ss( &u[2][0][1] );
    xmm6 = _mm_load_ss( &u[1][2][1] );
    xmm7 = _mm_load_ss( &u[2][1][1] );
    xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0x0);
    xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
    xmm0 = _mm_mul_ps(xmm2, xmm0);
    xmm6 = _mm_mul_ps(xmm1, xmm6);
    xmm7 = _mm_mul_ps(xmm2, xmm7);
    xmm3 = _mm_add_ps(xmm0, xmm3);
    xmm5 = _mm_add_ps(xmm6, xmm5);
    xmm4 = _mm_add_ps(xmm7, xmm4);

    /* Reconstruction: Upper components just go
       Reshuffle spin and color indices so we can use movaps
       to store - it is aligned and faster... */
    xmm0 = _mm_load_ps(&upper_sum[0][0][0]);
    xmm1 = _mm_load_ps(&upper_sum[1][0][0]);
    xmm2 = _mm_load_ps(&upper_sum[2][0][0]);

    xmm0 = _mm_add_ps(xmm0, xmm3);
    xmm1 = _mm_add_ps(xmm1, xmm4);
    xmm2 = _mm_add_ps(xmm2, xmm5);

    _mm_store_ps(&upper_sum[0][0][0], xmm0);
    _mm_store_ps(&upper_sum[1][0][0], xmm1);
    _mm_store_ps(&upper_sum[2][0][0], xmm2);
    
    /* Lower components - do projection */
    xmm0 = _mm_load_ps(&lower_sum[0][0][0]);
    xmm1 = _mm_load_ps(&lower_sum[1][0][0]);
    xmm2 = _mm_load_ps(&lower_sum[2][0][0]);

    /* Gamma_minus 1 reconstruction */

    xmm3 = _mm_shuffle_ps( xmm3, xmm3, 0xb1);
    xmm4 = _mm_shuffle_ps( xmm4, xmm4, 0xb1);
    xmm5 = _mm_shuffle_ps( xmm5, xmm5, 0xb1);

    xmm3 = _mm_xor_ps(signs14.vector, xmm3);
    xmm4 = _mm_xor_ps(signs14.vector, xmm4);
    xmm5 = _mm_xor_ps(signs14.vector, xmm5);

    /* Accumulate */
    xmm0 = _mm_add_ps(xmm0, xmm3);
    xmm1 = _mm_add_ps(xmm1, xmm4);
    xmm2 = _mm_add_ps(xmm2, xmm5);

    /* Store */
    _mm_store_ps(&lower_sum[0][0][0], xmm0);
    _mm_store_ps(&lower_sum[1][0][0], xmm1);
    _mm_store_ps(&lower_sum[2][0][0], xmm2);


  }

  inline
  void dslash_plus_dir2_backward_add( FourSpinor  spinor_in,
				      GaugeMatrix  u,
				      HalfSpinor  upper_sum,
				      HalfSpinor  lower_sum   )
  {

    SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
    SSESign signs12 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x80000000, 0x00000000, 0x00000000 }};
    SSESign signs14 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x00000000, 0x80000000 }};
    SSESign signs24 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x00000000, 0x80000000 }};
    SSESign signs34 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x00000000, 0x80000000, 0x80000000 }};
    SSESign signs23 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x80000000, 0x00000000 }};

    __m128 xmm0 ALIGN;
    __m128 xmm1 ALIGN;
    __m128 xmm2 ALIGN;
    __m128 xmm3 ALIGN;
    __m128 xmm4 ALIGN;
    __m128 xmm5 ALIGN;
    __m128 xmm6 ALIGN;
    __m128 xmm7 ALIGN;



    /* Component 0 into the low 2 floats */
    xmm0 = _mm_loadl_pi(xmm0, (__m64 *)&spinor_in[0][0][0]);
    xmm1 = _mm_loadl_pi(xmm1, (__m64 *)&spinor_in[0][1][0]);
    xmm2 = _mm_loadl_pi(xmm2, (__m64 *)&spinor_in[0][2][0]);

    /* Component 1 into the high 2 floats */
    xmm0 = _mm_loadh_pi(xmm0, (__m64 *)&spinor_in[1][0][0]);
    xmm1 = _mm_loadh_pi(xmm1, (__m64 *)&spinor_in[1][1][0]);
    xmm2 = _mm_loadh_pi(xmm2, (__m64 *)&spinor_in[1][2][0]);
    
    /* Component 2 into low 2 floats */
    xmm3 = _mm_loadl_pi(xmm3, (__m64 *)&spinor_in[2][0][0]);
    xmm4 = _mm_loadl_pi(xmm4, (__m64 *)&spinor_in[2][1][0]);
    xmm5 = _mm_loadl_pi(xmm5, (__m64 *)&spinor_in[2][2][0]);

    /* Component 3 into the high 2 floats */
    xmm3 = _mm_loadh_pi(xmm3, (__m64 *)&spinor_in[3][0][0]);
    xmm4 = _mm_loadh_pi(xmm4, (__m64 *)&spinor_in[3][1][0]);
    xmm5 = _mm_loadh_pi(xmm5, (__m64 *)&spinor_in[3][2][0]);

    /* Spin Projection. Results into xmm0-xmm2 */
    /* gamma1 plus projection */
    xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0xb1);
    xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0xb1);
    xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0xb1);

    xmm3 = _mm_xor_ps(signs14.vector, xmm3);
    xmm4 = _mm_xor_ps(signs14.vector, xmm4);
    xmm5 = _mm_xor_ps(signs14.vector, xmm5);

    xmm0 = _mm_add_ps(xmm0, xmm3);
    xmm1 = _mm_add_ps(xmm1, xmm4);
    xmm2 = _mm_add_ps(xmm2, xmm5);


    /* Adj SU(3) multiply */

  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm6 = _mm_load_ss(&u[0][1][0]);
  xmm4 = _mm_load_ss(&u[1][0][0]);
  xmm7 = _mm_load_ss(&u[1][2][0]);
  xmm5 = _mm_load_ss(&u[2][0][0]);
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
  xmm3 = _mm_mul_ps(xmm0,xmm3);
  xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
  xmm6 = _mm_mul_ps(xmm1,xmm6);
  xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
  xmm4 = _mm_mul_ps(xmm0, xmm4);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_mul_ps(xmm0, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[2][1][0]);
  xmm7 = _mm_load_ss(&u[0][2][0]);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm3 = _mm_add_ps(xmm7, xmm3);
  xmm6 = _mm_load_ss(&u[1][1][0]);
  xmm7 = _mm_load_ss(&u[2][2][0]);
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm4 = _mm_add_ps(xmm6, xmm4);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_xor_ps(signs24.vector, xmm0);
  xmm1 = _mm_xor_ps(signs24.vector, xmm1);
  xmm2 = _mm_xor_ps(signs24.vector, xmm2);
  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);
  xmm4 = _mm_add_ps(xmm7,xmm4);
  xmm6 = _mm_load_ss( &u[2][2][1] );
  xmm7 = _mm_load_ss( &u[1][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm2, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[0][1][1] );
  xmm7 = _mm_load_ss(&u[2][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm0 = _mm_load_ss( &u[0][2][1] );
  xmm6 = _mm_load_ss( &u[2][1][1] );
  xmm7 = _mm_load_ss( &u[1][2][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm0 = _mm_mul_ps(xmm2, xmm0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm3 = _mm_add_ps(xmm0, xmm3);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);

  /* Result in      xmm3,4,5 */
  /* END MVV */

  /* Reconstruction */

  /* Load up upper partial sum */
  xmm0 = _mm_load_ps(&upper_sum[0][0][0]);
  xmm1 = _mm_load_ps(&upper_sum[1][0][0]);
  xmm2 = _mm_load_ps(&upper_sum[2][0][0]);

  /* Add upper component */
  xmm0 = _mm_add_ps(xmm0, xmm3);
  xmm1 = _mm_add_ps(xmm1, xmm4);
  xmm2 = _mm_add_ps(xmm2, xmm5);
  
  _mm_store_ps(&upper_sum[0][0][0], xmm0);
  _mm_store_ps(&upper_sum[1][0][0], xmm1);
  _mm_store_ps(&upper_sum[2][0][0], xmm2);
    
  

  /* Lower components - do projection */
  xmm0 = _mm_load_ps(&lower_sum[0][0][0]);
  xmm1 = _mm_load_ps(&lower_sum[1][0][0]);
  xmm2 = _mm_load_ps(&lower_sum[2][0][0]);
  

  xmm3 = _mm_shuffle_ps( xmm3, xmm3, 0xb1);
  xmm4 = _mm_shuffle_ps( xmm4, xmm4, 0xb1);
  xmm5 = _mm_shuffle_ps( xmm5, xmm5, 0xb1);

  xmm3 = _mm_xor_ps(signs23.vector, xmm3);
  xmm4 = _mm_xor_ps(signs23.vector, xmm4);
  xmm5 = _mm_xor_ps(signs23.vector, xmm5);

  xmm0 = _mm_add_ps(xmm0, xmm3);
  xmm1 = _mm_add_ps(xmm1, xmm4);
  xmm2 = _mm_add_ps(xmm2, xmm5);
  
  /* Store */
  _mm_store_ps(&lower_sum[0][0][0], xmm0);
  _mm_store_ps(&lower_sum[1][0][0], xmm1);
  _mm_store_ps(&lower_sum[2][0][0], xmm2);


  }

  inline
  void dslash_plus_dir2_backward_add_store( FourSpinor  spinor_in,
					    GaugeMatrix  u,
					    HalfSpinor  upper_sum,
					    HalfSpinor  lower_sum,
					    FourSpinor spinor_out)
  {
    SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
    SSESign signs12 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x80000000, 0x00000000, 0x00000000 }};
    SSESign signs14 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x00000000, 0x80000000 }};
    SSESign signs24 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x00000000, 0x80000000 }};
    SSESign signs34 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x00000000, 0x80000000, 0x80000000 }};
    SSESign signs23 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x80000000, 0x00000000 }};

    __m128 xmm0 ALIGN;
    __m128 xmm1 ALIGN;
    __m128 xmm2 ALIGN;
    __m128 xmm3 ALIGN;
    __m128 xmm4 ALIGN;
    __m128 xmm5 ALIGN;
    __m128 xmm6 ALIGN;
    __m128 xmm7 ALIGN;



    /* Component 0 into the low 2 floats */
    xmm0 = _mm_loadl_pi(xmm0, (__m64 *)&spinor_in[0][0][0]);
    xmm1 = _mm_loadl_pi(xmm1, (__m64 *)&spinor_in[0][1][0]);
    xmm2 = _mm_loadl_pi(xmm2, (__m64 *)&spinor_in[0][2][0]);

    /* Component 1 into the high 2 floats */
    xmm0 = _mm_loadh_pi(xmm0, (__m64 *)&spinor_in[1][0][0]);
    xmm1 = _mm_loadh_pi(xmm1, (__m64 *)&spinor_in[1][1][0]);
    xmm2 = _mm_loadh_pi(xmm2, (__m64 *)&spinor_in[1][2][0]);
    
    /* Component 2 into low 2 floats */
    xmm3 = _mm_loadl_pi(xmm3, (__m64 *)&spinor_in[2][0][0]);
    xmm4 = _mm_loadl_pi(xmm4, (__m64 *)&spinor_in[2][1][0]);
    xmm5 = _mm_loadl_pi(xmm5, (__m64 *)&spinor_in[2][2][0]);

    /* Component 3 into the high 2 floats */
    xmm3 = _mm_loadh_pi(xmm3, (__m64 *)&spinor_in[3][0][0]);
    xmm4 = _mm_loadh_pi(xmm4, (__m64 *)&spinor_in[3][1][0]);
    xmm5 = _mm_loadh_pi(xmm5, (__m64 *)&spinor_in[3][2][0]);

    /* Spin Projection. Results into xmm0-xmm2 */
    /* gamma1 plus projection */
    xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0xb1);
    xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0xb1);
    xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0xb1);

    xmm3 = _mm_xor_ps(signs14.vector, xmm3);
    xmm4 = _mm_xor_ps(signs14.vector, xmm4);
    xmm5 = _mm_xor_ps(signs14.vector, xmm5);

    xmm0 = _mm_add_ps(xmm0, xmm3);
    xmm1 = _mm_add_ps(xmm1, xmm4);
    xmm2 = _mm_add_ps(xmm2, xmm5);


    /* Adj SU(3) multiply */

  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm6 = _mm_load_ss(&u[0][1][0]);
  xmm4 = _mm_load_ss(&u[1][0][0]);
  xmm7 = _mm_load_ss(&u[1][2][0]);
  xmm5 = _mm_load_ss(&u[2][0][0]);
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
  xmm3 = _mm_mul_ps(xmm0,xmm3);
  xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
  xmm6 = _mm_mul_ps(xmm1,xmm6);
  xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
  xmm4 = _mm_mul_ps(xmm0, xmm4);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_mul_ps(xmm0, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[2][1][0]);
  xmm7 = _mm_load_ss(&u[0][2][0]);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm3 = _mm_add_ps(xmm7, xmm3);
  xmm6 = _mm_load_ss(&u[1][1][0]);
  xmm7 = _mm_load_ss(&u[2][2][0]);
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm4 = _mm_add_ps(xmm6, xmm4);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_xor_ps(signs24.vector, xmm0);
  xmm1 = _mm_xor_ps(signs24.vector, xmm1);
  xmm2 = _mm_xor_ps(signs24.vector, xmm2);
  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);
  xmm4 = _mm_add_ps(xmm7,xmm4);
  xmm6 = _mm_load_ss( &u[2][2][1] );
  xmm7 = _mm_load_ss( &u[1][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm2, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[0][1][1] );
  xmm7 = _mm_load_ss(&u[2][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm0 = _mm_load_ss( &u[0][2][1] );
  xmm6 = _mm_load_ss( &u[2][1][1] );
  xmm7 = _mm_load_ss( &u[1][2][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm0 = _mm_mul_ps(xmm2, xmm0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm3 = _mm_add_ps(xmm0, xmm3);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);

  /* Result in      xmm3,4,5 */
  /* END MVV */

  /* Reconstruction */

  /* Load up upper partial sum */
  xmm0 = _mm_load_ps(&upper_sum[0][0][0]);
  xmm1 = _mm_load_ps(&upper_sum[1][0][0]);
  xmm2 = _mm_load_ps(&upper_sum[2][0][0]);

  /* Add upper component */
  xmm0 = _mm_add_ps(xmm0, xmm3);
  xmm1 = _mm_add_ps(xmm1, xmm4);
  xmm2 = _mm_add_ps(xmm2, xmm5);
  
  /* Was */
  /*
  _mm_store_ps(&upper_sum[0][0][0], xmm0);
  _mm_store_ps(&upper_sum[1][0][0], xmm1);
  _mm_store_ps(&upper_sum[2][0][0], xmm2);
  */

  /* Pair store */
  _mm_storel_pi((__m64 *)&spinor_out[0][0][0], xmm0);
  _mm_storel_pi((__m64 *)&spinor_out[0][1][0], xmm1);
  _mm_storel_pi((__m64 *)&spinor_out[0][2][0], xmm2);
  _mm_storeh_pi((__m64 *)&spinor_out[1][0][0], xmm0);
  _mm_storeh_pi((__m64 *)&spinor_out[1][1][0], xmm1);
  _mm_storeh_pi((__m64 *)&spinor_out[1][2][0], xmm2);

  

  /* Lower components - do projection */
  xmm0 = _mm_load_ps(&lower_sum[0][0][0]);
  xmm1 = _mm_load_ps(&lower_sum[1][0][0]);
  xmm2 = _mm_load_ps(&lower_sum[2][0][0]);
  

  xmm3 = _mm_shuffle_ps( xmm3, xmm3, 0xb1);
  xmm4 = _mm_shuffle_ps( xmm4, xmm4, 0xb1);
  xmm5 = _mm_shuffle_ps( xmm5, xmm5, 0xb1);

  xmm3 = _mm_xor_ps(signs23.vector, xmm3);
  xmm4 = _mm_xor_ps(signs23.vector, xmm4);
  xmm5 = _mm_xor_ps(signs23.vector, xmm5);

  xmm0 = _mm_add_ps(xmm0, xmm3);
  xmm1 = _mm_add_ps(xmm1, xmm4);
  xmm2 = _mm_add_ps(xmm2, xmm5);
  
  
  /* Store */
  /* Was */
  /*
  _mm_store_ps(&lower_sum[0][0][0], xmm0);
  _mm_store_ps(&lower_sum[1][0][0], xmm1);
  _mm_store_ps(&lower_sum[2][0][0], xmm2);
  */
  /* Store */
  _mm_storel_pi((__m64 *)&spinor_out[2][0][0], xmm0);
  _mm_storel_pi((__m64 *)&spinor_out[2][1][0], xmm1);
  _mm_storel_pi((__m64 *)&spinor_out[2][2][0], xmm2);
  _mm_storeh_pi((__m64 *)&spinor_out[3][0][0], xmm0);
  _mm_storeh_pi((__m64 *)&spinor_out[3][1][0], xmm1);
  _mm_storeh_pi((__m64 *)&spinor_out[3][2][0], xmm2);


  }


  inline
  void dslash_plus_dir3_forward_add( FourSpinor  spinor_in,
				     GaugeMatrix  u,
				     HalfSpinor  upper_sum,
				     HalfSpinor  lower_sum)
  {
    SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
    SSESign signs12 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x80000000, 0x00000000, 0x00000000 }};
    SSESign signs14 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x00000000, 0x80000000 }};
    SSESign signs24 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x00000000, 0x80000000 }};
    SSESign signs34 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x00000000, 0x80000000, 0x80000000 }};
    SSESign signs23 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x80000000, 0x00000000 }};

    __m128 xmm0 ALIGN;
    __m128 xmm1 ALIGN;
    __m128 xmm2 ALIGN;
    __m128 xmm3 ALIGN;
    __m128 xmm4 ALIGN;
    __m128 xmm5 ALIGN;
    __m128 xmm6 ALIGN;
    __m128 xmm7 ALIGN;

   
    /* Component 0 into the low 2 floats */
    xmm0 = _mm_loadl_pi(xmm0, (__m64 *)&spinor_in[0][0][0]);
    xmm1 = _mm_loadl_pi(xmm1, (__m64 *)&spinor_in[0][1][0]);
    xmm2 = _mm_loadl_pi(xmm2, (__m64 *)&spinor_in[0][2][0]);

    /* Component 1 into the high 2 floats */
    xmm0 = _mm_loadh_pi(xmm0, (__m64 *)&spinor_in[1][0][0]);
    xmm1 = _mm_loadh_pi(xmm1, (__m64 *)&spinor_in[1][1][0]);
    xmm2 = _mm_loadh_pi(xmm2, (__m64 *)&spinor_in[1][2][0]);
    
    /* Component 2 into low 2 floats */
    xmm3 = _mm_loadl_pi(xmm3, (__m64 *)&spinor_in[2][0][0]);
    xmm4 = _mm_loadl_pi(xmm4, (__m64 *)&spinor_in[2][1][0]);
    xmm5 = _mm_loadl_pi(xmm5, (__m64 *)&spinor_in[2][2][0]);

    /* Component 3 into the high 2 floats */
    xmm3 = _mm_loadh_pi(xmm3, (__m64 *)&spinor_in[3][0][0]);
    xmm4 = _mm_loadh_pi(xmm4, (__m64 *)&spinor_in[3][1][0]);
    xmm5 = _mm_loadh_pi(xmm5, (__m64 *)&spinor_in[3][2][0]);

    /* Spin Projection. Results into xmm0-xmm2 */
    xmm0 = _mm_sub_ps(xmm0, xmm3);
    xmm1 = _mm_sub_ps(xmm1, xmm4);
    xmm2 = _mm_sub_ps(xmm2, xmm5);


    /* SU3 Multiply */
    xmm3 = _mm_load_ss(&u[0][0][0]);
    xmm6 = _mm_load_ss(&u[1][0][0]);
    xmm4 = _mm_load_ss(&u[0][1][0]);
    xmm7 = _mm_load_ss(&u[2][1][0]);
    xmm5 = _mm_load_ss(&u[0][2][0]);
    xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
    xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
    xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
    xmm3 = _mm_mul_ps(xmm0,xmm3);
    xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
    xmm6 = _mm_mul_ps(xmm1,xmm6);
    xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
    xmm4 = _mm_mul_ps(xmm0, xmm4);
    xmm3 = _mm_add_ps(xmm6, xmm3);
    xmm7 = _mm_mul_ps(xmm2, xmm7);
    xmm5 = _mm_mul_ps(xmm0, xmm5);
    xmm4 = _mm_add_ps(xmm7, xmm4);
    xmm6 = _mm_load_ss(&u[1][2][0]);
    xmm7 = _mm_load_ss(&u[2][0][0]);
    xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
    xmm6 = _mm_mul_ps(xmm1, xmm6);
    xmm7 = _mm_mul_ps(xmm2, xmm7);
    xmm5 = _mm_add_ps(xmm6, xmm5);
    xmm3 = _mm_add_ps(xmm7, xmm3);
    xmm6 = _mm_load_ss(&u[1][1][0]);
    xmm7 = _mm_load_ss(&u[2][2][0]);
    xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
    xmm6 = _mm_mul_ps(xmm1, xmm6);
    xmm7 = _mm_mul_ps(xmm2, xmm7);
    xmm4 = _mm_add_ps(xmm6, xmm4);
    xmm5 = _mm_add_ps(xmm7, xmm5);
    xmm6 = _mm_load_ss( &u[0][0][1] );
    xmm7 = _mm_load_ss( &u[1][1][1] );
    xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
    xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
    xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
    xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
    xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
    xmm0 = _mm_xor_ps(signs13.vector, xmm0);
    xmm1 = _mm_xor_ps(signs13.vector, xmm1);
    xmm2 = _mm_xor_ps(signs13.vector, xmm2);
    xmm6 = _mm_mul_ps(xmm0,xmm6);
    xmm7 = _mm_mul_ps(xmm1,xmm7);
    xmm3 = _mm_add_ps(xmm6,xmm3);
    xmm4 = _mm_add_ps(xmm7,xmm4);
    xmm6 = _mm_load_ss( &u[2][2][1] );
    xmm7 = _mm_load_ss( &u[0][1][1] );
    xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
    xmm6 = _mm_mul_ps(xmm2, xmm6);
    xmm7 = _mm_mul_ps(xmm0, xmm7);
    xmm5 = _mm_add_ps(xmm6, xmm5);
    xmm4 = _mm_add_ps(xmm7, xmm4);
    xmm6 = _mm_load_ss(&u[1][0][1] );
    xmm7 = _mm_load_ss(&u[0][2][1] );
    xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
    xmm6 = _mm_mul_ps(xmm1, xmm6);
    xmm7 = _mm_mul_ps(xmm0, xmm7);
    xmm3 = _mm_add_ps(xmm6, xmm3);
    xmm5 = _mm_add_ps(xmm7, xmm5);
    xmm0 = _mm_load_ss( &u[2][0][1] );
    xmm6 = _mm_load_ss( &u[1][2][1] );
    xmm7 = _mm_load_ss( &u[2][1][1] );
    xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0x0);
    xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
    xmm0 = _mm_mul_ps(xmm2, xmm0);
    xmm6 = _mm_mul_ps(xmm1, xmm6);
    xmm7 = _mm_mul_ps(xmm2, xmm7);
    xmm3 = _mm_add_ps(xmm0, xmm3);
    xmm5 = _mm_add_ps(xmm6, xmm5);
    xmm4 = _mm_add_ps(xmm7, xmm4);

    /* Reconstruction: Upper components just go
       Reshuffle spin and color indices so we can use movaps
       to store - it is aligned and faster... */
    xmm0 = _mm_load_ps(&upper_sum[0][0][0]);
    xmm1 = _mm_load_ps(&upper_sum[1][0][0]);
    xmm2 = _mm_load_ps(&upper_sum[2][0][0]);

    xmm0 = _mm_add_ps(xmm0, xmm3);
    xmm1 = _mm_add_ps(xmm1, xmm4);
    xmm2 = _mm_add_ps(xmm2, xmm5);

    _mm_store_ps(&upper_sum[0][0][0], xmm0);
    _mm_store_ps(&upper_sum[1][0][0], xmm1);
    _mm_store_ps(&upper_sum[2][0][0], xmm2);
    
    /* Lower components - do projection */
    xmm0 = _mm_load_ps(&lower_sum[0][0][0]);
    xmm1 = _mm_load_ps(&lower_sum[1][0][0]);
    xmm2 = _mm_load_ps(&lower_sum[2][0][0]);

    /* Accumulate */
    xmm0 = _mm_sub_ps(xmm0, xmm3);
    xmm1 = _mm_sub_ps(xmm1, xmm4);
    xmm2 = _mm_sub_ps(xmm2, xmm5);

    /* Store */
    _mm_store_ps(&lower_sum[0][0][0], xmm0);
    _mm_store_ps(&lower_sum[1][0][0], xmm1);
    _mm_store_ps(&lower_sum[2][0][0], xmm2);



  }

  inline
  void dslash_plus_dir3_backward_add_store( FourSpinor  spinor_in,
					    GaugeMatrix  u,
					    HalfSpinor  upper_sum,
					    HalfSpinor  lower_sum,
					    FourSpinor spinor_out)

  {
    SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
    SSESign signs12 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x80000000, 0x00000000, 0x00000000 }};
    SSESign signs14 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x00000000, 0x80000000 }};
    SSESign signs24 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x00000000, 0x80000000 }};
    SSESign signs34 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x00000000, 0x80000000, 0x80000000 }};
    SSESign signs23 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x80000000, 0x00000000 }};

    __m128 xmm0 ALIGN;
    __m128 xmm1 ALIGN;
    __m128 xmm2 ALIGN;
    __m128 xmm3 ALIGN;
    __m128 xmm4 ALIGN;
    __m128 xmm5 ALIGN;
    __m128 xmm6 ALIGN;
    __m128 xmm7 ALIGN;



    /* Component 0 into the low 2 floats */
    xmm0 = _mm_loadl_pi(xmm0, (__m64 *)&spinor_in[0][0][0]);
    xmm1 = _mm_loadl_pi(xmm1, (__m64 *)&spinor_in[0][1][0]);
    xmm2 = _mm_loadl_pi(xmm2, (__m64 *)&spinor_in[0][2][0]);

    /* Component 1 into the high 2 floats */
    xmm0 = _mm_loadh_pi(xmm0, (__m64 *)&spinor_in[1][0][0]);
    xmm1 = _mm_loadh_pi(xmm1, (__m64 *)&spinor_in[1][1][0]);
    xmm2 = _mm_loadh_pi(xmm2, (__m64 *)&spinor_in[1][2][0]);
    
    /* Component 2 into low 2 floats */
    xmm3 = _mm_loadl_pi(xmm3, (__m64 *)&spinor_in[2][0][0]);
    xmm4 = _mm_loadl_pi(xmm4, (__m64 *)&spinor_in[2][1][0]);
    xmm5 = _mm_loadl_pi(xmm5, (__m64 *)&spinor_in[2][2][0]);

    /* Component 3 into the high 2 floats */
    xmm3 = _mm_loadh_pi(xmm3, (__m64 *)&spinor_in[3][0][0]);
    xmm4 = _mm_loadh_pi(xmm4, (__m64 *)&spinor_in[3][1][0]);
    xmm5 = _mm_loadh_pi(xmm5, (__m64 *)&spinor_in[3][2][0]);

    /* Spin Projection. Results into xmm0-xmm2 */
    xmm0 = _mm_add_ps(xmm0, xmm3);
    xmm1 = _mm_add_ps(xmm1, xmm4);
    xmm2 = _mm_add_ps(xmm2, xmm5);


    /* Adj SU(3) multiply */

  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm6 = _mm_load_ss(&u[0][1][0]);
  xmm4 = _mm_load_ss(&u[1][0][0]);
  xmm7 = _mm_load_ss(&u[1][2][0]);
  xmm5 = _mm_load_ss(&u[2][0][0]);
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
  xmm3 = _mm_mul_ps(xmm0,xmm3);
  xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
  xmm6 = _mm_mul_ps(xmm1,xmm6);
  xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
  xmm4 = _mm_mul_ps(xmm0, xmm4);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_mul_ps(xmm0, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[2][1][0]);
  xmm7 = _mm_load_ss(&u[0][2][0]);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm3 = _mm_add_ps(xmm7, xmm3);
  xmm6 = _mm_load_ss(&u[1][1][0]);
  xmm7 = _mm_load_ss(&u[2][2][0]);
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm4 = _mm_add_ps(xmm6, xmm4);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_xor_ps(signs24.vector, xmm0);
  xmm1 = _mm_xor_ps(signs24.vector, xmm1);
  xmm2 = _mm_xor_ps(signs24.vector, xmm2);
  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);
  xmm4 = _mm_add_ps(xmm7,xmm4);
  xmm6 = _mm_load_ss( &u[2][2][1] );
  xmm7 = _mm_load_ss( &u[1][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm2, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[0][1][1] );
  xmm7 = _mm_load_ss(&u[2][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm0 = _mm_load_ss( &u[0][2][1] );
  xmm6 = _mm_load_ss( &u[2][1][1] );
  xmm7 = _mm_load_ss( &u[1][2][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm0 = _mm_mul_ps(xmm2, xmm0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm3 = _mm_add_ps(xmm0, xmm3);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);

  /* Result in      xmm3,4,5 */
  /* END MVV */

  /* Reconstruction */

  /* Load up upper partial sum */
  xmm0 = _mm_load_ps(&upper_sum[0][0][0]);
  xmm1 = _mm_load_ps(&upper_sum[1][0][0]);
  xmm2 = _mm_load_ps(&upper_sum[2][0][0]);

  /* Add upper component */
  xmm0 = _mm_add_ps(xmm0, xmm3);
  xmm1 = _mm_add_ps(xmm1, xmm4);
  xmm2 = _mm_add_ps(xmm2, xmm5);
  
  /* Pair store */
  _mm_storel_pi((__m64 *)&spinor_out[0][0][0], xmm0);
  _mm_storel_pi((__m64 *)&spinor_out[0][1][0], xmm1);
  _mm_storel_pi((__m64 *)&spinor_out[0][2][0], xmm2);
  _mm_storeh_pi((__m64 *)&spinor_out[1][0][0], xmm0);
  _mm_storeh_pi((__m64 *)&spinor_out[1][1][0], xmm1);
  _mm_storeh_pi((__m64 *)&spinor_out[1][2][0], xmm2);

  /* Lower components - do projection */
  xmm0 = _mm_load_ps(&lower_sum[0][0][0]);
  xmm1 = _mm_load_ps(&lower_sum[1][0][0]);
  xmm2 = _mm_load_ps(&lower_sum[2][0][0]);
  

  xmm0 = _mm_add_ps(xmm0, xmm3);
  xmm1 = _mm_add_ps(xmm1, xmm4);
  xmm2 = _mm_add_ps(xmm2, xmm5);
  
  /* Store */
  _mm_storel_pi((__m64 *)&spinor_out[2][0][0], xmm0);
  _mm_storel_pi((__m64 *)&spinor_out[2][1][0], xmm1);
  _mm_storel_pi((__m64 *)&spinor_out[2][2][0], xmm2);
  _mm_storeh_pi((__m64 *)&spinor_out[3][0][0], xmm0);
  _mm_storeh_pi((__m64 *)&spinor_out[3][1][0], xmm1);
  _mm_storeh_pi((__m64 *)&spinor_out[3][2][0], xmm2);

  }


  inline
  void dslash_minus_dir0_forward( FourSpinor spinor_in,
				  GaugeMatrix u,
				  HalfSpinor  upper_sum,
				  HalfSpinor  lower_sum   )
       
  {
    SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
    SSESign signs12 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x80000000, 0x00000000, 0x00000000 }};
    SSESign signs14 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x00000000, 0x80000000 }};
    SSESign signs24 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x00000000, 0x80000000 }};
    SSESign signs34 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x00000000, 0x80000000, 0x80000000 }};
    SSESign signs23 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x80000000, 0x00000000 }};

    __m128 xmm0 ALIGN;
    __m128 xmm1 ALIGN;
    __m128 xmm2 ALIGN;
    __m128 xmm3 ALIGN;
    __m128 xmm4 ALIGN;
    __m128 xmm5 ALIGN;
    __m128 xmm6 ALIGN;
    __m128 xmm7 ALIGN;



    /* Component 0 into the low 2 floats */
    xmm0 = _mm_loadl_pi(xmm0, (__m64 *)&spinor_in[0][0][0]);
    xmm1 = _mm_loadl_pi(xmm1, (__m64 *)&spinor_in[0][1][0]);
    xmm2 = _mm_loadl_pi(xmm2, (__m64 *)&spinor_in[0][2][0]);

    /* Component 1 into the high 2 floats */
    xmm0 = _mm_loadh_pi(xmm0, (__m64 *)&spinor_in[1][0][0]);
    xmm1 = _mm_loadh_pi(xmm1, (__m64 *)&spinor_in[1][1][0]);
    xmm2 = _mm_loadh_pi(xmm2, (__m64 *)&spinor_in[1][2][0]);
    
    /* Component 2 into low 2 floats */
    xmm3 = _mm_loadl_pi(xmm3, (__m64 *)&spinor_in[2][0][0]);
    xmm4 = _mm_loadl_pi(xmm4, (__m64 *)&spinor_in[2][1][0]);
    xmm5 = _mm_loadl_pi(xmm5, (__m64 *)&spinor_in[2][2][0]);

    /* Component 3 into the high 2 floats */
    xmm3 = _mm_loadh_pi(xmm3, (__m64 *)&spinor_in[3][0][0]);
    xmm4 = _mm_loadh_pi(xmm4, (__m64 *)&spinor_in[3][1][0]);
    xmm5 = _mm_loadh_pi(xmm5, (__m64 *)&spinor_in[3][2][0]);

    /* Spin Projection. Results into xmm0-xmm2 */
    /* gamma0 minus projection */
    xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x1b);
    xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x1b);
    xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x1b);

    xmm3 = _mm_xor_ps(signs13.vector, xmm3);
    xmm4 = _mm_xor_ps(signs13.vector, xmm4);
    xmm5 = _mm_xor_ps(signs13.vector, xmm5);

    xmm0 = _mm_add_ps(xmm0, xmm3);
    xmm1 = _mm_add_ps(xmm1, xmm4);
    xmm2 = _mm_add_ps(xmm2, xmm5);

    /* SU3 Multiply */
    xmm3 = _mm_load_ss(&u[0][0][0]);
    xmm6 = _mm_load_ss(&u[1][0][0]);
    xmm4 = _mm_load_ss(&u[0][1][0]);
    xmm7 = _mm_load_ss(&u[2][1][0]);
    xmm5 = _mm_load_ss(&u[0][2][0]);
    xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
    xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
    xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
    xmm3 = _mm_mul_ps(xmm0,xmm3);
    xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
    xmm6 = _mm_mul_ps(xmm1,xmm6);
    xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
    xmm4 = _mm_mul_ps(xmm0, xmm4);
    xmm3 = _mm_add_ps(xmm6, xmm3);
    xmm7 = _mm_mul_ps(xmm2, xmm7);
    xmm5 = _mm_mul_ps(xmm0, xmm5);
    xmm4 = _mm_add_ps(xmm7, xmm4);
    xmm6 = _mm_load_ss(&u[1][2][0]);
    xmm7 = _mm_load_ss(&u[2][0][0]);
    xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
    xmm6 = _mm_mul_ps(xmm1, xmm6);
    xmm7 = _mm_mul_ps(xmm2, xmm7);
    xmm5 = _mm_add_ps(xmm6, xmm5);
    xmm3 = _mm_add_ps(xmm7, xmm3);
    xmm6 = _mm_load_ss(&u[1][1][0]);
    xmm7 = _mm_load_ss(&u[2][2][0]);
    xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
    xmm6 = _mm_mul_ps(xmm1, xmm6);
    xmm7 = _mm_mul_ps(xmm2, xmm7);
    xmm4 = _mm_add_ps(xmm6, xmm4);
    xmm5 = _mm_add_ps(xmm7, xmm5);
    xmm6 = _mm_load_ss( &u[0][0][1] );
    xmm7 = _mm_load_ss( &u[1][1][1] );
    xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
    xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
    xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
    xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
    xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
    xmm0 = _mm_xor_ps(signs13.vector, xmm0);
    xmm1 = _mm_xor_ps(signs13.vector, xmm1);
    xmm2 = _mm_xor_ps(signs13.vector, xmm2);
    xmm6 = _mm_mul_ps(xmm0,xmm6);
    xmm7 = _mm_mul_ps(xmm1,xmm7);
    xmm3 = _mm_add_ps(xmm6,xmm3);
    xmm4 = _mm_add_ps(xmm7,xmm4);
    xmm6 = _mm_load_ss( &u[2][2][1] );
    xmm7 = _mm_load_ss( &u[0][1][1] );
    xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
    xmm6 = _mm_mul_ps(xmm2, xmm6);
    xmm7 = _mm_mul_ps(xmm0, xmm7);
    xmm5 = _mm_add_ps(xmm6, xmm5);
    xmm4 = _mm_add_ps(xmm7, xmm4);
    xmm6 = _mm_load_ss(&u[1][0][1] );
    xmm7 = _mm_load_ss(&u[0][2][1] );
    xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
    xmm6 = _mm_mul_ps(xmm1, xmm6);
    xmm7 = _mm_mul_ps(xmm0, xmm7);
    xmm3 = _mm_add_ps(xmm6, xmm3);
    xmm5 = _mm_add_ps(xmm7, xmm5);
    xmm0 = _mm_load_ss( &u[2][0][1] );
    xmm6 = _mm_load_ss( &u[1][2][1] );
    xmm7 = _mm_load_ss( &u[2][1][1] );
    xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0x0);
    xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
    xmm0 = _mm_mul_ps(xmm2, xmm0);
    xmm6 = _mm_mul_ps(xmm1, xmm6);
    xmm7 = _mm_mul_ps(xmm2, xmm7);
    xmm3 = _mm_add_ps(xmm0, xmm3);
    xmm5 = _mm_add_ps(xmm6, xmm5);
    xmm4 = _mm_add_ps(xmm7, xmm4);

    /* Reconstruction: Upper components just go
       Reshuffle spin and color indices so we can use movaps
       to store - it is aligned and faster... */

    _mm_store_ps(&upper_sum[0][0][0], xmm3);
    _mm_store_ps(&upper_sum[1][0][0], xmm4);
    _mm_store_ps(&upper_sum[2][0][0], xmm5);
    
    /* Lower components - do projection */
    xmm3 = _mm_shuffle_ps( xmm3, xmm3, 0x1b);
    xmm4 = _mm_shuffle_ps( xmm4, xmm4, 0x1b);
    xmm5 = _mm_shuffle_ps( xmm5, xmm5, 0x1b);

    xmm3 = _mm_xor_ps(signs24.vector, xmm3);
    xmm4 = _mm_xor_ps(signs24.vector, xmm4);
    xmm5 = _mm_xor_ps(signs24.vector, xmm5);

    /* Store */
    _mm_store_ps(&lower_sum[0][0][0], xmm3);
    _mm_store_ps(&lower_sum[1][0][0], xmm4);
    _mm_store_ps(&lower_sum[2][0][0], xmm5);


  }

  inline
  void dslash_minus_dir0_backward_add( FourSpinor spinor_in,
				       GaugeMatrix u,
				       HalfSpinor  upper_sum,
				       HalfSpinor  lower_sum   )

  {
    SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
    SSESign signs12 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x80000000, 0x00000000, 0x00000000 }};
    SSESign signs14 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x00000000, 0x80000000 }};
    SSESign signs24 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x00000000, 0x80000000 }};
    SSESign signs34 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x00000000, 0x80000000, 0x80000000 }};
    SSESign signs23 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x80000000, 0x00000000 }};

    __m128 xmm0 ALIGN;
    __m128 xmm1 ALIGN;
    __m128 xmm2 ALIGN;
    __m128 xmm3 ALIGN;
    __m128 xmm4 ALIGN;
    __m128 xmm5 ALIGN;
    __m128 xmm6 ALIGN;
    __m128 xmm7 ALIGN;

   

    /* Component 0 into the low 2 floats */
    xmm0 = _mm_loadl_pi(xmm0, (__m64 *)&spinor_in[0][0][0]);
    xmm1 = _mm_loadl_pi(xmm1, (__m64 *)&spinor_in[0][1][0]);
    xmm2 = _mm_loadl_pi(xmm2, (__m64 *)&spinor_in[0][2][0]);

    /* Component 1 into the high 2 floats */
    xmm0 = _mm_loadh_pi(xmm0, (__m64 *)&spinor_in[1][0][0]);
    xmm1 = _mm_loadh_pi(xmm1, (__m64 *)&spinor_in[1][1][0]);
    xmm2 = _mm_loadh_pi(xmm2, (__m64 *)&spinor_in[1][2][0]);
    
    /* Component 2 into low 2 floats */
    xmm3 = _mm_loadl_pi(xmm3, (__m64 *)&spinor_in[2][0][0]);
    xmm4 = _mm_loadl_pi(xmm4, (__m64 *)&spinor_in[2][1][0]);
    xmm5 = _mm_loadl_pi(xmm5, (__m64 *)&spinor_in[2][2][0]);

    /* Component 3 into the high 2 floats */
    xmm3 = _mm_loadh_pi(xmm3, (__m64 *)&spinor_in[3][0][0]);
    xmm4 = _mm_loadh_pi(xmm4, (__m64 *)&spinor_in[3][1][0]);
    xmm5 = _mm_loadh_pi(xmm5, (__m64 *)&spinor_in[3][2][0]);

    /* Spin Projection. Results into xmm0-xmm2 */
    /* gamma0 minus projection */
    xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x1b);
    xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x1b);
    xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x1b);

    xmm3 = _mm_xor_ps(signs24.vector, xmm3);
    xmm4 = _mm_xor_ps(signs24.vector, xmm4);
    xmm5 = _mm_xor_ps(signs24.vector, xmm5);

    xmm0 = _mm_add_ps(xmm0, xmm3);
    xmm1 = _mm_add_ps(xmm1, xmm4);
    xmm2 = _mm_add_ps(xmm2, xmm5);


    /* Adj SU(3) multiply */

  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm6 = _mm_load_ss(&u[0][1][0]);
  xmm4 = _mm_load_ss(&u[1][0][0]);
  xmm7 = _mm_load_ss(&u[1][2][0]);
  xmm5 = _mm_load_ss(&u[2][0][0]);
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
  xmm3 = _mm_mul_ps(xmm0,xmm3);
  xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
  xmm6 = _mm_mul_ps(xmm1,xmm6);
  xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
  xmm4 = _mm_mul_ps(xmm0, xmm4);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_mul_ps(xmm0, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[2][1][0]);
  xmm7 = _mm_load_ss(&u[0][2][0]);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm3 = _mm_add_ps(xmm7, xmm3);
  xmm6 = _mm_load_ss(&u[1][1][0]);
  xmm7 = _mm_load_ss(&u[2][2][0]);
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm4 = _mm_add_ps(xmm6, xmm4);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_xor_ps(signs24.vector, xmm0);
  xmm1 = _mm_xor_ps(signs24.vector, xmm1);
  xmm2 = _mm_xor_ps(signs24.vector, xmm2);
  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);
  xmm4 = _mm_add_ps(xmm7,xmm4);
  xmm6 = _mm_load_ss( &u[2][2][1] );
  xmm7 = _mm_load_ss( &u[1][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm2, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[0][1][1] );
  xmm7 = _mm_load_ss(&u[2][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm0 = _mm_load_ss( &u[0][2][1] );
  xmm6 = _mm_load_ss( &u[2][1][1] );
  xmm7 = _mm_load_ss( &u[1][2][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm0 = _mm_mul_ps(xmm2, xmm0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm3 = _mm_add_ps(xmm0, xmm3);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);

  /* Result in      xmm3,4,5 */
  /* END MVV */

  /* Reconstruction */

  /* Load up upper partial sum */
  xmm0 = _mm_load_ps(&upper_sum[0][0][0]);
  xmm1 = _mm_load_ps(&upper_sum[1][0][0]);
  xmm2 = _mm_load_ps(&upper_sum[2][0][0]);

  /* Add upper component */
  xmm0 = _mm_add_ps(xmm0, xmm3);
  xmm1 = _mm_add_ps(xmm1, xmm4);
  xmm2 = _mm_add_ps(xmm2, xmm5);
  
  _mm_store_ps(&upper_sum[0][0][0], xmm0);
  _mm_store_ps(&upper_sum[1][0][0], xmm1);
  _mm_store_ps(&upper_sum[2][0][0], xmm2);
    
  

  /* Lower components - do projection */
  xmm0 = _mm_load_ps(&lower_sum[0][0][0]);
  xmm1 = _mm_load_ps(&lower_sum[1][0][0]);
  xmm2 = _mm_load_ps(&lower_sum[2][0][0]);
  

  xmm3 = _mm_shuffle_ps( xmm3, xmm3, 0x1b);
  xmm4 = _mm_shuffle_ps( xmm4, xmm4, 0x1b);
  xmm5 = _mm_shuffle_ps( xmm5, xmm5, 0x1b);

  xmm3 = _mm_xor_ps(signs13.vector, xmm3);
  xmm4 = _mm_xor_ps(signs13.vector, xmm4);
  xmm5 = _mm_xor_ps(signs13.vector, xmm5);

  xmm0 = _mm_add_ps(xmm0, xmm3);
  xmm1 = _mm_add_ps(xmm1, xmm4);
  xmm2 = _mm_add_ps(xmm2, xmm5);
  
  /* Store */
  _mm_store_ps(&lower_sum[0][0][0], xmm0);
  _mm_store_ps(&lower_sum[1][0][0], xmm1);
  _mm_store_ps(&lower_sum[2][0][0], xmm2);


  }


  inline
  void dslash_minus_dir1_forward_add( FourSpinor  spinor_in,
				      GaugeMatrix  u,
				      HalfSpinor  upper_sum,
				      HalfSpinor  lower_sum   )
  {
    SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
    SSESign signs12 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x80000000, 0x00000000, 0x00000000 }};
    SSESign signs14 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x00000000, 0x80000000 }};
    SSESign signs24 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x00000000, 0x80000000 }};
    SSESign signs34 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x00000000, 0x80000000, 0x80000000 }};
    SSESign signs23 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x80000000, 0x00000000 }};

    __m128 xmm0 ALIGN;
    __m128 xmm1 ALIGN;
    __m128 xmm2 ALIGN;
    __m128 xmm3 ALIGN;
    __m128 xmm4 ALIGN;
    __m128 xmm5 ALIGN;
    __m128 xmm6 ALIGN;
    __m128 xmm7 ALIGN;
  
    /* Component 0 into the low 2 floats */
    xmm0 = _mm_loadl_pi(xmm0, (__m64 *)&spinor_in[0][0][0]);
    xmm1 = _mm_loadl_pi(xmm1, (__m64 *)&spinor_in[0][1][0]);
    xmm2 = _mm_loadl_pi(xmm2, (__m64 *)&spinor_in[0][2][0]);

    /* Component 1 into the high 2 floats */
    xmm0 = _mm_loadh_pi(xmm0, (__m64 *)&spinor_in[1][0][0]);
    xmm1 = _mm_loadh_pi(xmm1, (__m64 *)&spinor_in[1][1][0]);
    xmm2 = _mm_loadh_pi(xmm2, (__m64 *)&spinor_in[1][2][0]);
    
    /* Component 2 into low 2 floats */
    xmm3 = _mm_loadl_pi(xmm3, (__m64 *)&spinor_in[2][0][0]);
    xmm4 = _mm_loadl_pi(xmm4, (__m64 *)&spinor_in[2][1][0]);
    xmm5 = _mm_loadl_pi(xmm5, (__m64 *)&spinor_in[2][2][0]);

    /* Component 3 into the high 2 floats */
    xmm3 = _mm_loadh_pi(xmm3, (__m64 *)&spinor_in[3][0][0]);
    xmm4 = _mm_loadh_pi(xmm4, (__m64 *)&spinor_in[3][1][0]);
    xmm5 = _mm_loadh_pi(xmm5, (__m64 *)&spinor_in[3][2][0]);

    /* Spin Projection. Results into xmm0-xmm2 */
    /* gamma1 minus projection */
    xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x4e);
    xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x4e);
    xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x4e);

    xmm3 = _mm_xor_ps(signs12.vector, xmm3);
    xmm4 = _mm_xor_ps(signs12.vector, xmm4);
    xmm5 = _mm_xor_ps(signs12.vector, xmm5);

    xmm0 = _mm_add_ps(xmm0, xmm3);
    xmm1 = _mm_add_ps(xmm1, xmm4);
    xmm2 = _mm_add_ps(xmm2, xmm5);


    /* SU3 Multiply */
    xmm3 = _mm_load_ss(&u[0][0][0]);
    xmm6 = _mm_load_ss(&u[1][0][0]);
    xmm4 = _mm_load_ss(&u[0][1][0]);
    xmm7 = _mm_load_ss(&u[2][1][0]);
    xmm5 = _mm_load_ss(&u[0][2][0]);
    xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
    xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
    xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
    xmm3 = _mm_mul_ps(xmm0,xmm3);
    xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
    xmm6 = _mm_mul_ps(xmm1,xmm6);
    xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
    xmm4 = _mm_mul_ps(xmm0, xmm4);
    xmm3 = _mm_add_ps(xmm6, xmm3);
    xmm7 = _mm_mul_ps(xmm2, xmm7);
    xmm5 = _mm_mul_ps(xmm0, xmm5);
    xmm4 = _mm_add_ps(xmm7, xmm4);
    xmm6 = _mm_load_ss(&u[1][2][0]);
    xmm7 = _mm_load_ss(&u[2][0][0]);
    xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
    xmm6 = _mm_mul_ps(xmm1, xmm6);
    xmm7 = _mm_mul_ps(xmm2, xmm7);
    xmm5 = _mm_add_ps(xmm6, xmm5);
    xmm3 = _mm_add_ps(xmm7, xmm3);
    xmm6 = _mm_load_ss(&u[1][1][0]);
    xmm7 = _mm_load_ss(&u[2][2][0]);
    xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
    xmm6 = _mm_mul_ps(xmm1, xmm6);
    xmm7 = _mm_mul_ps(xmm2, xmm7);
    xmm4 = _mm_add_ps(xmm6, xmm4);
    xmm5 = _mm_add_ps(xmm7, xmm5);
    xmm6 = _mm_load_ss( &u[0][0][1] );
    xmm7 = _mm_load_ss( &u[1][1][1] );
    xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
    xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
    xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
    xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
    xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
    xmm0 = _mm_xor_ps(signs13.vector, xmm0);
    xmm1 = _mm_xor_ps(signs13.vector, xmm1);
    xmm2 = _mm_xor_ps(signs13.vector, xmm2);
    xmm6 = _mm_mul_ps(xmm0,xmm6);
    xmm7 = _mm_mul_ps(xmm1,xmm7);
    xmm3 = _mm_add_ps(xmm6,xmm3);
    xmm4 = _mm_add_ps(xmm7,xmm4);
    xmm6 = _mm_load_ss( &u[2][2][1] );
    xmm7 = _mm_load_ss( &u[0][1][1] );
    xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
    xmm6 = _mm_mul_ps(xmm2, xmm6);
    xmm7 = _mm_mul_ps(xmm0, xmm7);
    xmm5 = _mm_add_ps(xmm6, xmm5);
    xmm4 = _mm_add_ps(xmm7, xmm4);
    xmm6 = _mm_load_ss(&u[1][0][1] );
    xmm7 = _mm_load_ss(&u[0][2][1] );
    xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
    xmm6 = _mm_mul_ps(xmm1, xmm6);
    xmm7 = _mm_mul_ps(xmm0, xmm7);
    xmm3 = _mm_add_ps(xmm6, xmm3);
    xmm5 = _mm_add_ps(xmm7, xmm5);
    xmm0 = _mm_load_ss( &u[2][0][1] );
    xmm6 = _mm_load_ss( &u[1][2][1] );
    xmm7 = _mm_load_ss( &u[2][1][1] );
    xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0x0);
    xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
    xmm0 = _mm_mul_ps(xmm2, xmm0);
    xmm6 = _mm_mul_ps(xmm1, xmm6);
    xmm7 = _mm_mul_ps(xmm2, xmm7);
    xmm3 = _mm_add_ps(xmm0, xmm3);
    xmm5 = _mm_add_ps(xmm6, xmm5);
    xmm4 = _mm_add_ps(xmm7, xmm4);

    /* Reconstruction: Upper components just go
       Reshuffle spin and color indices so we can use movaps
       to store - it is aligned and faster... */
    xmm0 = _mm_load_ps(&upper_sum[0][0][0]);
    xmm1 = _mm_load_ps(&upper_sum[1][0][0]);
    xmm2 = _mm_load_ps(&upper_sum[2][0][0]);

    xmm0 = _mm_add_ps(xmm0, xmm3);
    xmm1 = _mm_add_ps(xmm1, xmm4);
    xmm2 = _mm_add_ps(xmm2, xmm5);

    _mm_store_ps(&upper_sum[0][0][0], xmm0);
    _mm_store_ps(&upper_sum[1][0][0], xmm1);
    _mm_store_ps(&upper_sum[2][0][0], xmm2);
    
    /* Lower components - do projection */
    xmm0 = _mm_load_ps(&lower_sum[0][0][0]);
    xmm1 = _mm_load_ps(&lower_sum[1][0][0]);
    xmm2 = _mm_load_ps(&lower_sum[2][0][0]);

    /* Gamma_minus 1 reconstruction */

    xmm3 = _mm_shuffle_ps( xmm3, xmm3, 0x4e);
    xmm4 = _mm_shuffle_ps( xmm4, xmm4, 0x4e);
    xmm5 = _mm_shuffle_ps( xmm5, xmm5, 0x4e);

    xmm3 = _mm_xor_ps(signs34.vector, xmm3);
    xmm4 = _mm_xor_ps(signs34.vector, xmm4);
    xmm5 = _mm_xor_ps(signs34.vector, xmm5);

    /* Accumulate */
    xmm0 = _mm_add_ps(xmm0, xmm3);
    xmm1 = _mm_add_ps(xmm1, xmm4);
    xmm2 = _mm_add_ps(xmm2, xmm5);

    /* Store */
    _mm_store_ps(&lower_sum[0][0][0], xmm0);
    _mm_store_ps(&lower_sum[1][0][0], xmm1);
    _mm_store_ps(&lower_sum[2][0][0], xmm2);

  }

  inline
  void dslash_minus_dir1_backward_add( FourSpinor  spinor_in,
				       GaugeMatrix  u,
				       HalfSpinor  upper_sum,
				       HalfSpinor  lower_sum   )
  {
    SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
    SSESign signs12 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x80000000, 0x00000000, 0x00000000 }};
    SSESign signs14 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x00000000, 0x80000000 }};
    SSESign signs24 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x00000000, 0x80000000 }};
    SSESign signs34 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x00000000, 0x80000000, 0x80000000 }};
    SSESign signs23 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x80000000, 0x00000000 }};

    __m128 xmm0 ALIGN;
    __m128 xmm1 ALIGN;
    __m128 xmm2 ALIGN;
    __m128 xmm3 ALIGN;
    __m128 xmm4 ALIGN;
    __m128 xmm5 ALIGN;
    __m128 xmm6 ALIGN;
    __m128 xmm7 ALIGN;


    /* Component 0 into the low 2 floats */
    xmm0 = _mm_loadl_pi(xmm0, (__m64 *)&spinor_in[0][0][0]);
    xmm1 = _mm_loadl_pi(xmm1, (__m64 *)&spinor_in[0][1][0]);
    xmm2 = _mm_loadl_pi(xmm2, (__m64 *)&spinor_in[0][2][0]);

    /* Component 1 into the high 2 floats */
    xmm0 = _mm_loadh_pi(xmm0, (__m64 *)&spinor_in[1][0][0]);
    xmm1 = _mm_loadh_pi(xmm1, (__m64 *)&spinor_in[1][1][0]);
    xmm2 = _mm_loadh_pi(xmm2, (__m64 *)&spinor_in[1][2][0]);
    
    /* Component 2 into low 2 floats */
    xmm3 = _mm_loadl_pi(xmm3, (__m64 *)&spinor_in[2][0][0]);
    xmm4 = _mm_loadl_pi(xmm4, (__m64 *)&spinor_in[2][1][0]);
    xmm5 = _mm_loadl_pi(xmm5, (__m64 *)&spinor_in[2][2][0]);

    /* Component 3 into the high 2 floats */
    xmm3 = _mm_loadh_pi(xmm3, (__m64 *)&spinor_in[3][0][0]);
    xmm4 = _mm_loadh_pi(xmm4, (__m64 *)&spinor_in[3][1][0]);
    xmm5 = _mm_loadh_pi(xmm5, (__m64 *)&spinor_in[3][2][0]);

    /* Spin Projection. Results into xmm0-xmm2 */
    /* gamma1 minus projection */
    xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x4e);
    xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x4e);
    xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x4e);

    xmm3 = _mm_xor_ps(signs34.vector, xmm3);
    xmm4 = _mm_xor_ps(signs34.vector, xmm4);
    xmm5 = _mm_xor_ps(signs34.vector, xmm5);

    xmm0 = _mm_add_ps(xmm0, xmm3);
    xmm1 = _mm_add_ps(xmm1, xmm4);
    xmm2 = _mm_add_ps(xmm2, xmm5);


    /* Adj SU(3) multiply */

  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm6 = _mm_load_ss(&u[0][1][0]);
  xmm4 = _mm_load_ss(&u[1][0][0]);
  xmm7 = _mm_load_ss(&u[1][2][0]);
  xmm5 = _mm_load_ss(&u[2][0][0]);
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
  xmm3 = _mm_mul_ps(xmm0,xmm3);
  xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
  xmm6 = _mm_mul_ps(xmm1,xmm6);
  xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
  xmm4 = _mm_mul_ps(xmm0, xmm4);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_mul_ps(xmm0, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[2][1][0]);
  xmm7 = _mm_load_ss(&u[0][2][0]);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm3 = _mm_add_ps(xmm7, xmm3);
  xmm6 = _mm_load_ss(&u[1][1][0]);
  xmm7 = _mm_load_ss(&u[2][2][0]);
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm4 = _mm_add_ps(xmm6, xmm4);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_xor_ps(signs24.vector, xmm0);
  xmm1 = _mm_xor_ps(signs24.vector, xmm1);
  xmm2 = _mm_xor_ps(signs24.vector, xmm2);
  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);
  xmm4 = _mm_add_ps(xmm7,xmm4);
  xmm6 = _mm_load_ss( &u[2][2][1] );
  xmm7 = _mm_load_ss( &u[1][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm2, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[0][1][1] );
  xmm7 = _mm_load_ss(&u[2][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm0 = _mm_load_ss( &u[0][2][1] );
  xmm6 = _mm_load_ss( &u[2][1][1] );
  xmm7 = _mm_load_ss( &u[1][2][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm0 = _mm_mul_ps(xmm2, xmm0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm3 = _mm_add_ps(xmm0, xmm3);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);

  /* Result in      xmm3,4,5 */
  /* END MVV */

  /* Reconstruction */

  /* Load up upper partial sum */
  xmm0 = _mm_load_ps(&upper_sum[0][0][0]);
  xmm1 = _mm_load_ps(&upper_sum[1][0][0]);
  xmm2 = _mm_load_ps(&upper_sum[2][0][0]);

  /* Add upper component */
  xmm0 = _mm_add_ps(xmm0, xmm3);
  xmm1 = _mm_add_ps(xmm1, xmm4);
  xmm2 = _mm_add_ps(xmm2, xmm5);
  
  _mm_store_ps(&upper_sum[0][0][0], xmm0);
  _mm_store_ps(&upper_sum[1][0][0], xmm1);
  _mm_store_ps(&upper_sum[2][0][0], xmm2);
    
  

  /* Lower components - do projection */
  xmm0 = _mm_load_ps(&lower_sum[0][0][0]);
  xmm1 = _mm_load_ps(&lower_sum[1][0][0]);
  xmm2 = _mm_load_ps(&lower_sum[2][0][0]);
  

  xmm3 = _mm_shuffle_ps( xmm3, xmm3, 0x4e);
  xmm4 = _mm_shuffle_ps( xmm4, xmm4, 0x4e);
  xmm5 = _mm_shuffle_ps( xmm5, xmm5, 0x4e);

  xmm3 = _mm_xor_ps(signs12.vector, xmm3);
  xmm4 = _mm_xor_ps(signs12.vector, xmm4);
  xmm5 = _mm_xor_ps(signs12.vector, xmm5);

  xmm0 = _mm_add_ps(xmm0, xmm3);
  xmm1 = _mm_add_ps(xmm1, xmm4);
  xmm2 = _mm_add_ps(xmm2, xmm5);
  
  /* Store */
  _mm_store_ps(&lower_sum[0][0][0], xmm0);
  _mm_store_ps(&lower_sum[1][0][0], xmm1);
  _mm_store_ps(&lower_sum[2][0][0], xmm2);


  }

  inline
  void dslash_minus_dir2_forward_add( FourSpinor  spinor_in,
				      GaugeMatrix  u,
				      HalfSpinor  upper_sum,
				      HalfSpinor  lower_sum   )
  {
    SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
    SSESign signs12 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x80000000, 0x00000000, 0x00000000 }};
    SSESign signs14 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x00000000, 0x80000000 }};
    SSESign signs24 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x00000000, 0x80000000 }};
    SSESign signs34 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x00000000, 0x80000000, 0x80000000 }};
    SSESign signs23 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x80000000, 0x00000000 }};

    __m128 xmm0 ALIGN;
    __m128 xmm1 ALIGN;
    __m128 xmm2 ALIGN;
    __m128 xmm3 ALIGN;
    __m128 xmm4 ALIGN;
    __m128 xmm5 ALIGN;
    __m128 xmm6 ALIGN;
    __m128 xmm7 ALIGN;
  
    /* Component 0 into the low 2 floats */
    xmm0 = _mm_loadl_pi(xmm0, (__m64 *)&spinor_in[0][0][0]);
    xmm1 = _mm_loadl_pi(xmm1, (__m64 *)&spinor_in[0][1][0]);
    xmm2 = _mm_loadl_pi(xmm2, (__m64 *)&spinor_in[0][2][0]);

    /* Component 1 into the high 2 floats */
    xmm0 = _mm_loadh_pi(xmm0, (__m64 *)&spinor_in[1][0][0]);
    xmm1 = _mm_loadh_pi(xmm1, (__m64 *)&spinor_in[1][1][0]);
    xmm2 = _mm_loadh_pi(xmm2, (__m64 *)&spinor_in[1][2][0]);
    
    /* Component 2 into low 2 floats */
    xmm3 = _mm_loadl_pi(xmm3, (__m64 *)&spinor_in[2][0][0]);
    xmm4 = _mm_loadl_pi(xmm4, (__m64 *)&spinor_in[2][1][0]);
    xmm5 = _mm_loadl_pi(xmm5, (__m64 *)&spinor_in[2][2][0]);

    /* Component 3 into the high 2 floats */
    xmm3 = _mm_loadh_pi(xmm3, (__m64 *)&spinor_in[3][0][0]);
    xmm4 = _mm_loadh_pi(xmm4, (__m64 *)&spinor_in[3][1][0]);
    xmm5 = _mm_loadh_pi(xmm5, (__m64 *)&spinor_in[3][2][0]);

    /* Spin Projection. Results into xmm0-xmm2 */
    /* gamma1 minus projection */
    xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0xb1);
    xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0xb1);
    xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0xb1);

    xmm3 = _mm_xor_ps(signs14.vector, xmm3);
    xmm4 = _mm_xor_ps(signs14.vector, xmm4);
    xmm5 = _mm_xor_ps(signs14.vector, xmm5);

    xmm0 = _mm_add_ps(xmm0, xmm3);
    xmm1 = _mm_add_ps(xmm1, xmm4);
    xmm2 = _mm_add_ps(xmm2, xmm5);


    /* SU3 Multiply */
    xmm3 = _mm_load_ss(&u[0][0][0]);
    xmm6 = _mm_load_ss(&u[1][0][0]);
    xmm4 = _mm_load_ss(&u[0][1][0]);
    xmm7 = _mm_load_ss(&u[2][1][0]);
    xmm5 = _mm_load_ss(&u[0][2][0]);
    xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
    xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
    xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
    xmm3 = _mm_mul_ps(xmm0,xmm3);
    xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
    xmm6 = _mm_mul_ps(xmm1,xmm6);
    xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
    xmm4 = _mm_mul_ps(xmm0, xmm4);
    xmm3 = _mm_add_ps(xmm6, xmm3);
    xmm7 = _mm_mul_ps(xmm2, xmm7);
    xmm5 = _mm_mul_ps(xmm0, xmm5);
    xmm4 = _mm_add_ps(xmm7, xmm4);
    xmm6 = _mm_load_ss(&u[1][2][0]);
    xmm7 = _mm_load_ss(&u[2][0][0]);
    xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
    xmm6 = _mm_mul_ps(xmm1, xmm6);
    xmm7 = _mm_mul_ps(xmm2, xmm7);
    xmm5 = _mm_add_ps(xmm6, xmm5);
    xmm3 = _mm_add_ps(xmm7, xmm3);
    xmm6 = _mm_load_ss(&u[1][1][0]);
    xmm7 = _mm_load_ss(&u[2][2][0]);
    xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
    xmm6 = _mm_mul_ps(xmm1, xmm6);
    xmm7 = _mm_mul_ps(xmm2, xmm7);
    xmm4 = _mm_add_ps(xmm6, xmm4);
    xmm5 = _mm_add_ps(xmm7, xmm5);
    xmm6 = _mm_load_ss( &u[0][0][1] );
    xmm7 = _mm_load_ss( &u[1][1][1] );
    xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
    xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
    xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
    xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
    xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
    xmm0 = _mm_xor_ps(signs13.vector, xmm0);
    xmm1 = _mm_xor_ps(signs13.vector, xmm1);
    xmm2 = _mm_xor_ps(signs13.vector, xmm2);
    xmm6 = _mm_mul_ps(xmm0,xmm6);
    xmm7 = _mm_mul_ps(xmm1,xmm7);
    xmm3 = _mm_add_ps(xmm6,xmm3);
    xmm4 = _mm_add_ps(xmm7,xmm4);
    xmm6 = _mm_load_ss( &u[2][2][1] );
    xmm7 = _mm_load_ss( &u[0][1][1] );
    xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
    xmm6 = _mm_mul_ps(xmm2, xmm6);
    xmm7 = _mm_mul_ps(xmm0, xmm7);
    xmm5 = _mm_add_ps(xmm6, xmm5);
    xmm4 = _mm_add_ps(xmm7, xmm4);
    xmm6 = _mm_load_ss(&u[1][0][1] );
    xmm7 = _mm_load_ss(&u[0][2][1] );
    xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
    xmm6 = _mm_mul_ps(xmm1, xmm6);
    xmm7 = _mm_mul_ps(xmm0, xmm7);
    xmm3 = _mm_add_ps(xmm6, xmm3);
    xmm5 = _mm_add_ps(xmm7, xmm5);
    xmm0 = _mm_load_ss( &u[2][0][1] );
    xmm6 = _mm_load_ss( &u[1][2][1] );
    xmm7 = _mm_load_ss( &u[2][1][1] );
    xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0x0);
    xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
    xmm0 = _mm_mul_ps(xmm2, xmm0);
    xmm6 = _mm_mul_ps(xmm1, xmm6);
    xmm7 = _mm_mul_ps(xmm2, xmm7);
    xmm3 = _mm_add_ps(xmm0, xmm3);
    xmm5 = _mm_add_ps(xmm6, xmm5);
    xmm4 = _mm_add_ps(xmm7, xmm4);

    /* Reconstruction: Upper components just go
       Reshuffle spin and color indices so we can use movaps
       to store - it is aligned and faster... */
    xmm0 = _mm_load_ps(&upper_sum[0][0][0]);
    xmm1 = _mm_load_ps(&upper_sum[1][0][0]);
    xmm2 = _mm_load_ps(&upper_sum[2][0][0]);

    xmm0 = _mm_add_ps(xmm0, xmm3);
    xmm1 = _mm_add_ps(xmm1, xmm4);
    xmm2 = _mm_add_ps(xmm2, xmm5);

    _mm_store_ps(&upper_sum[0][0][0], xmm0);
    _mm_store_ps(&upper_sum[1][0][0], xmm1);
    _mm_store_ps(&upper_sum[2][0][0], xmm2);
    
    /* Lower components - do projection */
    xmm0 = _mm_load_ps(&lower_sum[0][0][0]);
    xmm1 = _mm_load_ps(&lower_sum[1][0][0]);
    xmm2 = _mm_load_ps(&lower_sum[2][0][0]);

    /* Gamma_minus 1 reconstruction */

    xmm3 = _mm_shuffle_ps( xmm3, xmm3, 0xb1);
    xmm4 = _mm_shuffle_ps( xmm4, xmm4, 0xb1);
    xmm5 = _mm_shuffle_ps( xmm5, xmm5, 0xb1);

    xmm3 = _mm_xor_ps(signs23.vector, xmm3);
    xmm4 = _mm_xor_ps(signs23.vector, xmm4);
    xmm5 = _mm_xor_ps(signs23.vector, xmm5);

    /* Accumulate */
    xmm0 = _mm_add_ps(xmm0, xmm3);
    xmm1 = _mm_add_ps(xmm1, xmm4);
    xmm2 = _mm_add_ps(xmm2, xmm5);

    /* Store */
    _mm_store_ps(&lower_sum[0][0][0], xmm0);
    _mm_store_ps(&lower_sum[1][0][0], xmm1);
    _mm_store_ps(&lower_sum[2][0][0], xmm2);


  }

  inline
  void dslash_minus_dir2_backward_add( FourSpinor  spinor_in,
				       GaugeMatrix  u,
				       HalfSpinor  upper_sum,
				       HalfSpinor  lower_sum   )
  {
    SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
    SSESign signs12 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x80000000, 0x00000000, 0x00000000 }};
    SSESign signs14 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x00000000, 0x80000000 }};
    SSESign signs24 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x00000000, 0x80000000 }};
    SSESign signs34 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x00000000, 0x80000000, 0x80000000 }};
    SSESign signs23 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x80000000, 0x00000000 }};

    __m128 xmm0 ALIGN;
    __m128 xmm1 ALIGN;
    __m128 xmm2 ALIGN;
    __m128 xmm3 ALIGN;
    __m128 xmm4 ALIGN;
    __m128 xmm5 ALIGN;
    __m128 xmm6 ALIGN;
    __m128 xmm7 ALIGN;


    /* Component 0 into the low 2 floats */
    xmm0 = _mm_loadl_pi(xmm0, (__m64 *)&spinor_in[0][0][0]);
    xmm1 = _mm_loadl_pi(xmm1, (__m64 *)&spinor_in[0][1][0]);
    xmm2 = _mm_loadl_pi(xmm2, (__m64 *)&spinor_in[0][2][0]);

    /* Component 1 into the high 2 floats */
    xmm0 = _mm_loadh_pi(xmm0, (__m64 *)&spinor_in[1][0][0]);
    xmm1 = _mm_loadh_pi(xmm1, (__m64 *)&spinor_in[1][1][0]);
    xmm2 = _mm_loadh_pi(xmm2, (__m64 *)&spinor_in[1][2][0]);
    
    /* Component 2 into low 2 floats */
    xmm3 = _mm_loadl_pi(xmm3, (__m64 *)&spinor_in[2][0][0]);
    xmm4 = _mm_loadl_pi(xmm4, (__m64 *)&spinor_in[2][1][0]);
    xmm5 = _mm_loadl_pi(xmm5, (__m64 *)&spinor_in[2][2][0]);

    /* Component 3 into the high 2 floats */
    xmm3 = _mm_loadh_pi(xmm3, (__m64 *)&spinor_in[3][0][0]);
    xmm4 = _mm_loadh_pi(xmm4, (__m64 *)&spinor_in[3][1][0]);
    xmm5 = _mm_loadh_pi(xmm5, (__m64 *)&spinor_in[3][2][0]);

    /* Spin Projection. Results into xmm0-xmm2 */
    /* gamma1 minus projection */
    xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0xb1);
    xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0xb1);
    xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0xb1);

    xmm3 = _mm_xor_ps(signs23.vector, xmm3);
    xmm4 = _mm_xor_ps(signs23.vector, xmm4);
    xmm5 = _mm_xor_ps(signs23.vector, xmm5);

    xmm0 = _mm_add_ps(xmm0, xmm3);
    xmm1 = _mm_add_ps(xmm1, xmm4);
    xmm2 = _mm_add_ps(xmm2, xmm5);


    /* Adj SU(3) multiply */

  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm6 = _mm_load_ss(&u[0][1][0]);
  xmm4 = _mm_load_ss(&u[1][0][0]);
  xmm7 = _mm_load_ss(&u[1][2][0]);
  xmm5 = _mm_load_ss(&u[2][0][0]);
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
  xmm3 = _mm_mul_ps(xmm0,xmm3);
  xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
  xmm6 = _mm_mul_ps(xmm1,xmm6);
  xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
  xmm4 = _mm_mul_ps(xmm0, xmm4);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_mul_ps(xmm0, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[2][1][0]);
  xmm7 = _mm_load_ss(&u[0][2][0]);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm3 = _mm_add_ps(xmm7, xmm3);
  xmm6 = _mm_load_ss(&u[1][1][0]);
  xmm7 = _mm_load_ss(&u[2][2][0]);
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm4 = _mm_add_ps(xmm6, xmm4);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_xor_ps(signs24.vector, xmm0);
  xmm1 = _mm_xor_ps(signs24.vector, xmm1);
  xmm2 = _mm_xor_ps(signs24.vector, xmm2);
  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);
  xmm4 = _mm_add_ps(xmm7,xmm4);
  xmm6 = _mm_load_ss( &u[2][2][1] );
  xmm7 = _mm_load_ss( &u[1][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm2, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[0][1][1] );
  xmm7 = _mm_load_ss(&u[2][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm0 = _mm_load_ss( &u[0][2][1] );
  xmm6 = _mm_load_ss( &u[2][1][1] );
  xmm7 = _mm_load_ss( &u[1][2][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm0 = _mm_mul_ps(xmm2, xmm0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm3 = _mm_add_ps(xmm0, xmm3);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);

  /* Result in      xmm3,4,5 */
  /* END MVV */

  /* Reconstruction */

  /* Load up upper partial sum */
  xmm0 = _mm_load_ps(&upper_sum[0][0][0]);
  xmm1 = _mm_load_ps(&upper_sum[1][0][0]);
  xmm2 = _mm_load_ps(&upper_sum[2][0][0]);

  /* Add upper component */
  xmm0 = _mm_add_ps(xmm0, xmm3);
  xmm1 = _mm_add_ps(xmm1, xmm4);
  xmm2 = _mm_add_ps(xmm2, xmm5);
  
  _mm_store_ps(&upper_sum[0][0][0], xmm0);
  _mm_store_ps(&upper_sum[1][0][0], xmm1);
  _mm_store_ps(&upper_sum[2][0][0], xmm2);
    
  

  /* Lower components - do projection */
  xmm0 = _mm_load_ps(&lower_sum[0][0][0]);
  xmm1 = _mm_load_ps(&lower_sum[1][0][0]);
  xmm2 = _mm_load_ps(&lower_sum[2][0][0]);
  

  xmm3 = _mm_shuffle_ps( xmm3, xmm3, 0xb1);
  xmm4 = _mm_shuffle_ps( xmm4, xmm4, 0xb1);
  xmm5 = _mm_shuffle_ps( xmm5, xmm5, 0xb1);

  xmm3 = _mm_xor_ps(signs14.vector, xmm3);
  xmm4 = _mm_xor_ps(signs14.vector, xmm4);
  xmm5 = _mm_xor_ps(signs14.vector, xmm5);

  xmm0 = _mm_add_ps(xmm0, xmm3);
  xmm1 = _mm_add_ps(xmm1, xmm4);
  xmm2 = _mm_add_ps(xmm2, xmm5);
  
  /* Store */
  _mm_store_ps(&lower_sum[0][0][0], xmm0);
  _mm_store_ps(&lower_sum[1][0][0], xmm1);
  _mm_store_ps(&lower_sum[2][0][0], xmm2);


  }

  inline
  void dslash_minus_dir2_backward_add_store( FourSpinor  spinor_in,
					     GaugeMatrix  u,
					     HalfSpinor  upper_sum,
					     HalfSpinor  lower_sum,
					     FourSpinor spinor_out)
  {
    SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
    SSESign signs12 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x80000000, 0x00000000, 0x00000000 }};
    SSESign signs14 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x00000000, 0x80000000 }};
    SSESign signs24 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x00000000, 0x80000000 }};
    SSESign signs34 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x00000000, 0x80000000, 0x80000000 }};
    SSESign signs23 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x80000000, 0x00000000 }};

    __m128 xmm0 ALIGN;
    __m128 xmm1 ALIGN;
    __m128 xmm2 ALIGN;
    __m128 xmm3 ALIGN;
    __m128 xmm4 ALIGN;
    __m128 xmm5 ALIGN;
    __m128 xmm6 ALIGN;
    __m128 xmm7 ALIGN;


    /* Component 0 into the low 2 floats */
    xmm0 = _mm_loadl_pi(xmm0, (__m64 *)&spinor_in[0][0][0]);
    xmm1 = _mm_loadl_pi(xmm1, (__m64 *)&spinor_in[0][1][0]);
    xmm2 = _mm_loadl_pi(xmm2, (__m64 *)&spinor_in[0][2][0]);

    /* Component 1 into the high 2 floats */
    xmm0 = _mm_loadh_pi(xmm0, (__m64 *)&spinor_in[1][0][0]);
    xmm1 = _mm_loadh_pi(xmm1, (__m64 *)&spinor_in[1][1][0]);
    xmm2 = _mm_loadh_pi(xmm2, (__m64 *)&spinor_in[1][2][0]);
    
    /* Component 2 into low 2 floats */
    xmm3 = _mm_loadl_pi(xmm3, (__m64 *)&spinor_in[2][0][0]);
    xmm4 = _mm_loadl_pi(xmm4, (__m64 *)&spinor_in[2][1][0]);
    xmm5 = _mm_loadl_pi(xmm5, (__m64 *)&spinor_in[2][2][0]);

    /* Component 3 into the high 2 floats */
    xmm3 = _mm_loadh_pi(xmm3, (__m64 *)&spinor_in[3][0][0]);
    xmm4 = _mm_loadh_pi(xmm4, (__m64 *)&spinor_in[3][1][0]);
    xmm5 = _mm_loadh_pi(xmm5, (__m64 *)&spinor_in[3][2][0]);

    /* Spin Projection. Results into xmm0-xmm2 */
    /* gamma1 minus projection */
    xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0xb1);
    xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0xb1);
    xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0xb1);

    xmm3 = _mm_xor_ps(signs23.vector, xmm3);
    xmm4 = _mm_xor_ps(signs23.vector, xmm4);
    xmm5 = _mm_xor_ps(signs23.vector, xmm5);

    xmm0 = _mm_add_ps(xmm0, xmm3);
    xmm1 = _mm_add_ps(xmm1, xmm4);
    xmm2 = _mm_add_ps(xmm2, xmm5);


    /* Adj SU(3) multiply */

  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm6 = _mm_load_ss(&u[0][1][0]);
  xmm4 = _mm_load_ss(&u[1][0][0]);
  xmm7 = _mm_load_ss(&u[1][2][0]);
  xmm5 = _mm_load_ss(&u[2][0][0]);
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
  xmm3 = _mm_mul_ps(xmm0,xmm3);
  xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
  xmm6 = _mm_mul_ps(xmm1,xmm6);
  xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
  xmm4 = _mm_mul_ps(xmm0, xmm4);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_mul_ps(xmm0, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[2][1][0]);
  xmm7 = _mm_load_ss(&u[0][2][0]);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm3 = _mm_add_ps(xmm7, xmm3);
  xmm6 = _mm_load_ss(&u[1][1][0]);
  xmm7 = _mm_load_ss(&u[2][2][0]);
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm4 = _mm_add_ps(xmm6, xmm4);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_xor_ps(signs24.vector, xmm0);
  xmm1 = _mm_xor_ps(signs24.vector, xmm1);
  xmm2 = _mm_xor_ps(signs24.vector, xmm2);
  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);
  xmm4 = _mm_add_ps(xmm7,xmm4);
  xmm6 = _mm_load_ss( &u[2][2][1] );
  xmm7 = _mm_load_ss( &u[1][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm2, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[0][1][1] );
  xmm7 = _mm_load_ss(&u[2][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm0 = _mm_load_ss( &u[0][2][1] );
  xmm6 = _mm_load_ss( &u[2][1][1] );
  xmm7 = _mm_load_ss( &u[1][2][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm0 = _mm_mul_ps(xmm2, xmm0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm3 = _mm_add_ps(xmm0, xmm3);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);

  /* Result in      xmm3,4,5 */
  /* END MVV */

  /* Reconstruction */

  /* Load up upper partial sum */
  xmm0 = _mm_load_ps(&upper_sum[0][0][0]);
  xmm1 = _mm_load_ps(&upper_sum[1][0][0]);
  xmm2 = _mm_load_ps(&upper_sum[2][0][0]);

  /* Add upper component */
  xmm0 = _mm_add_ps(xmm0, xmm3);
  xmm1 = _mm_add_ps(xmm1, xmm4);
  xmm2 = _mm_add_ps(xmm2, xmm5);
  
  /* Was */
  /*
  _mm_store_ps(&upper_sum[0][0][0], xmm0);
  _mm_store_ps(&upper_sum[1][0][0], xmm1);
  _mm_store_ps(&upper_sum[2][0][0], xmm2);
  */

  /* Pair store */
  _mm_storel_pi((__m64 *)&spinor_out[0][0][0], xmm0);
  _mm_storel_pi((__m64 *)&spinor_out[0][1][0], xmm1);
  _mm_storel_pi((__m64 *)&spinor_out[0][2][0], xmm2);
  _mm_storeh_pi((__m64 *)&spinor_out[1][0][0], xmm0);
  _mm_storeh_pi((__m64 *)&spinor_out[1][1][0], xmm1);
  _mm_storeh_pi((__m64 *)&spinor_out[1][2][0], xmm2);


  /* Lower components - do projection */
  xmm0 = _mm_load_ps(&lower_sum[0][0][0]);
  xmm1 = _mm_load_ps(&lower_sum[1][0][0]);
  xmm2 = _mm_load_ps(&lower_sum[2][0][0]);
  

  xmm3 = _mm_shuffle_ps( xmm3, xmm3, 0xb1);
  xmm4 = _mm_shuffle_ps( xmm4, xmm4, 0xb1);
  xmm5 = _mm_shuffle_ps( xmm5, xmm5, 0xb1);

  xmm3 = _mm_xor_ps(signs14.vector, xmm3);
  xmm4 = _mm_xor_ps(signs14.vector, xmm4);
  xmm5 = _mm_xor_ps(signs14.vector, xmm5);

  xmm0 = _mm_add_ps(xmm0, xmm3);
  xmm1 = _mm_add_ps(xmm1, xmm4);
  xmm2 = _mm_add_ps(xmm2, xmm5);
  
  /* Store */
  /* Was */
  /*
  _mm_store_ps(&lower_sum[0][0][0], xmm0);
  _mm_store_ps(&lower_sum[1][0][0], xmm1);
  _mm_store_ps(&lower_sum[2][0][0], xmm2);
  */

  _mm_storel_pi((__m64 *)&spinor_out[2][0][0], xmm0);
  _mm_storel_pi((__m64 *)&spinor_out[2][1][0], xmm1);
  _mm_storel_pi((__m64 *)&spinor_out[2][2][0], xmm2);
  _mm_storeh_pi((__m64 *)&spinor_out[3][0][0], xmm0);
  _mm_storeh_pi((__m64 *)&spinor_out[3][1][0], xmm1);
  _mm_storeh_pi((__m64 *)&spinor_out[3][2][0], xmm2);

  }

  inline
  void dslash_minus_dir3_forward_add( FourSpinor  spinor_in,
				      GaugeMatrix  u,
				      HalfSpinor  upper_sum,
				      HalfSpinor  lower_sum)
  {
    SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
    SSESign signs12 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x80000000, 0x00000000, 0x00000000 }};
    SSESign signs14 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x00000000, 0x80000000 }};
    SSESign signs24 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x00000000, 0x80000000 }};
    SSESign signs34 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x00000000, 0x80000000, 0x80000000 }};
    SSESign signs23 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x80000000, 0x00000000 }};

   __m128 xmm0 ALIGN;
    __m128 xmm1 ALIGN;
    __m128 xmm2 ALIGN;
    __m128 xmm3 ALIGN;
    __m128 xmm4 ALIGN;
    __m128 xmm5 ALIGN;
    __m128 xmm6 ALIGN;
    __m128 xmm7 ALIGN;
  
    /* Component 0 into the low 2 floats */
    xmm0 = _mm_loadl_pi(xmm0, (__m64 *)&spinor_in[0][0][0]);
    xmm1 = _mm_loadl_pi(xmm1, (__m64 *)&spinor_in[0][1][0]);
    xmm2 = _mm_loadl_pi(xmm2, (__m64 *)&spinor_in[0][2][0]);

    /* Component 1 into the high 2 floats */
    xmm0 = _mm_loadh_pi(xmm0, (__m64 *)&spinor_in[1][0][0]);
    xmm1 = _mm_loadh_pi(xmm1, (__m64 *)&spinor_in[1][1][0]);
    xmm2 = _mm_loadh_pi(xmm2, (__m64 *)&spinor_in[1][2][0]);
    
    /* Component 2 into low 2 floats */
    xmm3 = _mm_loadl_pi(xmm3, (__m64 *)&spinor_in[2][0][0]);
    xmm4 = _mm_loadl_pi(xmm4, (__m64 *)&spinor_in[2][1][0]);
    xmm5 = _mm_loadl_pi(xmm5, (__m64 *)&spinor_in[2][2][0]);

    /* Component 3 into the high 2 floats */
    xmm3 = _mm_loadh_pi(xmm3, (__m64 *)&spinor_in[3][0][0]);
    xmm4 = _mm_loadh_pi(xmm4, (__m64 *)&spinor_in[3][1][0]);
    xmm5 = _mm_loadh_pi(xmm5, (__m64 *)&spinor_in[3][2][0]);

    /* Spin Projection. Results into xmm0-xmm2 */
    xmm0 = _mm_add_ps(xmm0, xmm3);
    xmm1 = _mm_add_ps(xmm1, xmm4);
    xmm2 = _mm_add_ps(xmm2, xmm5);


    /* SU3 Multiply */
    xmm3 = _mm_load_ss(&u[0][0][0]);
    xmm6 = _mm_load_ss(&u[1][0][0]);
    xmm4 = _mm_load_ss(&u[0][1][0]);
    xmm7 = _mm_load_ss(&u[2][1][0]);
    xmm5 = _mm_load_ss(&u[0][2][0]);
    xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
    xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
    xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
    xmm3 = _mm_mul_ps(xmm0,xmm3);
    xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
    xmm6 = _mm_mul_ps(xmm1,xmm6);
    xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
    xmm4 = _mm_mul_ps(xmm0, xmm4);
    xmm3 = _mm_add_ps(xmm6, xmm3);
    xmm7 = _mm_mul_ps(xmm2, xmm7);
    xmm5 = _mm_mul_ps(xmm0, xmm5);
    xmm4 = _mm_add_ps(xmm7, xmm4);
    xmm6 = _mm_load_ss(&u[1][2][0]);
    xmm7 = _mm_load_ss(&u[2][0][0]);
    xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
    xmm6 = _mm_mul_ps(xmm1, xmm6);
    xmm7 = _mm_mul_ps(xmm2, xmm7);
    xmm5 = _mm_add_ps(xmm6, xmm5);
    xmm3 = _mm_add_ps(xmm7, xmm3);
    xmm6 = _mm_load_ss(&u[1][1][0]);
    xmm7 = _mm_load_ss(&u[2][2][0]);
    xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
    xmm6 = _mm_mul_ps(xmm1, xmm6);
    xmm7 = _mm_mul_ps(xmm2, xmm7);
    xmm4 = _mm_add_ps(xmm6, xmm4);
    xmm5 = _mm_add_ps(xmm7, xmm5);
    xmm6 = _mm_load_ss( &u[0][0][1] );
    xmm7 = _mm_load_ss( &u[1][1][1] );
    xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
    xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
    xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
    xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
    xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
    xmm0 = _mm_xor_ps(signs13.vector, xmm0);
    xmm1 = _mm_xor_ps(signs13.vector, xmm1);
    xmm2 = _mm_xor_ps(signs13.vector, xmm2);
    xmm6 = _mm_mul_ps(xmm0,xmm6);
    xmm7 = _mm_mul_ps(xmm1,xmm7);
    xmm3 = _mm_add_ps(xmm6,xmm3);
    xmm4 = _mm_add_ps(xmm7,xmm4);
    xmm6 = _mm_load_ss( &u[2][2][1] );
    xmm7 = _mm_load_ss( &u[0][1][1] );
    xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
    xmm6 = _mm_mul_ps(xmm2, xmm6);
    xmm7 = _mm_mul_ps(xmm0, xmm7);
    xmm5 = _mm_add_ps(xmm6, xmm5);
    xmm4 = _mm_add_ps(xmm7, xmm4);
    xmm6 = _mm_load_ss(&u[1][0][1] );
    xmm7 = _mm_load_ss(&u[0][2][1] );
    xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
    xmm6 = _mm_mul_ps(xmm1, xmm6);
    xmm7 = _mm_mul_ps(xmm0, xmm7);
    xmm3 = _mm_add_ps(xmm6, xmm3);
    xmm5 = _mm_add_ps(xmm7, xmm5);
    xmm0 = _mm_load_ss( &u[2][0][1] );
    xmm6 = _mm_load_ss( &u[1][2][1] );
    xmm7 = _mm_load_ss( &u[2][1][1] );
    xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0x0);
    xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
    xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
    xmm0 = _mm_mul_ps(xmm2, xmm0);
    xmm6 = _mm_mul_ps(xmm1, xmm6);
    xmm7 = _mm_mul_ps(xmm2, xmm7);
    xmm3 = _mm_add_ps(xmm0, xmm3);
    xmm5 = _mm_add_ps(xmm6, xmm5);
    xmm4 = _mm_add_ps(xmm7, xmm4);

    /* Reconstruction: Upper components just go
       Reshuffle spin and color indices so we can use movaps
       to store - it is aligned and faster... */
    xmm0 = _mm_load_ps(&upper_sum[0][0][0]);
    xmm1 = _mm_load_ps(&upper_sum[1][0][0]);
    xmm2 = _mm_load_ps(&upper_sum[2][0][0]);

    xmm0 = _mm_add_ps(xmm0, xmm3);
    xmm1 = _mm_add_ps(xmm1, xmm4);
    xmm2 = _mm_add_ps(xmm2, xmm5);

    _mm_store_ps(&upper_sum[0][0][0], xmm0);
    _mm_store_ps(&upper_sum[1][0][0], xmm1);
    _mm_store_ps(&upper_sum[2][0][0], xmm2);
    
    /* Lower components - do projection */
    xmm0 = _mm_load_ps(&lower_sum[0][0][0]);
    xmm1 = _mm_load_ps(&lower_sum[1][0][0]);
    xmm2 = _mm_load_ps(&lower_sum[2][0][0]);

    /* Accumulate */
    xmm0 = _mm_add_ps(xmm0, xmm3);
    xmm1 = _mm_add_ps(xmm1, xmm4);
    xmm2 = _mm_add_ps(xmm2, xmm5);

    /* Store */
    _mm_store_ps(&lower_sum[0][0][0], xmm0);
    _mm_store_ps(&lower_sum[1][0][0], xmm1);
    _mm_store_ps(&lower_sum[2][0][0], xmm2);



  }

  inline
  void dslash_minus_dir3_backward_add_store( FourSpinor  spinor_in,
					     GaugeMatrix  u,
					     HalfSpinor  upper_sum,
					     HalfSpinor  lower_sum,
					     FourSpinor spinor_out)
       
  {
    SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
    SSESign signs12 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x80000000, 0x00000000, 0x00000000 }};
    SSESign signs14 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x00000000, 0x80000000 }};
    SSESign signs24 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x00000000, 0x80000000 }};
    SSESign signs34 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x00000000, 0x80000000, 0x80000000 }};
    SSESign signs23 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x80000000, 0x00000000 }};

    __m128 xmm0 ALIGN;
    __m128 xmm1 ALIGN;
    __m128 xmm2 ALIGN;
    __m128 xmm3 ALIGN;
    __m128 xmm4 ALIGN;
    __m128 xmm5 ALIGN;
    __m128 xmm6 ALIGN;
    __m128 xmm7 ALIGN;



    /* Component 0 into the low 2 floats */
    xmm0 = _mm_loadl_pi(xmm0, (__m64 *)&spinor_in[0][0][0]);
    xmm1 = _mm_loadl_pi(xmm1, (__m64 *)&spinor_in[0][1][0]);
    xmm2 = _mm_loadl_pi(xmm2, (__m64 *)&spinor_in[0][2][0]);

    /* Component 1 into the high 2 floats */
    xmm0 = _mm_loadh_pi(xmm0, (__m64 *)&spinor_in[1][0][0]);
    xmm1 = _mm_loadh_pi(xmm1, (__m64 *)&spinor_in[1][1][0]);
    xmm2 = _mm_loadh_pi(xmm2, (__m64 *)&spinor_in[1][2][0]);
    
    /* Component 2 into low 2 floats */
    xmm3 = _mm_loadl_pi(xmm3, (__m64 *)&spinor_in[2][0][0]);
    xmm4 = _mm_loadl_pi(xmm4, (__m64 *)&spinor_in[2][1][0]);
    xmm5 = _mm_loadl_pi(xmm5, (__m64 *)&spinor_in[2][2][0]);

    /* Component 3 into the high 2 floats */
    xmm3 = _mm_loadh_pi(xmm3, (__m64 *)&spinor_in[3][0][0]);
    xmm4 = _mm_loadh_pi(xmm4, (__m64 *)&spinor_in[3][1][0]);
    xmm5 = _mm_loadh_pi(xmm5, (__m64 *)&spinor_in[3][2][0]);

    /* Spin Projection. Results into xmm0-xmm2 */
    xmm0 = _mm_sub_ps(xmm0, xmm3);
    xmm1 = _mm_sub_ps(xmm1, xmm4);
    xmm2 = _mm_sub_ps(xmm2, xmm5);


    /* Adj SU(3) multiply */

  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm6 = _mm_load_ss(&u[0][1][0]);
  xmm4 = _mm_load_ss(&u[1][0][0]);
  xmm7 = _mm_load_ss(&u[1][2][0]);
  xmm5 = _mm_load_ss(&u[2][0][0]);
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
  xmm3 = _mm_mul_ps(xmm0,xmm3);
  xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
  xmm6 = _mm_mul_ps(xmm1,xmm6);
  xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
  xmm4 = _mm_mul_ps(xmm0, xmm4);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_mul_ps(xmm0, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[2][1][0]);
  xmm7 = _mm_load_ss(&u[0][2][0]);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm3 = _mm_add_ps(xmm7, xmm3);
  xmm6 = _mm_load_ss(&u[1][1][0]);
  xmm7 = _mm_load_ss(&u[2][2][0]);
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm4 = _mm_add_ps(xmm6, xmm4);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_xor_ps(signs24.vector, xmm0);
  xmm1 = _mm_xor_ps(signs24.vector, xmm1);
  xmm2 = _mm_xor_ps(signs24.vector, xmm2);
  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);
  xmm4 = _mm_add_ps(xmm7,xmm4);
  xmm6 = _mm_load_ss( &u[2][2][1] );
  xmm7 = _mm_load_ss( &u[1][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm2, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[0][1][1] );
  xmm7 = _mm_load_ss(&u[2][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm0 = _mm_load_ss( &u[0][2][1] );
  xmm6 = _mm_load_ss( &u[2][1][1] );
  xmm7 = _mm_load_ss( &u[1][2][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm0 = _mm_mul_ps(xmm2, xmm0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm3 = _mm_add_ps(xmm0, xmm3);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);

  /* Result in      xmm3,4,5 */
  /* END MVV */

  /* Reconstruction */

  /* Load up upper partial sum */
  xmm0 = _mm_load_ps(&upper_sum[0][0][0]);
  xmm1 = _mm_load_ps(&upper_sum[1][0][0]);
  xmm2 = _mm_load_ps(&upper_sum[2][0][0]);

  /* Add upper component */
  xmm0 = _mm_add_ps(xmm0, xmm3);
  xmm1 = _mm_add_ps(xmm1, xmm4);
  xmm2 = _mm_add_ps(xmm2, xmm5);
  
  /* Pair store */
  _mm_storel_pi((__m64 *)&spinor_out[0][0][0], xmm0);
  _mm_storel_pi((__m64 *)&spinor_out[0][1][0], xmm1);
  _mm_storel_pi((__m64 *)&spinor_out[0][2][0], xmm2);
  _mm_storeh_pi((__m64 *)&spinor_out[1][0][0], xmm0);
  _mm_storeh_pi((__m64 *)&spinor_out[1][1][0], xmm1);
  _mm_storeh_pi((__m64 *)&spinor_out[1][2][0], xmm2);

  /* Lower components - do projection */
  xmm0 = _mm_load_ps(&lower_sum[0][0][0]);
  xmm1 = _mm_load_ps(&lower_sum[1][0][0]);
  xmm2 = _mm_load_ps(&lower_sum[2][0][0]);
  

  xmm0 = _mm_sub_ps(xmm0, xmm3);
  xmm1 = _mm_sub_ps(xmm1, xmm4);
  xmm2 = _mm_sub_ps(xmm2, xmm5);
  
  /* Store */
  _mm_storel_pi((__m64 *)&spinor_out[2][0][0], xmm0);
  _mm_storel_pi((__m64 *)&spinor_out[2][1][0], xmm1);
  _mm_storel_pi((__m64 *)&spinor_out[2][2][0], xmm2);
  _mm_storeh_pi((__m64 *)&spinor_out[3][0][0], xmm0);
  _mm_storeh_pi((__m64 *)&spinor_out[3][1][0], xmm1);
  _mm_storeh_pi((__m64 *)&spinor_out[3][2][0], xmm2);

  }
  }
}
#endif

    




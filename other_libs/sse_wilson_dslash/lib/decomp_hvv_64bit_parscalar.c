#include <sse_config.h>

#include "decomp_hvv.h"
#include "sse_align.h"
#include "xmmintrin.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef  union { 
  unsigned int a[4];
  __m128d vector;
} SSEMask;
				
#ifdef SSEDSLASH_SLOPPY
#define SLOPPY_REGS				\
  __m128 xmm8 ALIGN;				\
  __m128 xmm9 ALIGN;				\
  __m128 xmm10 ALIGN

#define STORE(r)					\
  xmm8 = _mm_cvtpd_ps(xmm3);				\
  _mm_storel_pi((__m64*)(&dst[(r)][0][0]), xmm8);	\
  xmm9 = _mm_cvtpd_ps(xmm4);				\
  _mm_storel_pi((__m64*)(&dst[(r)][1][0]), xmm9);	\
  xmm10 = _mm_cvtpd_ps(xmm5);				\
  _mm_storel_pi((__m64*)(&dst[(r)][2][0]), xmm10)
#else
#define SLOPPY_REGS

#define STORE(r)				\
  _mm_store_pd(&dst[(r)][0][0], xmm3);		\
  _mm_store_pd(&dst[(r)][1][0], xmm4);		\
  _mm_store_pd(&dst[(r)][2][0], xmm5)
#endif

void decomp_hvv_gamma0_plus( spinor_array src, 
			     u_mat_array u,
			    halfspinor_array dst)
{
  __m128d xmm0 ALIGN;
  __m128d xmm1 ALIGN;
  __m128d xmm2 ALIGN;
  __m128d xmm3 ALIGN;
  __m128d xmm4 ALIGN;
  __m128d xmm5 ALIGN;
  __m128d xmm6 ALIGN;
  __m128d xmm7 ALIGN;

  SLOPPY_REGS;

  SSEMask sse_sgn = {{0x0, 0x80000000, 0x0,0x0 }};

  /* Projection: Munge components 0 & 3 */

  /* Load spinor component 0 */
  xmm0 = _mm_load_pd(&src[0][0][0]);
  xmm1 = _mm_load_pd(&src[0][1][0]);
  xmm2 = _mm_load_pd(&src[0][2][0]);

  xmm3 = _mm_load_pd(&src[3][0][0]);
  xmm4 = _mm_load_pd(&src[3][1][0]);
  xmm5 = _mm_load_pd(&src[3][2][0]);

  /* Shuffle the spinor components */
  xmm3 = _mm_shuffle_pd( xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd( xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd( xmm5, xmm5, 0x1);
  
  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_add_pd(xmm3, xmm0);
  xmm1 = _mm_add_pd(xmm4, xmm1);
  xmm2 = _mm_add_pd(xmm5, xmm2);

  /* Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0]);  /* _c11_re */
  xmm6 = _mm_load_sd(&u[0][1][0]);  /* _c21_re */
  xmm4 = _mm_load_sd(&u[1][0][0]);  /* _c12_re */
  xmm7 = _mm_load_sd(&u[1][2][0]);  /* _c32_re */
  xmm5 = _mm_load_sd(&u[2][0][0]);  /* _c13_re */

  xmm3 = _mm_unpacklo_pd(xmm3, xmm3); 
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm4 = _mm_unpacklo_pd(xmm4, xmm4);
  
  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd(xmm5, xmm5);
  
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);
  
  xmm6 = _mm_load_sd(&u[2][1][0]); /* _c23_re */
  xmm7 = _mm_load_sd(&u[0][2][0]); /* _c31_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0]); /* _c22_re */
  xmm7 = _mm_load_sd(&u[2][2][0]); /* _c33_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1]); /* _c11_im */
  xmm7 = _mm_load_sd(&u[1][1][1]); /* _c22_im */

  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[2][2][1] ); /* c33im */
  xmm7 = _mm_load_sd( &u[1][0][1] ); /* c12im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm2, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[0][1][1] ); /* c21im */
  xmm7 = _mm_load_sd( &u[2][0][1] ); /* c13im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd( &u[0][2][1] ); /* c31im */
  xmm6 = _mm_load_sd( &u[2][1][1] ); /* c23im */
  xmm7 = _mm_load_sd( &u[1][2][1] ); /* c32im */

  xmm0 = _mm_unpacklo_pd(xmm0, xmm0);  
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_mul_pd(xmm2, xmm0);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);

  xmm3 = _mm_add_pd(xmm0, xmm3);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  /* Store component 0  */
  STORE(0);


  /* Components 1 & 2 */
  xmm0 = _mm_load_pd(&src[1][0][0]);
  xmm1 = _mm_load_pd(&src[1][1][0]);
  xmm2 = _mm_load_pd(&src[1][2][0]);

  xmm3 = _mm_load_pd(&src[2][0][0]);
  xmm4 = _mm_load_pd(&src[2][1][0]);
  xmm5 = _mm_load_pd(&src[2][2][0]);

  xmm3 = _mm_shuffle_pd( xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd( xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd( xmm5, xmm5, 0x1);
  
  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);


  /* Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0]);  /* _c11_re */
  xmm6 = _mm_load_sd(&u[0][1][0]);  /* _c21_re */
  xmm4 = _mm_load_sd(&u[1][0][0]);  /* _c12_re */
  xmm7 = _mm_load_sd(&u[1][2][0]);  /* _c32_re */
  xmm5 = _mm_load_sd(&u[2][0][0]);  /* _c13_re */

  xmm3 = _mm_unpacklo_pd(xmm3, xmm3); 
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm4 = _mm_unpacklo_pd(xmm4, xmm4);
  
  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd(xmm5, xmm5);
  
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);
  
  xmm6 = _mm_load_sd(&u[2][1][0]); /* _c23_re */
  xmm7 = _mm_load_sd(&u[0][2][0]); /* _c31_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0]); /* _c22_re */
  xmm7 = _mm_load_sd(&u[2][2][0]); /* _c33_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1]); /* _c11_im */
  xmm7 = _mm_load_sd(&u[1][1][1]); /* _c22_im */

  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[2][2][1] ); /* c33im */
  xmm7 = _mm_load_sd( &u[1][0][1] ); /* c12im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm2, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[0][1][1] ); /* c21im */
  xmm7 = _mm_load_sd( &u[2][0][1] ); /* c13im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd( &u[0][2][1] ); /* c31im */
  xmm6 = _mm_load_sd( &u[2][1][1] ); /* c23im */
  xmm7 = _mm_load_sd( &u[1][2][1] ); /* c32im */

  xmm0 = _mm_unpacklo_pd(xmm0, xmm0);  
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_mul_pd(xmm2, xmm0);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);

  xmm3 = _mm_add_pd(xmm0, xmm3);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  STORE(1);

}

void decomp_hvv_gamma1_plus( spinor_array src, 
			     u_mat_array u,
			    halfspinor_array dst)
{
  __m128d xmm0 ALIGN;
  __m128d xmm1 ALIGN;
  __m128d xmm2 ALIGN;
  __m128d xmm3 ALIGN;
  __m128d xmm4 ALIGN;
  __m128d xmm5 ALIGN;
  __m128d xmm6 ALIGN;
  __m128d xmm7 ALIGN;
  
  SLOPPY_REGS;

  SSEMask sse_sgn = {{0x0, 0x80000000, 0x0,0x0 }};

  
  /* Projection: Munge components 0 & 3 */

  /* Load spinor component 0 */
  xmm0 = _mm_load_pd(&src[0][0][0]);
  xmm1 = _mm_load_pd(&src[0][1][0]);
  xmm2 = _mm_load_pd(&src[0][2][0]);

  xmm3 = _mm_load_pd(&src[3][0][0]);
  xmm4 = _mm_load_pd(&src[3][1][0]);
  xmm5 = _mm_load_pd(&src[3][2][0]);

  /* Shuffle the spinor components */

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);


  /* Multiply */
  /* Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0]);  /* _c11_re */
  xmm6 = _mm_load_sd(&u[0][1][0]);  /* _c21_re */
  xmm4 = _mm_load_sd(&u[1][0][0]);  /* _c12_re */
  xmm7 = _mm_load_sd(&u[1][2][0]);  /* _c32_re */
  xmm5 = _mm_load_sd(&u[2][0][0]);  /* _c13_re */

  xmm3 = _mm_unpacklo_pd(xmm3, xmm3); 
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm4 = _mm_unpacklo_pd(xmm4, xmm4);
  
  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd(xmm5, xmm5);
  
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);
  
  xmm6 = _mm_load_sd(&u[2][1][0]); /* _c23_re */
  xmm7 = _mm_load_sd(&u[0][2][0]); /* _c31_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0]); /* _c22_re */
  xmm7 = _mm_load_sd(&u[2][2][0]); /* _c33_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1]); /* _c11_im */
  xmm7 = _mm_load_sd(&u[1][1][1]); /* _c22_im */

  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[2][2][1] ); /* c33im */
  xmm7 = _mm_load_sd( &u[1][0][1] ); /* c12im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm2, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[0][1][1] ); /* c21im */
  xmm7 = _mm_load_sd( &u[2][0][1] ); /* c13im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd( &u[0][2][1] ); /* c31im */
  xmm6 = _mm_load_sd( &u[2][1][1] ); /* c23im */
  xmm7 = _mm_load_sd( &u[1][2][1] ); /* c32im */

  xmm0 = _mm_unpacklo_pd(xmm0, xmm0);  
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_mul_pd(xmm2, xmm0);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);

  xmm3 = _mm_add_pd(xmm0, xmm3);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  /* Store component 0  */
  STORE(0);

  /* Components 1 & 2 */
  xmm0 = _mm_load_pd(&src[1][0][0]);
  xmm1 = _mm_load_pd(&src[1][1][0]);
  xmm2 = _mm_load_pd(&src[1][2][0]);

  xmm3 = _mm_load_pd(&src[2][0][0]);
  xmm4 = _mm_load_pd(&src[2][1][0]);
  xmm5 = _mm_load_pd(&src[2][2][0]);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  /* Multiply */
  /* Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0]);  /* _c11_re */
  xmm6 = _mm_load_sd(&u[0][1][0]);  /* _c21_re */
  xmm4 = _mm_load_sd(&u[1][0][0]);  /* _c12_re */
  xmm7 = _mm_load_sd(&u[1][2][0]);  /* _c32_re */
  xmm5 = _mm_load_sd(&u[2][0][0]);  /* _c13_re */

  xmm3 = _mm_unpacklo_pd(xmm3, xmm3); 
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm4 = _mm_unpacklo_pd(xmm4, xmm4);
  
  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd(xmm5, xmm5);
  
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);
  
  xmm6 = _mm_load_sd(&u[2][1][0]); /* _c23_re */
  xmm7 = _mm_load_sd(&u[0][2][0]); /* _c31_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0]); /* _c22_re */
  xmm7 = _mm_load_sd(&u[2][2][0]); /* _c33_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1]); /* _c11_im */
  xmm7 = _mm_load_sd(&u[1][1][1]); /* _c22_im */

  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[2][2][1] ); /* c33im */
  xmm7 = _mm_load_sd( &u[1][0][1] ); /* c12im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm2, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[0][1][1] ); /* c21im */
  xmm7 = _mm_load_sd( &u[2][0][1] ); /* c13im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd( &u[0][2][1] ); /* c31im */
  xmm6 = _mm_load_sd( &u[2][1][1] ); /* c23im */
  xmm7 = _mm_load_sd( &u[1][2][1] ); /* c32im */

  xmm0 = _mm_unpacklo_pd(xmm0, xmm0);  
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_mul_pd(xmm2, xmm0);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);

  xmm3 = _mm_add_pd(xmm0, xmm3);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  /* Store component */

  STORE(1);

}

void decomp_hvv_gamma2_plus( spinor_array src, 
			     u_mat_array u,
			    halfspinor_array dst)
{
  __m128d xmm0 ALIGN;
  __m128d xmm1 ALIGN;
  __m128d xmm2 ALIGN;
  __m128d xmm3 ALIGN;
  __m128d xmm4 ALIGN;
  __m128d xmm5 ALIGN;
  __m128d xmm6 ALIGN;
  __m128d xmm7 ALIGN;

  SLOPPY_REGS;

  SSEMask sse_sgn = {{0x0, 0x80000000, 0x0 ,0x0 }};

  /* Projection: Munge components 0 & 2 */

  /* Load spinor component 0 */
  xmm0 = _mm_load_pd(&src[0][0][0]);
  xmm1 = _mm_load_pd(&src[0][1][0]);
  xmm2 = _mm_load_pd(&src[0][2][0]);

  xmm3 = _mm_load_pd(&src[2][0][0]);
  xmm4 = _mm_load_pd(&src[2][1][0]);
  xmm5 = _mm_load_pd(&src[2][2][0]);

  /* Shuffle the spinor components */
  xmm3 = _mm_shuffle_pd( xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd( xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd( xmm5, xmm5, 0x1);
  
  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  /* Multiply */
  /* Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0]);  /* _c11_re */
  xmm6 = _mm_load_sd(&u[0][1][0]);  /* _c21_re */
  xmm4 = _mm_load_sd(&u[1][0][0]);  /* _c12_re */
  xmm7 = _mm_load_sd(&u[1][2][0]);  /* _c32_re */
  xmm5 = _mm_load_sd(&u[2][0][0]);  /* _c13_re */

  xmm3 = _mm_unpacklo_pd(xmm3, xmm3); 
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm4 = _mm_unpacklo_pd(xmm4, xmm4);
  
  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd(xmm5, xmm5);
  
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);
  
  xmm6 = _mm_load_sd(&u[2][1][0]); /* _c23_re */
  xmm7 = _mm_load_sd(&u[0][2][0]); /* _c31_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0]); /* _c22_re */
  xmm7 = _mm_load_sd(&u[2][2][0]); /* _c33_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1]); /* _c11_im */
  xmm7 = _mm_load_sd(&u[1][1][1]); /* _c22_im */

  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[2][2][1] ); /* c33im */
  xmm7 = _mm_load_sd( &u[1][0][1] ); /* c12im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm2, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[0][1][1] ); /* c21im */
  xmm7 = _mm_load_sd( &u[2][0][1] ); /* c13im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd( &u[0][2][1] ); /* c31im */
  xmm6 = _mm_load_sd( &u[2][1][1] ); /* c23im */
  xmm7 = _mm_load_sd( &u[1][2][1] ); /* c32im */

  xmm0 = _mm_unpacklo_pd(xmm0, xmm0);  
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_mul_pd(xmm2, xmm0);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);

  xmm3 = _mm_add_pd(xmm0, xmm3);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  /* Store component 0  */
  STORE(0);

  /* Components 1 & 3 */
  xmm0 = _mm_load_pd(&src[1][0][0]);
  xmm1 = _mm_load_pd(&src[1][1][0]);
  xmm2 = _mm_load_pd(&src[1][2][0]);

  xmm3 = _mm_load_pd(&src[3][0][0]);
  xmm4 = _mm_load_pd(&src[3][1][0]);
  xmm5 = _mm_load_pd(&src[3][2][0]);

  xmm3 = _mm_shuffle_pd( xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd( xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd( xmm5, xmm5, 0x1);
  
  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  /* Multiply */
  /* Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0]);  /* _c11_re */
  xmm6 = _mm_load_sd(&u[0][1][0]);  /* _c21_re */
  xmm4 = _mm_load_sd(&u[1][0][0]);  /* _c12_re */
  xmm7 = _mm_load_sd(&u[1][2][0]);  /* _c32_re */
  xmm5 = _mm_load_sd(&u[2][0][0]);  /* _c13_re */

  xmm3 = _mm_unpacklo_pd(xmm3, xmm3); 
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm4 = _mm_unpacklo_pd(xmm4, xmm4);
  
  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd(xmm5, xmm5);
  
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);
  
  xmm6 = _mm_load_sd(&u[2][1][0]); /* _c23_re */
  xmm7 = _mm_load_sd(&u[0][2][0]); /* _c31_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0]); /* _c22_re */
  xmm7 = _mm_load_sd(&u[2][2][0]); /* _c33_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1]); /* _c11_im */
  xmm7 = _mm_load_sd(&u[1][1][1]); /* _c22_im */

  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[2][2][1] ); /* c33im */
  xmm7 = _mm_load_sd( &u[1][0][1] ); /* c12im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm2, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[0][1][1] ); /* c21im */
  xmm7 = _mm_load_sd( &u[2][0][1] ); /* c13im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd( &u[0][2][1] ); /* c31im */
  xmm6 = _mm_load_sd( &u[2][1][1] ); /* c23im */
  xmm7 = _mm_load_sd( &u[1][2][1] ); /* c32im */

  xmm0 = _mm_unpacklo_pd(xmm0, xmm0);  
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_mul_pd(xmm2, xmm0);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);

  xmm3 = _mm_add_pd(xmm0, xmm3);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  STORE(1);

}


void decomp_hvv_gamma3_plus( spinor_array src, 
			     u_mat_array u,
			    halfspinor_array dst)
{
  __m128d xmm0 ALIGN;
  __m128d xmm1 ALIGN;
  __m128d xmm2 ALIGN;
  __m128d xmm3 ALIGN;
  __m128d xmm4 ALIGN;
  __m128d xmm5 ALIGN;
  __m128d xmm6 ALIGN;
  __m128d xmm7 ALIGN;

  SLOPPY_REGS;

  SSEMask sse_sgn = {{0x0, 0x80000000, 0x0,0x0 }};

   /* Load spinor component 0 */
  xmm0 = _mm_load_pd(&src[0][0][0]);
  xmm1 = _mm_load_pd(&src[0][1][0]);
  xmm2 = _mm_load_pd(&src[0][2][0]);

  xmm3 = _mm_load_pd(&src[2][0][0]);
  xmm4 = _mm_load_pd(&src[2][1][0]);
  xmm5 = _mm_load_pd(&src[2][2][0]);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  /* Multiply */
  /* Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0]);  /* _c11_re */
  xmm6 = _mm_load_sd(&u[0][1][0]);  /* _c21_re */
  xmm4 = _mm_load_sd(&u[1][0][0]);  /* _c12_re */
  xmm7 = _mm_load_sd(&u[1][2][0]);  /* _c32_re */
  xmm5 = _mm_load_sd(&u[2][0][0]);  /* _c13_re */

  xmm3 = _mm_unpacklo_pd(xmm3, xmm3); 
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm4 = _mm_unpacklo_pd(xmm4, xmm4);
  
  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd(xmm5, xmm5);
  
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);
  
  xmm6 = _mm_load_sd(&u[2][1][0]); /* _c23_re */
  xmm7 = _mm_load_sd(&u[0][2][0]); /* _c31_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0]); /* _c22_re */
  xmm7 = _mm_load_sd(&u[2][2][0]); /* _c33_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1]); /* _c11_im */
  xmm7 = _mm_load_sd(&u[1][1][1]); /* _c22_im */

  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[2][2][1] ); /* c33im */
  xmm7 = _mm_load_sd( &u[1][0][1] ); /* c12im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm2, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[0][1][1] ); /* c21im */
  xmm7 = _mm_load_sd( &u[2][0][1] ); /* c13im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd( &u[0][2][1] ); /* c31im */
  xmm6 = _mm_load_sd( &u[2][1][1] ); /* c23im */
  xmm7 = _mm_load_sd( &u[1][2][1] ); /* c32im */

  xmm0 = _mm_unpacklo_pd(xmm0, xmm0);  
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_mul_pd(xmm2, xmm0);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);

  xmm3 = _mm_add_pd(xmm0, xmm3);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  /* Store component 0  */
  STORE(0);


  /* Load spinor component 0 */
  xmm0 = _mm_load_pd(&src[1][0][0]);
  xmm1 = _mm_load_pd(&src[1][1][0]);
  xmm2 = _mm_load_pd(&src[1][2][0]);

  xmm3 = _mm_load_pd(&src[3][0][0]);
  xmm4 = _mm_load_pd(&src[3][1][0]);
  xmm5 = _mm_load_pd(&src[3][2][0]);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);


  /* Multiply */
  /* Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0]);  /* _c11_re */
  xmm6 = _mm_load_sd(&u[0][1][0]);  /* _c21_re */
  xmm4 = _mm_load_sd(&u[1][0][0]);  /* _c12_re */
  xmm7 = _mm_load_sd(&u[1][2][0]);  /* _c32_re */
  xmm5 = _mm_load_sd(&u[2][0][0]);  /* _c13_re */

  xmm3 = _mm_unpacklo_pd(xmm3, xmm3); 
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm4 = _mm_unpacklo_pd(xmm4, xmm4);
  
  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd(xmm5, xmm5);
  
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);
  
  xmm6 = _mm_load_sd(&u[2][1][0]); /* _c23_re */
  xmm7 = _mm_load_sd(&u[0][2][0]); /* _c31_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0]); /* _c22_re */
  xmm7 = _mm_load_sd(&u[2][2][0]); /* _c33_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1]); /* _c11_im */
  xmm7 = _mm_load_sd(&u[1][1][1]); /* _c22_im */

  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[2][2][1] ); /* c33im */
  xmm7 = _mm_load_sd( &u[1][0][1] ); /* c12im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm2, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[0][1][1] ); /* c21im */
  xmm7 = _mm_load_sd( &u[2][0][1] ); /* c13im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd( &u[0][2][1] ); /* c31im */
  xmm6 = _mm_load_sd( &u[2][1][1] ); /* c23im */
  xmm7 = _mm_load_sd( &u[1][2][1] ); /* c32im */

  xmm0 = _mm_unpacklo_pd(xmm0, xmm0);  
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_mul_pd(xmm2, xmm0);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);

  xmm3 = _mm_add_pd(xmm0, xmm3);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);
  
  /* Store component */
  STORE(1);
}


void decomp_hvv_gamma0_minus( spinor_array src, 
			     u_mat_array u,
			    halfspinor_array dst)
{
  __m128d xmm0 ALIGN;
  __m128d xmm1 ALIGN;
  __m128d xmm2 ALIGN;
  __m128d xmm3 ALIGN;
  __m128d xmm4 ALIGN;
  __m128d xmm5 ALIGN;
  __m128d xmm6 ALIGN;
  __m128d xmm7 ALIGN;
 
  SLOPPY_REGS;


  SSEMask sse_sgn = {{0x0, 0x80000000, 0x0,0x0 }};

  /* Projection: Munge components 0 & 3 */

  /* Load spinor component 0 */
  xmm0 = _mm_load_pd(&src[0][0][0]);
  xmm1 = _mm_load_pd(&src[0][1][0]);
  xmm2 = _mm_load_pd(&src[0][2][0]);

  xmm3 = _mm_load_pd(&src[3][0][0]);
  xmm4 = _mm_load_pd(&src[3][1][0]);
  xmm5 = _mm_load_pd(&src[3][2][0]);

  /* Shuffle the spinor components */
  xmm3 = _mm_shuffle_pd( xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd( xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd( xmm5, xmm5, 0x1);
  
  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  /* Multiply */
  /* Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0]);  /* _c11_re */
  xmm6 = _mm_load_sd(&u[0][1][0]);  /* _c21_re */
  xmm4 = _mm_load_sd(&u[1][0][0]);  /* _c12_re */
  xmm7 = _mm_load_sd(&u[1][2][0]);  /* _c32_re */
  xmm5 = _mm_load_sd(&u[2][0][0]);  /* _c13_re */

  xmm3 = _mm_unpacklo_pd(xmm3, xmm3); 
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm4 = _mm_unpacklo_pd(xmm4, xmm4);
  
  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd(xmm5, xmm5);
  
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);
  
  xmm6 = _mm_load_sd(&u[2][1][0]); /* _c23_re */
  xmm7 = _mm_load_sd(&u[0][2][0]); /* _c31_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0]); /* _c22_re */
  xmm7 = _mm_load_sd(&u[2][2][0]); /* _c33_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1]); /* _c11_im */
  xmm7 = _mm_load_sd(&u[1][1][1]); /* _c22_im */

  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[2][2][1] ); /* c33im */
  xmm7 = _mm_load_sd( &u[1][0][1] ); /* c12im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm2, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[0][1][1] ); /* c21im */
  xmm7 = _mm_load_sd( &u[2][0][1] ); /* c13im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd( &u[0][2][1] ); /* c31im */
  xmm6 = _mm_load_sd( &u[2][1][1] ); /* c23im */
  xmm7 = _mm_load_sd( &u[1][2][1] ); /* c32im */

  xmm0 = _mm_unpacklo_pd(xmm0, xmm0);  
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_mul_pd(xmm2, xmm0);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);

  xmm3 = _mm_add_pd(xmm0, xmm3);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  STORE(0);

  /* Components 1 & 2 */
  xmm0 = _mm_load_pd(&src[1][0][0]);
  xmm1 = _mm_load_pd(&src[1][1][0]);
  xmm2 = _mm_load_pd(&src[1][2][0]);

  xmm3 = _mm_load_pd(&src[2][0][0]);
  xmm4 = _mm_load_pd(&src[2][1][0]);
  xmm5 = _mm_load_pd(&src[2][2][0]);

  xmm3 = _mm_shuffle_pd( xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd( xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd( xmm5, xmm5, 0x1);
  
  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  /* Multiply */
  /* Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0]);  /* _c11_re */
  xmm6 = _mm_load_sd(&u[0][1][0]);  /* _c21_re */
  xmm4 = _mm_load_sd(&u[1][0][0]);  /* _c12_re */
  xmm7 = _mm_load_sd(&u[1][2][0]);  /* _c32_re */
  xmm5 = _mm_load_sd(&u[2][0][0]);  /* _c13_re */

  xmm3 = _mm_unpacklo_pd(xmm3, xmm3); 
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm4 = _mm_unpacklo_pd(xmm4, xmm4);
  
  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd(xmm5, xmm5);
  
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);
  
  xmm6 = _mm_load_sd(&u[2][1][0]); /* _c23_re */
  xmm7 = _mm_load_sd(&u[0][2][0]); /* _c31_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0]); /* _c22_re */
  xmm7 = _mm_load_sd(&u[2][2][0]); /* _c33_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1]); /* _c11_im */
  xmm7 = _mm_load_sd(&u[1][1][1]); /* _c22_im */

  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[2][2][1] ); /* c33im */
  xmm7 = _mm_load_sd( &u[1][0][1] ); /* c12im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm2, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[0][1][1] ); /* c21im */
  xmm7 = _mm_load_sd( &u[2][0][1] ); /* c13im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd( &u[0][2][1] ); /* c31im */
  xmm6 = _mm_load_sd( &u[2][1][1] ); /* c23im */
  xmm7 = _mm_load_sd( &u[1][2][1] ); /* c32im */

  xmm0 = _mm_unpacklo_pd(xmm0, xmm0);  
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_mul_pd(xmm2, xmm0);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);

  xmm3 = _mm_add_pd(xmm0, xmm3);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  /* Store component 1 */
  STORE(1);

}

void decomp_hvv_gamma1_minus( spinor_array src, 
			     u_mat_array u,
			    halfspinor_array dst)
{
  __m128d xmm0 ALIGN;
  __m128d xmm1 ALIGN;
  __m128d xmm2 ALIGN;
  __m128d xmm3 ALIGN;
  __m128d xmm4 ALIGN;
  __m128d xmm5 ALIGN;
  __m128d xmm6 ALIGN;
  __m128d xmm7 ALIGN;
  
  SLOPPY_REGS;


  SSEMask sse_sgn = {{0x0, 0x80000000, 0x0,0x0 }};
  
  /* Projection: Munge components 0 & 3 */

  /* Load spinor component 0 */
  xmm0 = _mm_load_pd(&src[0][0][0]);
  xmm1 = _mm_load_pd(&src[0][1][0]);
  xmm2 = _mm_load_pd(&src[0][2][0]);

  xmm3 = _mm_load_pd(&src[3][0][0]);
  xmm4 = _mm_load_pd(&src[3][1][0]);
  xmm5 = _mm_load_pd(&src[3][2][0]);

  /* Shuffle the spinor components */

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  /* Multiply */
  /* Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0]);  /* _c11_re */
  xmm6 = _mm_load_sd(&u[0][1][0]);  /* _c21_re */
  xmm4 = _mm_load_sd(&u[1][0][0]);  /* _c12_re */
  xmm7 = _mm_load_sd(&u[1][2][0]);  /* _c32_re */
  xmm5 = _mm_load_sd(&u[2][0][0]);  /* _c13_re */

  xmm3 = _mm_unpacklo_pd(xmm3, xmm3); 
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm4 = _mm_unpacklo_pd(xmm4, xmm4);
  
  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd(xmm5, xmm5);
  
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);
  
  xmm6 = _mm_load_sd(&u[2][1][0]); /* _c23_re */
  xmm7 = _mm_load_sd(&u[0][2][0]); /* _c31_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0]); /* _c22_re */
  xmm7 = _mm_load_sd(&u[2][2][0]); /* _c33_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1]); /* _c11_im */
  xmm7 = _mm_load_sd(&u[1][1][1]); /* _c22_im */

  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[2][2][1] ); /* c33im */
  xmm7 = _mm_load_sd( &u[1][0][1] ); /* c12im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm2, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[0][1][1] ); /* c21im */
  xmm7 = _mm_load_sd( &u[2][0][1] ); /* c13im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd( &u[0][2][1] ); /* c31im */
  xmm6 = _mm_load_sd( &u[2][1][1] ); /* c23im */
  xmm7 = _mm_load_sd( &u[1][2][1] ); /* c32im */

  xmm0 = _mm_unpacklo_pd(xmm0, xmm0);  
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_mul_pd(xmm2, xmm0);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);

  xmm3 = _mm_add_pd(xmm0, xmm3);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);


  /* Store component 0  */
  STORE(0);

  /* Components 1 & 2 */
  xmm0 = _mm_load_pd(&src[1][0][0]);
  xmm1 = _mm_load_pd(&src[1][1][0]);
  xmm2 = _mm_load_pd(&src[1][2][0]);

  xmm3 = _mm_load_pd(&src[2][0][0]);
  xmm4 = _mm_load_pd(&src[2][1][0]);
  xmm5 = _mm_load_pd(&src[2][2][0]);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  /* Multiply */
  /* Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0]);  /* _c11_re */
  xmm6 = _mm_load_sd(&u[0][1][0]);  /* _c21_re */
  xmm4 = _mm_load_sd(&u[1][0][0]);  /* _c12_re */
  xmm7 = _mm_load_sd(&u[1][2][0]);  /* _c32_re */
  xmm5 = _mm_load_sd(&u[2][0][0]);  /* _c13_re */

  xmm3 = _mm_unpacklo_pd(xmm3, xmm3); 
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm4 = _mm_unpacklo_pd(xmm4, xmm4);
  
  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd(xmm5, xmm5);
  
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);
  
  xmm6 = _mm_load_sd(&u[2][1][0]); /* _c23_re */
  xmm7 = _mm_load_sd(&u[0][2][0]); /* _c31_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0]); /* _c22_re */
  xmm7 = _mm_load_sd(&u[2][2][0]); /* _c33_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1]); /* _c11_im */
  xmm7 = _mm_load_sd(&u[1][1][1]); /* _c22_im */

  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[2][2][1] ); /* c33im */
  xmm7 = _mm_load_sd( &u[1][0][1] ); /* c12im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm2, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[0][1][1] ); /* c21im */
  xmm7 = _mm_load_sd( &u[2][0][1] ); /* c13im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd( &u[0][2][1] ); /* c31im */
  xmm6 = _mm_load_sd( &u[2][1][1] ); /* c23im */
  xmm7 = _mm_load_sd( &u[1][2][1] ); /* c32im */

  xmm0 = _mm_unpacklo_pd(xmm0, xmm0);  
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_mul_pd(xmm2, xmm0);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);

  xmm3 = _mm_add_pd(xmm0, xmm3);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  /* Store component 1 */
  STORE(1);

}

void decomp_hvv_gamma2_minus( spinor_array src, 
			     u_mat_array u,
			    halfspinor_array dst)
{
  __m128d xmm0 ALIGN;
  __m128d xmm1 ALIGN;
  __m128d xmm2 ALIGN;

  __m128d xmm3 ALIGN;
  __m128d xmm4 ALIGN;
  __m128d xmm5 ALIGN;
  __m128d xmm6 ALIGN;
  __m128d xmm7 ALIGN;
 
  SLOPPY_REGS;

  SSEMask sse_sgn = {{0x0, 0x80000000, 0x0, 0x0 }};

  /* Projection: Munge components 0 & 2 */

  /* Load spinor component 0 */
  xmm0 = _mm_load_pd(&src[0][0][0]);
  xmm1 = _mm_load_pd(&src[0][1][0]);
  xmm2 = _mm_load_pd(&src[0][2][0]);

  xmm3 = _mm_load_pd(&src[2][0][0]);
  xmm4 = _mm_load_pd(&src[2][1][0]);
  xmm5 = _mm_load_pd(&src[2][2][0]);

  /* Shuffle the spinor components */
  xmm3 = _mm_shuffle_pd( xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd( xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd( xmm5, xmm5, 0x1);
  
  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  /* Multiply */
  /* Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0]);  /* _c11_re */
  xmm6 = _mm_load_sd(&u[0][1][0]);  /* _c21_re */
  xmm4 = _mm_load_sd(&u[1][0][0]);  /* _c12_re */
  xmm7 = _mm_load_sd(&u[1][2][0]);  /* _c32_re */
  xmm5 = _mm_load_sd(&u[2][0][0]);  /* _c13_re */

  xmm3 = _mm_unpacklo_pd(xmm3, xmm3); 
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm4 = _mm_unpacklo_pd(xmm4, xmm4);
  
  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd(xmm5, xmm5);
  
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);
  
  xmm6 = _mm_load_sd(&u[2][1][0]); /* _c23_re */
  xmm7 = _mm_load_sd(&u[0][2][0]); /* _c31_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0]); /* _c22_re */
  xmm7 = _mm_load_sd(&u[2][2][0]); /* _c33_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1]); /* _c11_im */
  xmm7 = _mm_load_sd(&u[1][1][1]); /* _c22_im */

  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[2][2][1] ); /* c33im */
  xmm7 = _mm_load_sd( &u[1][0][1] ); /* c12im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm2, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[0][1][1] ); /* c21im */
  xmm7 = _mm_load_sd( &u[2][0][1] ); /* c13im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd( &u[0][2][1] ); /* c31im */
  xmm6 = _mm_load_sd( &u[2][1][1] ); /* c23im */
  xmm7 = _mm_load_sd( &u[1][2][1] ); /* c32im */

  xmm0 = _mm_unpacklo_pd(xmm0, xmm0);  
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_mul_pd(xmm2, xmm0);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);

  xmm3 = _mm_add_pd(xmm0, xmm3);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  /* Store component 0  */
  STORE(0);

  /* Components 1 & 3 */
  xmm0 = _mm_load_pd(&src[1][0][0]);
  xmm1 = _mm_load_pd(&src[1][1][0]);
  xmm2 = _mm_load_pd(&src[1][2][0]);

  xmm3 = _mm_load_pd(&src[3][0][0]);
  xmm4 = _mm_load_pd(&src[3][1][0]);
  xmm5 = _mm_load_pd(&src[3][2][0]);

  xmm3 = _mm_shuffle_pd( xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd( xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd( xmm5, xmm5, 0x1);
  
  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  /* Multiply */
  /* Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0]);  /* _c11_re */
  xmm6 = _mm_load_sd(&u[0][1][0]);  /* _c21_re */
  xmm4 = _mm_load_sd(&u[1][0][0]);  /* _c12_re */
  xmm7 = _mm_load_sd(&u[1][2][0]);  /* _c32_re */
  xmm5 = _mm_load_sd(&u[2][0][0]);  /* _c13_re */

  xmm3 = _mm_unpacklo_pd(xmm3, xmm3); 
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm4 = _mm_unpacklo_pd(xmm4, xmm4);
  
  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd(xmm5, xmm5);
  
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);
  
  xmm6 = _mm_load_sd(&u[2][1][0]); /* _c23_re */
  xmm7 = _mm_load_sd(&u[0][2][0]); /* _c31_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0]); /* _c22_re */
  xmm7 = _mm_load_sd(&u[2][2][0]); /* _c33_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1]); /* _c11_im */
  xmm7 = _mm_load_sd(&u[1][1][1]); /* _c22_im */

  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[2][2][1] ); /* c33im */
  xmm7 = _mm_load_sd( &u[1][0][1] ); /* c12im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm2, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[0][1][1] ); /* c21im */
  xmm7 = _mm_load_sd( &u[2][0][1] ); /* c13im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd( &u[0][2][1] ); /* c31im */
  xmm6 = _mm_load_sd( &u[2][1][1] ); /* c23im */
  xmm7 = _mm_load_sd( &u[1][2][1] ); /* c32im */

  xmm0 = _mm_unpacklo_pd(xmm0, xmm0);  
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_mul_pd(xmm2, xmm0);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);

  xmm3 = _mm_add_pd(xmm0, xmm3);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  /* Store component */
  STORE(1);

}


void decomp_hvv_gamma3_minus( spinor_array src, 
			     u_mat_array u,
			    halfspinor_array dst)
{
  __m128d xmm0 ALIGN;
  __m128d xmm1 ALIGN;
  __m128d xmm2 ALIGN;
  __m128d xmm3 ALIGN;
  __m128d xmm4 ALIGN;
  __m128d xmm5 ALIGN;
  __m128d xmm6 ALIGN;
  __m128d xmm7 ALIGN;
  
  SLOPPY_REGS;

  SSEMask sse_sgn = {{0x0, 0x80000000, 0x0,0x0 }};
  
   /* Load spinor component 0 */
  xmm0 = _mm_load_pd(&src[0][0][0]);
  xmm1 = _mm_load_pd(&src[0][1][0]);
  xmm2 = _mm_load_pd(&src[0][2][0]);

  xmm3 = _mm_load_pd(&src[2][0][0]);
  xmm4 = _mm_load_pd(&src[2][1][0]);
  xmm5 = _mm_load_pd(&src[2][2][0]);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  /* Multiply */
  /* Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0]);  /* _c11_re */
  xmm6 = _mm_load_sd(&u[0][1][0]);  /* _c21_re */
  xmm4 = _mm_load_sd(&u[1][0][0]);  /* _c12_re */
  xmm7 = _mm_load_sd(&u[1][2][0]);  /* _c32_re */
  xmm5 = _mm_load_sd(&u[2][0][0]);  /* _c13_re */

  xmm3 = _mm_unpacklo_pd(xmm3, xmm3); 
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm4 = _mm_unpacklo_pd(xmm4, xmm4);
  
  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd(xmm5, xmm5);
  
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);
  
  xmm6 = _mm_load_sd(&u[2][1][0]); /* _c23_re */
  xmm7 = _mm_load_sd(&u[0][2][0]); /* _c31_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0]); /* _c22_re */
  xmm7 = _mm_load_sd(&u[2][2][0]); /* _c33_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1]); /* _c11_im */
  xmm7 = _mm_load_sd(&u[1][1][1]); /* _c22_im */

  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[2][2][1] ); /* c33im */
  xmm7 = _mm_load_sd( &u[1][0][1] ); /* c12im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm2, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[0][1][1] ); /* c21im */
  xmm7 = _mm_load_sd( &u[2][0][1] ); /* c13im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd( &u[0][2][1] ); /* c31im */
  xmm6 = _mm_load_sd( &u[2][1][1] ); /* c23im */
  xmm7 = _mm_load_sd( &u[1][2][1] ); /* c32im */

  xmm0 = _mm_unpacklo_pd(xmm0, xmm0);  
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_mul_pd(xmm2, xmm0);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);

  xmm3 = _mm_add_pd(xmm0, xmm3);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  
  /* Store component 0 */
  STORE(0);

  /* Load spinor component 0 */
  xmm0 = _mm_load_pd(&src[1][0][0]);
  xmm1 = _mm_load_pd(&src[1][1][0]);
  xmm2 = _mm_load_pd(&src[1][2][0]);

  xmm3 = _mm_load_pd(&src[3][0][0]);
  xmm4 = _mm_load_pd(&src[3][1][0]);
  xmm5 = _mm_load_pd(&src[3][2][0]);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);


  /* Multiply */
  /* Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0]);  /* _c11_re */
  xmm6 = _mm_load_sd(&u[0][1][0]);  /* _c21_re */
  xmm4 = _mm_load_sd(&u[1][0][0]);  /* _c12_re */
  xmm7 = _mm_load_sd(&u[1][2][0]);  /* _c32_re */
  xmm5 = _mm_load_sd(&u[2][0][0]);  /* _c13_re */

  xmm3 = _mm_unpacklo_pd(xmm3, xmm3); 
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm4 = _mm_unpacklo_pd(xmm4, xmm4);
  
  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd(xmm5, xmm5);
  
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);
  
  xmm6 = _mm_load_sd(&u[2][1][0]); /* _c23_re */
  xmm7 = _mm_load_sd(&u[0][2][0]); /* _c31_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0]); /* _c22_re */
  xmm7 = _mm_load_sd(&u[2][2][0]); /* _c33_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1]); /* _c11_im */
  xmm7 = _mm_load_sd(&u[1][1][1]); /* _c22_im */

  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[2][2][1] ); /* c33im */
  xmm7 = _mm_load_sd( &u[1][0][1] ); /* c12im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm2, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[0][1][1] ); /* c21im */
  xmm7 = _mm_load_sd( &u[2][0][1] ); /* c13im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd( &u[0][2][1] ); /* c31im */
  xmm6 = _mm_load_sd( &u[2][1][1] ); /* c23im */
  xmm7 = _mm_load_sd( &u[1][2][1] ); /* c32im */

  xmm0 = _mm_unpacklo_pd(xmm0, xmm0);  
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_mul_pd(xmm2, xmm0);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);

  xmm3 = _mm_add_pd(xmm0, xmm3);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);


  /* Store component 1 */
  STORE(1);

}

#undef SLOPPY_REGS
#undef STORE

#ifdef __cplusplus
};
#endif

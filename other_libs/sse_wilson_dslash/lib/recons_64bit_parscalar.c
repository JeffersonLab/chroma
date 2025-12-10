#include "recons.h"
#include "xmmintrin.h"
#include "sse_align.h"
#include "types64.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef  union { 
    unsigned int a[4];
    __m128d vector;
  } SSEMask;


#ifdef SSEDSLASH_SLOPPY
#define SLOPPY_REGS   \
  __m128 xmm6 ALIGN; \
  __m128 xmm7 ALIGN; \
  __m128 xmm8 ALIGN

#else
#define SLOPPY_REGS 
#endif

#ifdef SSEDSLASH_SLOPPY
#define LOAD( var, row )					\
  xmm6 = _mm_loadl_pi(xmm6,(__m64*)(&(var)[(row)][0][0]) );	\
  xmm3 = _mm_cvtps_pd(xmm6);					\
  								\
  xmm7 = _mm_loadl_pi(xmm7,(__m64*)(&(var)[(row)][1][0]) );	\
  xmm4 = _mm_cvtps_pd(xmm7);				\
  							\
  xmm8 = _mm_loadl_pi(xmm8,(__m64*)(&(var)[(row)][2][0]) );	\
  xmm5 = _mm_cvtps_pd(xmm8)
#else
#define LOAD( var, row ) \
  xmm3 = _mm_load_pd( &(var)[(row)][0][0] );	\
  xmm4 = _mm_load_pd( &(var)[(row)][1][0] );	\
  xmm5 = _mm_load_pd( &(var)[(row)][2][0] )
#endif

void recons_4dir_plus( halfspinor_array hs0,
		       halfspinor_array hs1,
		       halfspinor_array hs2,
		       halfspinor_array hs3,
		      spinor_array spinor)
{
  __m128d xmm0 ALIGN;
  __m128d xmm1 ALIGN;  
  __m128d xmm2 ALIGN;
  __m128d xmm3 ALIGN;
  __m128d xmm4 ALIGN;
  __m128d xmm5 ALIGN;

  SSEMask sse_sgn ALIGN = {{0x0, 0x80000000, 0x0,0x0 }};

  SLOPPY_REGS ;

  /* Component 0 */
  LOAD(hs0, 0);

  xmm0 = _mm_load_pd( &spinor[0][0][0] );
  xmm1 = _mm_load_pd( &spinor[0][1][0] );
  xmm2 = _mm_load_pd( &spinor[0][2][0] );

  
  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  LOAD(hs1, 0);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  LOAD(hs2, 0);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  LOAD(hs3, 0);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  _mm_store_pd(&spinor[0][0][0], xmm0);
  _mm_store_pd(&spinor[0][1][0], xmm1);
  _mm_store_pd(&spinor[0][2][0], xmm2);

  /* Component 1 */
  LOAD(hs0,1);

  xmm0 = _mm_load_pd( &spinor[1][0][0] );
  xmm1 = _mm_load_pd( &spinor[1][1][0] );
  xmm2 = _mm_load_pd( &spinor[1][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  LOAD(hs1,1);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  LOAD(hs2, 1);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  LOAD(hs3, 1);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  _mm_store_pd(&spinor[1][0][0], xmm0);
  _mm_store_pd(&spinor[1][1][0], xmm1);
  _mm_store_pd(&spinor[1][2][0], xmm2);

  /* Component 2 */

  LOAD(hs0, 1);

  xmm0 = _mm_load_pd( &spinor[2][0][0] );
  xmm1 = _mm_load_pd( &spinor[2][1][0] );
  xmm2 = _mm_load_pd( &spinor[2][2][0] );


  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1); 

  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  LOAD(hs1, 1);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  LOAD(hs2, 0);

  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1);

  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  LOAD(hs3, 0);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);


  _mm_store_pd(&spinor[2][0][0], xmm0);
  _mm_store_pd(&spinor[2][1][0], xmm1);
  _mm_store_pd(&spinor[2][2][0], xmm2);


  LOAD(hs0, 0);

  xmm0 = _mm_load_pd( &spinor[3][0][0] );
  xmm1 = _mm_load_pd( &spinor[3][1][0] );
  xmm2 = _mm_load_pd( &spinor[3][2][0] );

  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1);


  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);


  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  LOAD(hs1, 0);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  LOAD(hs2, 1);

  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1);

  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_add_pd(xmm3, xmm0);
  xmm1 = _mm_add_pd(xmm4, xmm1);
  xmm2 = _mm_add_pd(xmm5, xmm2);

  LOAD(hs3, 1);

  xmm0 = _mm_add_pd(xmm3, xmm0);
  xmm1 = _mm_add_pd(xmm4, xmm1);
  xmm2 = _mm_add_pd(xmm5, xmm2);

  _mm_store_pd(&spinor[3][0][0], xmm0);
  _mm_store_pd(&spinor[3][1][0], xmm1);
  _mm_store_pd(&spinor[3][2][0], xmm2);
 


}

void recons_3dir_plus( halfspinor_array hs0,
		       halfspinor_array hs1,
		       halfspinor_array hs2,
		       spinor_array spinor)
{
  __m128d xmm0 ALIGN;
  __m128d xmm1 ALIGN;  
  __m128d xmm2 ALIGN;
  __m128d xmm3 ALIGN;
  __m128d xmm4 ALIGN;
  __m128d xmm5 ALIGN;

  SSEMask sse_sgn ALIGN = {{0x0, 0x80000000, 0x0,0x0 }};
  SLOPPY_REGS ; 

  /* Component 0 */
  LOAD( hs0, 0 );

  xmm0 = _mm_load_pd( &spinor[0][0][0] );
  xmm1 = _mm_load_pd( &spinor[0][1][0] );
  xmm2 = _mm_load_pd( &spinor[0][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  LOAD(hs1, 0);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  LOAD(hs2, 0);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  _mm_store_pd(&spinor[0][0][0], xmm0);
  _mm_store_pd(&spinor[0][1][0], xmm1);
  _mm_store_pd(&spinor[0][2][0], xmm2);

  /* Component 1 */
  LOAD( hs0, 1 );

  xmm0 = _mm_load_pd( &spinor[1][0][0] );
  xmm1 = _mm_load_pd( &spinor[1][1][0] );
  xmm2 = _mm_load_pd( &spinor[1][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  LOAD( hs1, 1 );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  LOAD( hs2, 1 );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);


  _mm_store_pd(&spinor[1][0][0], xmm0);
  _mm_store_pd(&spinor[1][1][0], xmm1);
  _mm_store_pd(&spinor[1][2][0], xmm2);

  /* Component 2 */

  LOAD( hs0, 1 );

  xmm0 = _mm_load_pd( &spinor[2][0][0] );
  xmm1 = _mm_load_pd( &spinor[2][1][0] );
  xmm2 = _mm_load_pd( &spinor[2][2][0] );


  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1); 

  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  LOAD( hs1, 1);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);


  LOAD( hs2, 0 );

  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1);

  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  _mm_store_pd(&spinor[2][0][0], xmm0);
  _mm_store_pd(&spinor[2][1][0], xmm1);
  _mm_store_pd(&spinor[2][2][0], xmm2);

  LOAD(hs0, 0);

  xmm0 = _mm_load_pd( &spinor[3][0][0] );
  xmm1 = _mm_load_pd( &spinor[3][1][0] );
  xmm2 = _mm_load_pd( &spinor[3][2][0] );

  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1);


  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);


  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  LOAD(hs1, 0);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  LOAD(hs2, 1);

  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1);

  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_add_pd(xmm3, xmm0);
  xmm1 = _mm_add_pd(xmm4, xmm1);
  xmm2 = _mm_add_pd(xmm5, xmm2);


  _mm_store_pd(&spinor[3][0][0], xmm0);
  _mm_store_pd(&spinor[3][1][0], xmm1);
  _mm_store_pd(&spinor[3][2][0], xmm2);

}

void recons_4dir_minus( halfspinor_array hs0,
		        halfspinor_array hs1,
		        halfspinor_array hs2,
		        halfspinor_array hs3,
		       spinor_array spinor)
{
  __m128d xmm0 ALIGN;
  __m128d xmm1 ALIGN;  
  __m128d xmm2 ALIGN;
  __m128d xmm3 ALIGN;
  __m128d xmm4 ALIGN;
  __m128d xmm5 ALIGN;


  SSEMask sse_sgn ALIGN = {{0x0, 0x80000000, 0x0,0x0 }};

  SLOPPY_REGS; 
  /* Component 0 */
  LOAD( hs0, 0 );

  xmm0 = _mm_load_pd( &spinor[0][0][0] );
  xmm1 = _mm_load_pd( &spinor[0][1][0] );
  xmm2 = _mm_load_pd( &spinor[0][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  LOAD( hs1, 0 );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  LOAD( hs2, 0 );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  LOAD( hs3, 0 );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  _mm_store_pd(&spinor[0][0][0], xmm0);
  _mm_store_pd(&spinor[0][1][0], xmm1);
  _mm_store_pd(&spinor[0][2][0], xmm2);

  /* Component 1 */

  LOAD( hs0, 1 );

  xmm0 = _mm_load_pd( &spinor[1][0][0] );
  xmm1 = _mm_load_pd( &spinor[1][1][0] );
  xmm2 = _mm_load_pd( &spinor[1][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  LOAD( hs1, 1 );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  LOAD( hs2, 1 );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  LOAD( hs3, 1 );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  _mm_store_pd(&spinor[1][0][0], xmm0);
  _mm_store_pd(&spinor[1][1][0], xmm1);
  _mm_store_pd(&spinor[1][2][0], xmm2);

  /* Component 2 */

  LOAD( hs0, 1 );

  xmm0 = _mm_load_pd( &spinor[2][0][0] );
  xmm1 = _mm_load_pd( &spinor[2][1][0] );
  xmm2 = _mm_load_pd( &spinor[2][2][0] );


  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1); 

  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  LOAD( hs1, 1 );

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  LOAD( hs2, 0 );

  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1);

  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  LOAD( hs3, 0);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);


  _mm_store_pd(&spinor[2][0][0], xmm0);
  _mm_store_pd(&spinor[2][1][0], xmm1);
  _mm_store_pd(&spinor[2][2][0], xmm2);

  /* Component 3 */
  LOAD( hs0, 0 );

  xmm0 = _mm_load_pd( &spinor[3][0][0] );
  xmm1 = _mm_load_pd( &spinor[3][1][0] );
  xmm2 = _mm_load_pd( &spinor[3][2][0] );

  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1);


  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);


  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  LOAD( hs1, 0 );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  LOAD( hs2, 1 );

  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1);

  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  LOAD( hs3, 1 );

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  _mm_store_pd(&spinor[3][0][0], xmm0);
  _mm_store_pd(&spinor[3][1][0], xmm1);
  _mm_store_pd(&spinor[3][2][0], xmm2);
 


}

void recons_3dir_minus( halfspinor_array hs0,
		        halfspinor_array hs1,
		        halfspinor_array hs2,
			spinor_array spinor)
{
  __m128d xmm0 ALIGN;
  __m128d xmm1 ALIGN;  
  __m128d xmm2 ALIGN;
  __m128d xmm3 ALIGN;
  __m128d xmm4 ALIGN;
  __m128d xmm5 ALIGN;


  SSEMask sse_sgn ALIGN = {{0x0, 0x80000000, 0x0,0x0 }};

  SLOPPY_REGS;

  /* Component 0 */
  LOAD(hs0, 0);

  xmm0 = _mm_load_pd( &spinor[0][0][0] );
  xmm1 = _mm_load_pd( &spinor[0][1][0] );
  xmm2 = _mm_load_pd( &spinor[0][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  LOAD(hs1, 0);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  LOAD(hs2, 0);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);


  _mm_store_pd(&spinor[0][0][0], xmm0);
  _mm_store_pd(&spinor[0][1][0], xmm1);
  _mm_store_pd(&spinor[0][2][0], xmm2);

  /* Component 1 */
  LOAD(hs0, 1);

  xmm0 = _mm_load_pd( &spinor[1][0][0] );
  xmm1 = _mm_load_pd( &spinor[1][1][0] );
  xmm2 = _mm_load_pd( &spinor[1][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  LOAD(hs1, 1);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  LOAD(hs2, 1);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);


  _mm_store_pd(&spinor[1][0][0], xmm0);
  _mm_store_pd(&spinor[1][1][0], xmm1);
  _mm_store_pd(&spinor[1][2][0], xmm2);

  /* Component 2 */

  LOAD(hs0, 1);

  xmm0 = _mm_load_pd( &spinor[2][0][0] );
  xmm1 = _mm_load_pd( &spinor[2][1][0] );
  xmm2 = _mm_load_pd( &spinor[2][2][0] );


  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1); 

  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  LOAD(hs1, 1);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  LOAD(hs2, 0); 

  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1);

  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);


  _mm_store_pd(&spinor[2][0][0], xmm0);
  _mm_store_pd(&spinor[2][1][0], xmm1);
  _mm_store_pd(&spinor[2][2][0], xmm2);

  /* Component 3 */
  LOAD(hs0, 0);

  xmm0 = _mm_load_pd( &spinor[3][0][0] );
  xmm1 = _mm_load_pd( &spinor[3][1][0] );
  xmm2 = _mm_load_pd( &spinor[3][2][0] );

  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1);


  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);


  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  LOAD(hs1, 0);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  LOAD(hs2, 1);


  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1);

  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);


  _mm_store_pd(&spinor[3][0][0], xmm0);
  _mm_store_pd(&spinor[3][1][0], xmm1);
  _mm_store_pd(&spinor[3][2][0], xmm2);
 


}

#undef SLOPPY_REGS
#undef LOAD

#ifdef __cplusplus
};
#endif

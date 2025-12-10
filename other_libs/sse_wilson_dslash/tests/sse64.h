#include <sse_config.h>

/*******************************************************************************
*
* File sse64.h
*
* Macros for Dirac spinors, SU(3) vectors and SU(3) matrices using
* inline assembly SSE and SSE2 instructions
*
* Needs gcc version 2.95.2 or later, and binutils snapshot 010122 or later
* if the SSE2 instructions are used
*
* Version: 1.1
* Author: Martin Luescher <luscher@mail.desy.de>
* Modified by Chris McClendon <cmcclend@jlab.org> for LHPC Collaboration
* Date: 17.03.2001
*
*******************************************************************************/

#ifdef __cplusplus
extern "C" { 
#endif

typedef struct
{
   int c1,c2,c3,c4;
} sse_int __attribute__ ((aligned (16)));

typedef struct
{
   float c1,c2,c3,c4;
} sse_float __attribute__ ((aligned (16)));

typedef struct
{
   double c1,c2;
} sse_double __attribute__ ((aligned (16)));


/*******************************************************************************
*
* Cache manipulation macros
*
*******************************************************************************/


#define mfence() \
__asm__ __volatile__ ("mfence" \
                      : )


#define _prefetch_single(addr) \
__asm__ __volatile__ ("prefetcht0 %0 \n\t" \
                       : \
                       : \
                       "m" (*(((char*)(((unsigned long)(addr))&~0x7f)))))

#define _prefetch_spinor(addr) \
__asm__ __volatile__ ("prefetcht0 %0 \n\t" \
                      "prefetcht0 %1" \
                      : \
                      : \
                      "m" (*(((char*)(((unsigned long)(addr))&~0x7f)))), \
                      "m" (*(((char*)(((unsigned long)(addr))&~0x7f))+128)))

#define _prefetch_nta_spinor(addr) \
__asm__ __volatile__ ("prefetchnta %0 \n\t" \
                      "prefetchnta %1" \
                      : \
                      : \
                      "m" (*(((char*)(((unsigned long)(addr))&~0x7f)))), \
                      "m" (*(((char*)(((unsigned long)(addr))&~0x7f))+128)))

#define _prefetch_su3(addr) \
__asm__ __volatile__ ("prefetcht0 %0 \n\t" \
                      "prefetcht0 %1" \
                      : \
                      : \
                      "m" (*(((char*)(((unsigned long)(addr))&~0x7f)))), \
                      "m" (*(((char*)(((unsigned long)(addr))&~0x7f))+128)))



static sse_int _sse_sgn __attribute__ ((unused)) ={0x0,0x80000000,0x0,0x0};
static sse_double _minus_one __attribute__ ((unused)) = {-1.0, -1.0};
static sse_double _conj      __attribute__ ((unused)) = {1.0, -1.0};
/*******************************************************************************
*
* Macros for su3 vectors used in D_psi version 2.0
*
* Most of these macros operate on su3 vectors that are stored
* in  xmm0,xmm1,xmm2 or xmm3,xmm4,xmm5. For example,
*
* xmm0 -> s.c1.re,s.c1.im
* xmm1 -> s.c2.re,s.c2.im
* xmm2 -> s.c3.re,s.c3.im
*
* where s is of type su3_vector
*
*******************************************************************************/

/*
* Loads an su3 std::vector s to xmm0,xmm1,xmm2
*/

#define _sse_load(s) \
__asm__ __volatile__ ("movapd %0, %%xmm0 \n\t" \
                      "movapd %1, %%xmm1 \n\t" \
                      "movapd %2, %%xmm2" \
                      : \
                      : \
                      "m" ((s)[0][0]), \
                      "m" ((s)[1][0]), \
                      "m" ((s)[2][0]))

     /*
      * Loads an su3 std::vector s to xmm3,xmm4,xmm5
      */

#define _sse_load_up(s) \
__asm__ __volatile__ ("movapd %0, %%xmm3 \n\t" \
                      "movapd %1, %%xmm4 \n\t" \
                      "movapd %2, %%xmm5" \
                      : \
                      : \
                      "m" ((s)[0][0]), \
                      "m" ((s)[1][0]), \
                      "m" ((s)[2][0]))



/*
* Stores xmm0,xmm1,xmm2 to the components r.c1,r.c2,r.c3 of an su3 std::vector
*/


#define _sse_store(r) \
__asm__ __volatile__ ("movapd %%xmm0, %0 \n\t" \
                      "movapd %%xmm1, %1 \n\t" \
                      "movapd %%xmm2, %2" \
                      : \
                      "=m" ((r)[0]), \
                      "=m" ((r)[1]), \
                      "=m" ((r)[2]))

/*
* Stores xmm3,xmm4,xmm5 to the components r.c1,r.c2,r.c3 of an su3 std::vector
*/

#define _sse_store_up(r) \
__asm__ __volatile__ ("movapd %%xmm3, %0 \n\t" \
                      "movapd %%xmm4, %1 \n\t" \
                      "movapd %%xmm5, %2" \
                      : \
                      "=m" ((r)[0]), \
                      "=m" ((r)[1]), \
                      "=m" ((r)[2]))





/*
* Multiplies xmm0,xmm1,xmm2 with a constant sse_double c
*/

#define _sse_vector_mul(c) \
__asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t" \
                      "mulpd %0, %%xmm1 \n\t" \
                      "mulpd %0, %%xmm2" \
                      : \
                      : \
                      "m" (c))

#define _sse_vector_mul_up(c) \
__asm__ __volatile__ ("mulpd %0, %%xmm3 \n\t" \
                      "mulpd %0, %%xmm4 \n\t" \
                      "mulpd %0, %%xmm5" \
                      : \
                      : \
                      "m" (c))

/*
* Adds xmm3,xmm4,xmm5 to xmm1,xmm2,xmm3
*/

#define _sse_vector_add() \
__asm__ __volatile__ ("addpd %%xmm3, %%xmm0 \n\t" \
                      "addpd %%xmm4, %%xmm1 \n\t" \
                      "addpd %%xmm5, %%xmm2" \
                      : \
                      :)


/*
* Subtracts xmm3,xmm4,xmm5 from xmm1,xmm2,xmm3
*/

#define _sse_vector_sub() \
__asm__ __volatile__ ("subpd %%xmm3, %%xmm0 \n\t" \
                      "subpd %%xmm4, %%xmm1 \n\t" \
                      "subpd %%xmm5, %%xmm2" \
                      : \
                      :)

/*
* Multiplies xmm3,xmm4,xmm5 with i
*/
/*note for compatibility these next two are no different */
#define _sse_vector_i_mul() \
__asm__ __volatile__ ("shufpd $0x1, %%xmm3, %%xmm3 \n\t" \
                      "shufpd $0x1, %%xmm4, %%xmm4 \n\t" \
                      "shufpd $0x1, %%xmm5, %%xmm5 \n\t" \
                      "xorpd %0, %%xmm3 \n\t" \
                      "xorpd %0, %%xmm4 \n\t" \
                      "xorpd %0, %%xmm5" \
                      : \
                      : \
                      "m" (_sse_sgn))


#define _sse_vector_i_mul_up() \
__asm__ __volatile__ ("shufpd $0x1, %%xmm3, %%xmm3 \n\t" \
                      "shufpd $0x1, %%xmm4, %%xmm4 \n\t" \
                      "shufpd $0x1, %%xmm5, %%xmm5 \n\t" \
                      "xorpd %0, %%xmm3 \n\t" \
                      "xorpd %0, %%xmm4 \n\t" \
                      "xorpd %0, %%xmm5" \
                      : \
                      : \
                      "m" (_sse_sgn))
/*needs to be fixed the neg part */
#define _sse_vector_i_mul_neg_up() \
__asm__ __volatile__ ("shufpd $0x1, %%xmm3, %%xmm3 \n\t" \
                      "shufpd $0x1, %%xmm4, %%xmm4 \n\t" \
                      "shufpd $0x1, %%xmm5, %%xmm5 \n\t" \
                      "mulpd %0, %%xmm3 \n\t" \
                      "mulpd %0, %%xmm4 \n\t" \
                      "mulpd %0, %%xmm5" \
                      : \
                      : \
                      "m" (_conj))


/*
* Multiplies an su3 std::vector s with an su3 matrix u, assuming s is
* stored in  xmm0,xmm1,xmm2
*
* On output the result is in xmm3,xmm4,xmm5 and the registers 
* xmm0,xmm1,xmm2 are changed
*/

/* first some macros */


/* macro overlay time */

#define COL_ROW_ORDER


/*use array indexing for SZIN */

#define _c11 [0][0][0]  /*note for asm instructions we have to get all the way to the data block */
#define _c12 [1][0][0]
#define _c13 [2][0][0]
#define _c21 [0][1][0]
#define _c22 [1][1][0]
#define _c23 [2][1][0]
#define _c31 [0][2][0]
#define _c32 [1][2][0]
#define _c33 [2][2][0]

#define _c11re [0][0][0]
#define _c12re [1][0][0]
#define _c13re [2][0][0]
#define _c21re [0][1][0]
#define _c22re [1][1][0]
#define _c23re [2][1][0]
#define _c31re [0][2][0]
#define _c32re [1][2][0]
#define _c33re [2][2][0]

#define _c11im [0][0][1]
#define _c12im [1][0][1]
#define _c13im [2][0][1]
#define _c21im [0][1][1]
#define _c22im [1][1][1]
#define _c23im [2][1][1]
#define _c31im [0][2][1]
#define _c32im [1][2][1]
#define _c33im [2][2][1]




#define _sse_su3_multiply(u) \
__asm__ __volatile__ ("movsd %0, %%xmm3 \n\t" \
                      "movsd %1, %%xmm6 \n\t" \
                      "movsd %2, %%xmm4 \n\t" \
                      "movsd %3, %%xmm7 \n\t" \
                      "movsd %4, %%xmm5 \n\t" \
                      "unpcklpd %%xmm3, %%xmm3 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm4, %%xmm4 \n\t" \
                      "mulpd %%xmm0, %%xmm3 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "unpcklpd %%xmm5, %%xmm5 \n\t" \
                      "mulpd %%xmm0, %%xmm4 \n\t" \
                      "addpd %%xmm6, %%xmm3 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "mulpd %%xmm0, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "movsd %5, %%xmm6 \n\t" \
                      "movsd %6, %%xmm7 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm3 \n\t" \
                      "movsd %7, %%xmm6 \n\t" \
                      "movsd %8, %%xmm7 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm4 \n\t" \
                      "addpd %%xmm7, %%xmm5" \
                      : \
                      : \
                      "m" ((u)_c11re), \
                      "m" ((u)_c12re), \
                      "m" ((u)_c21re), \
                      "m" ((u)_c23re), \
                      "m" ((u)_c31re), \
                      "m" ((u)_c32re), \
                      "m" ((u)_c13re), \
                      "m" ((u)_c22re), \
                      "m" ((u)_c33re)); \
__asm__ __volatile__ ("movsd %0, %%xmm6 \n\t" \
                      "movsd %1, %%xmm7 \n\t" \
                      "shufpd $0x1, %%xmm0, %%xmm0 \n\t" \
                      "shufpd $0x1, %%xmm1, %%xmm1 \n\t" \
                      "shufpd $0x1, %%xmm2, %%xmm2 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "xorpd %9, %%xmm0 \n\t" \
                      "xorpd %9, %%xmm1 \n\t" \
                      "xorpd %9, %%xmm2 \n\t" \
                      "mulpd %%xmm0, %%xmm6 \n\t" \
                      "mulpd %%xmm1, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm3 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "movsd %2, %%xmm6 \n\t" \
                      "movsd %3, %%xmm7 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm2, %%xmm6 \n\t" \
                      "mulpd %%xmm0, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "movsd %4, %%xmm6 \n\t" \
                      "movsd %5, %%xmm7 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm0, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm3 \n\t" \
                      "addpd %%xmm7, %%xmm5 \n\t" \
                      "movsd %6, %%xmm0 \n\t" \
                      "movsd %7, %%xmm6 \n\t" \
                      "movsd %8, %%xmm7 \n\t" \
                      "unpcklpd %%xmm0, %%xmm0 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm2, %%xmm0 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "addpd %%xmm0, %%xmm3 \n\t" \
                      "addpd %%xmm6, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm4" \
                      : \
                      : \
                      "m" ((u)_c11im), \
                      "m" ((u)_c22im), \
                      "m" ((u)_c33im), \
                      "m" ((u)_c21im), \
                      "m" ((u)_c12im), \
                      "m" ((u)_c31im), \
                      "m" ((u)_c13im), \
                      "m" ((u)_c32im), \
                      "m" ((u)_c23im), \
                      "m" (_sse_sgn))

/*
* Multiplies an su3 std::vector s with an su3 matrix u^dagger, assuming s is
* stored in  xmm0,xmm1,xmm2
*
* On output the result is in xmm3,xmm4,xmm5 and the registers 
* xmm0,xmm1,xmm2 are changed
*/

#define _sse_su3_inverse_multiply(u) \
__asm__ __volatile__ ("movsd %0, %%xmm3 \n\t" \
                      "movsd %1, %%xmm6 \n\t" \
                      "movsd %2, %%xmm4 \n\t" \
                      "movsd %3, %%xmm7 \n\t" \
                      "movsd %4, %%xmm5 \n\t" \
                      "unpcklpd %%xmm3, %%xmm3 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm4, %%xmm4 \n\t" \
                      "mulpd %%xmm0, %%xmm3 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "unpcklpd %%xmm5, %%xmm5 \n\t" \
                      "mulpd %%xmm0, %%xmm4 \n\t" \
                      "addpd %%xmm6, %%xmm3 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "mulpd %%xmm0, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "movsd %5, %%xmm6 \n\t" \
                      "movsd %6, %%xmm7 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm3 \n\t" \
                      "movsd %7, %%xmm6 \n\t" \
                      "movsd %8, %%xmm7 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm4 \n\t" \
                      "addpd %%xmm7, %%xmm5" \
                      : \
                      : \
                      "m" ((u)_c11re), \
                      "m" ((u)_c21re), \
                      "m" ((u)_c12re), \
                      "m" ((u)_c32re), \
                      "m" ((u)_c13re), \
                      "m" ((u)_c23re), \
                      "m" ((u)_c31re), \
                      "m" ((u)_c22re), \
                      "m" ((u)_c33re)); \
__asm__ __volatile__ ("movsd %0, %%xmm6 \n\t" \
                      "movsd %1, %%xmm7 \n\t" \
                      "xorpd %9, %%xmm0 \n\t" \
                      "xorpd %9, %%xmm1 \n\t" \
                      "xorpd %9, %%xmm2 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "shufpd $0x1, %%xmm0, %%xmm0 \n\t" \
                      "shufpd $0x1, %%xmm1, %%xmm1 \n\t" \
                      "shufpd $0x1, %%xmm2, %%xmm2 \n\t" \
                      "mulpd %%xmm0, %%xmm6 \n\t" \
                      "mulpd %%xmm1, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm3 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "movsd %2, %%xmm6 \n\t" \
                      "movsd %3, %%xmm7 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm2, %%xmm6 \n\t" \
                      "mulpd %%xmm0, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "movsd %4, %%xmm6 \n\t" \
                      "movsd %5, %%xmm7 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm0, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm3 \n\t" \
                      "addpd %%xmm7, %%xmm5 \n\t" \
                      "movsd %6, %%xmm0 \n\t" \
                      "movsd %7, %%xmm6 \n\t" \
                      "movsd %8, %%xmm7 \n\t" \
                      "unpcklpd %%xmm0, %%xmm0 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm2, %%xmm0 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "addpd %%xmm0, %%xmm3 \n\t" \
                      "addpd %%xmm6, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm4" \
                      : \
                      : \
                      "m" ((u)_c11im), \
                      "m" ((u)_c22im), \
                      "m" ((u)_c33im), \
                      "m" ((u)_c12im), \
                      "m" ((u)_c21im), \
                      "m" ((u)_c13im), \
                      "m" ((u)_c31im), \
                      "m" ((u)_c23im), \
                      "m" ((u)_c32im), \
                      "m" (_sse_sgn));





#ifdef __cplusplus
}
#endif

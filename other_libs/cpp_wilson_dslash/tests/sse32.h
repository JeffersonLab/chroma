/*******************************************************************************
* $Id: sse32.h,v 1.1 2008-08-26 13:24:29 bjoo Exp $
*
* File sse32.h
*
* Macros for Dirac spinors, SU(3) vectors and SU(3) matrices using
* inline assembly SSE and SSE2 instructions
*
* Needs gcc version 2.95.2 or later, and binutils snapshot 010122 or later
* if the SSE2 instructions are used
*
* Version: 2.1
* Author: Martin Luescher <luscher@mail.desy.de>
* Date: 15.03.2001
*
*******************************************************************************/

/* I have adapted this so that it is equivalent to old -DSSE2 -DSZIN -DP4
   which is how we always used it... */

#ifdef __cplusplus
extern "C" {
#endif

  /*! sse_float: A struct of 4 floats */
  typedef struct {
    float c1,c2,c3,c4;
  } sse_float __attribute__ ((aligned (16)));

  /*! sse_vector: A struct of 3 sse_floats */
  typedef struct {
    sse_float c1,c2,c3;
  } sse_vector __attribute__ ((aligned (16)));

  /*! sse_vector_wide sse_vector of 4 sse_floats */
  typedef struct { 
    sse_float c1,c2,c3,c4; 
  } sse_vector_wide __attribute__ ((aligned (16))); 

  /*! sign macros for shufs... */
  static sse_float _sse_internal_sgn12 __attribute__ ((unused)) ={-1.0f,-1.0f,1.0f,1.0f};
  static sse_float _sse_internal_sgn13 __attribute__ ((unused)) ={-1.0f,1.0f,-1.0f,1.0f};
  static sse_float _sse_internal_sgn14 __attribute__ ((unused)) ={-1.0f,1.0f,1.0f,-1.0f};
  static sse_float _sse_internal_sgn23 __attribute__ ((unused)) ={1.0f,-1.0f,-1.0f,1.0f};
  static sse_float _sse_internal_sgn24 __attribute__ ((unused)) ={1.0f,-1.0f,1.0f,-1.0f};
  static sse_float _sse_internal_sgn34 __attribute__ ((unused)) ={1.0f,1.0f,-1.0f,-1.0f};
  static sse_float _sse_internal_sgn1234 __attribute__ ((unused)) = {-1.0f,-1.0f,-1.0f,-1.0f};

  /*******************************************************************************
   *
   * Cache manipulation macros
   *
   *******************************************************************************/


#define _prefetch_spinor(addr) \
__asm__ __volatile__ ("prefetcht0 %0 \n\t" \
                      "prefetcht0 %1" \
                      : \
                      : \
                      "m" (*(((char*)(((unsigned long)(addr))&~0x7f)))), \
                      "m" (*(((char*)(((unsigned long)(addr))&~0x7f))+128)))


#define _prefetch_spinor_nta(addr) \
__asm__ __volatile__ ("prefetchnta %0 \n\t" \
                      "prefetchnta %1" \
                      : \
                      : \
                      "m" (*(((char*)(((unsigned long)(addr))&~0x7f)))), \
                      "m" (*(((char*)(((unsigned long)(addr))&~0x7f))+128)))


/*prefetches 4 directions for 2 sites' worth */
#define _prefetch_hspinor(addr) \
__asm__ __volatile__ ("prefetcht0 %0 \n\t" \
                      "prefetcht0 %1 \n\t" \
                      "prefetcht0 %2 \n\t" \
                      : \
		      : \
                      "m" (*(((char*)(((unsigned long)(addr))&~0x7f)))), \
                      "m" (*(((char*)(((unsigned long)(addr))&~0x7f))+128)), \
                      "m" (*(((char*)(((unsigned long)(addr))&~0x7f))+256)))
                     
#define _prefetch_single(addr) \
__asm__ __volatile__ ("prefetcht0 %0 \n\t" \
                      : \
                      : \
                      "m" (*(((char*)(((unsigned long)(addr))&~0x7f)))))

#define _prefetch_single_nta(addr) \
__asm__ __volatile__ ("prefetchnta %0 \n\t" \
                      : \
                      : \
                      "m" (*(((char*)(((unsigned long)(addr))&~0x7f)))))

 


#define _prefetch_su3(addr) \
__asm__ __volatile__ ("prefetcht0 %0 \n\t" \
                      "prefetcht0 %1" \
                      : \
                      : \
                      "m" (*(((char*)(((unsigned long)(addr))&~0x7f)))), \
                      "m" (*(((char*)(((unsigned long)(addr))&~0x7f))+128)))



/*******************************************************************************
*
* Macros for su3 vectors used in D_psi version 2.1
*
* Most of these macros operate on pairs of su3 vectors that are stored
* in the low and high words of xmm0,xmm1,xmm2 or xmm3,xmm4,xmm5. For example,
*
* xmm0 -> sl.c1.re,sl.c1.im,sh.c1.re,sh.c1.im
* xmm1 -> sl.c2.re,sl.c2.im,sh.c2.re,sh.c2.im
* xmm2 -> sl.c3.re,sl.c3.im,sh.c3.re,sh.c3.im
*
* (where sl and sh are of type su3_vector). This can also be interpreted as
* an sse_vector s that is stored in these registers according to
*
* xmm0 -> s.c1.c1,s.c1.c2,s.c1.c3,s.c1.c4
* xmm1 -> s.c2.c1,s.c2.c2,s.c2.c3,s.c2.c4
* xmm2 -> s.c3.c1,s.c3.c2,s.c3.c3,s.c3.c4
*
* The load and store macros can be used to move data in either format
* from and to the xmm registers
*
*******************************************************************************/

/*
* Loads two su3 vectors sl and sh to the low and high words of xmm0,xmm1,xmm2
*/

#define mfence() \
__asm__ __volatile__ ("mfence" \
                      : )

/* use array indices instead of structures for SZIN*/
/* note we must use enough brackets to hit the physical memory and not a pointer to the physical memory */
#define _sse_pair_load(sl,sh) \
__asm__ __volatile__ ("movsd %0, %%xmm0 \n\t" \
                      "movsd %1, %%xmm1 \n\t" \
                      "movsd %2, %%xmm2 \n\t" \
                      "movhps %3, %%xmm0 \n\t" \
                      "movhps %4, %%xmm1 \n\t" \
                      "movhps %5, %%xmm2" \
                      : \
                      : \
                      "m" ((sl)[0][0]), \
                      "m" ((sl)[1][0]), \
                      "m" ((sl)[2][0]), \
                      "m" ((sh)[0][0]), \
                      "m" ((sh)[1][0]), \
                      "m" ((sh)[2][0]))


#define _sse_pair_load_to_xmm(sh,sh2) \
__asm__ __volatile__ ("movhps %0, %%xmm0 \n\t" \
                      "movhps %1, %%xmm1 \n\t" \
                      "movhps %2, %%xmm2 \n\t" \
                      "movhps %3, %%xmm3 \n\t" \
                      "movhps %4, %%xmm4 \n\t" \
                      "movhps %5, %%xmm5" \
                      : \
                      : \
                      "m" ((sh).c1), \
                      "m" ((sh).c2), \
                      "m" ((sh).c3), \
                      "m" ((sh2).c1), \
                      "m" ((sh2).c2), \
                      "m" ((sh2).c3))

#define _sse_clear_mmx() \
__asm__ __volatile__ ("emms" \
                      : )
      
    


#define _sse_pair_load_to_mmx(sl,sl2) \
__asm__ __volatile__ ("movq %0, %%mm2 \n\t" \
                      "movq %1, %%mm3 \n\t" \
                      "movq %2, %%mm4 \n\t" \
                      "movq %3, %%mm5 \n\t" \
                      "movq %4, %%mm6 \n\t" \
                      "movq %5, %%mm7" \
                      : \
                      : \
                      "m" ((sl).c1), \
                      "m" ((sl).c2), \
                      "m" ((sl).c3), \
                      "m" ((sl2).c1), \
                      "m" ((sl2).c2), \
                      "m" ((sl2).c3))



#define _sse_pair_load_from_mmx(sl,sh) \
__asm__ __volatile__ ("movdq2q %%mm2, %%xmm0 \n\t" \
                      "movdq2q %%mm3, %%xmm1 \n\t" \
                      "movdq2q %%mm4, %%xmm2 \n\t" \
                      "movdq2q %%mm5, %%xmm3 \n\t" \
                      "movdq2q %%mm6, %%xmm4 \n\t" \
                      "movdq2q %%mm7, %%xmm5" \
                      : \
                      : )


/*
* Loads two su3 vectors sl and sh to the low and high words of xmm3,xmm4,xmm5
*/  




/* use array indexes for SZIN */
#define _sse_pair_load_up(sl,sh) \
__asm__ __volatile__ ("movsd %0, %%xmm3 \n\t" \
                      "movsd %1, %%xmm4 \n\t" \
                      "movsd %2, %%xmm5 \n\t" \
                      "movhps %3, %%xmm3 \n\t" \
                      "movhps %4, %%xmm4 \n\t" \
                      "movhps %5, %%xmm5" \
                      : \
                      : \
                      "m" ((sl)[0][0]), \
                      "m" ((sl)[1][0]), \
                      "m" ((sl)[2][0]), \
                      "m" ((sh)[0][0]), \
                      "m" ((sh)[1][0]), \
                      "m" ((sh)[2][0]))



/*
* Stores the low and high words of xmm0,xmm1,xmm2 to the su3 vectors rl and rh
*/

/* use array indexing for SZIN */
#define _sse_pair_store(rl,rh) \
__asm__ __volatile__ ("movlps %%xmm0, %0 \n\t" \
                      "movlps %%xmm1, %1 \n\t" \
                      "movlps %%xmm2, %2 \n\t" \
                      "movhps %%xmm0, %3 \n\t" \
                      "movhps %%xmm1, %4 \n\t" \
                      "movhps %%xmm2, %5" \
                      : \
                      "=m" ((rl)[0]), \
                      "=m" ((rl)[1]), \
                      "=m" ((rl)[2]), \
                      "=m" ((rh)[0]), \
                      "=m" ((rh)[1]), \
                      "=m" ((rh)[2]))

#define _sse_vector_store_mmx(s) \
__asm__ __volatile__ ("movlps %%xmm0, %%mm0 \n\t" \
                      "movlps %%xmm1, %%mm1 \n\t" \
                      "movlps %%xmm2, %%mm2 \n\t" \
                      "movhps %%xmm0, %%mm3 \n\t" \
                      "movhps %%xmm1, %%mm4 \n\t" \
                      "movhps %%xmm2, %%mm5 " \
                      : \
                      : )

#define _sse_vector_load_mmx(s) \
__asm__ __volatile__ ("movlps %%mm0, %%xmm0 \n\t" \
                      "movlps %%mm1, %%xmm1 \n\t" \
                      "movlps %%mm2, %%xmm2 \n\t" \
                      "movhps %%mm0, %%xmm3 \n\t" \
                      "movhps %%mm1, %%xmm4 \n\t" \
                      "movhps %%mm2, %%xmm5 " \
                      : \
                      : )

/*
* Stores the low and high words of xmm3,xmm4,xmm5 to the su3 vectors rl and rh
*/

/* use array indexing for SZIN */

#define _sse_pair_store_up(rl,rh) \
__asm__ __volatile__ ("movlps %%xmm3, %0 \n\t" \
                      "movlps %%xmm4, %1 \n\t" \
                      "movlps %%xmm5, %2 \n\t" \
                      "movhps %%xmm3, %3 \n\t" \
                      "movhps %%xmm4, %4 \n\t" \
                      "movhps %%xmm5, %5" \
                      : \
                      "=m" ((rl)[0]), \
                      "=m" ((rl)[1]), \
                      "=m" ((rl)[2]), \
                      "=m" ((rh)[0]), \
                      "=m" ((rh)[1]), \
                      "=m" ((rh)[2]))



/*
* Loads the components s.c1,s.c2,s.c3 of an _sse_vector s to xmm0,xmm1,xmm2
*/



#define _sse_vector_load(s) \
__asm__ __volatile__ ("movaps %0, %%xmm0 \n\t" \
                      "movaps %1, %%xmm1 \n\t" \
                      "movaps %2, %%xmm2" \
                      : \
                      : \
                      "m" ((s)[0][0][0]), \
                      "m" ((s)[1][0][0]), \
                      "m" ((s)[2][0][0]))




/*
* Loads the components s.c1,s.c2,s.c3 of an _sse_vector s to xmm3,xmm4,xmm5
*/
#define _sse_vector_load_up(s) \
__asm__ __volatile__ ("movaps %0, %%xmm3 \n\t" \
                      "movaps %1, %%xmm4 \n\t" \
                      "movaps %2, %%xmm5" \
                      : \
                      : \
                      "m" ((s)[0][0][0]), \
                      "m" ((s)[1][0][0]), \
                      "m" ((s)[2][0][0]))




/*
* Stores xmm0,xmm1,xmm2 to the components r.c1,r.c2,r.c3 of an _sse_vector r 
*/


#define _sse_vector_store(r) \
__asm__ __volatile__ ("movaps %%xmm0, %0 \n\t" \
                      "movaps %%xmm1, %1 \n\t" \
                      "movaps %%xmm2, %2" \
                      : \
                      "=m" ((r)[0]), \
                      "=m" ((r)[1]), \
                      "=m" ((r)[2]))




/* movntps might only work on p4 */

#define _sse_vector_store_nta(r) \
__asm__ __volatile__ ("movntps %%xmm0, %0 \n\t" \
                      "movntps %%xmm1, %1 \n\t" \
                      "movntps %%xmm2, %2" \
                      : \
                      "=m" ((r)[0]), \
                      "=m" ((r)[1]), \
                      "=m" ((r)[2]))


/*
* Stores xmm3,xmm4,xmm5 to the components r.c1,r.c2,r.c3 of an _sse_vector r 
*/


#define _sse_vector_store_up(r) \
__asm__ __volatile__ ("movaps %%xmm3, %0 \n\t" \
                      "movaps %%xmm4, %1 \n\t" \
                      "movaps %%xmm5, %2" \
                      : \
                      "=m" ((r)[0]), \
                      "=m" ((r)[1]), \
                      "=m" ((r)[2]))





#define _sse_vector_store_up_nta(r) \
__asm__ __volatile__ ("movntps %%xmm3, %0 \n\t" \
                      "movntps %%xmm4, %1 \n\t" \
                      "movntps %%xmm5, %2" \
                      : \
                      "=m" ((r)[0]), \
                      "=m" ((r)[1]), \
                      "=m" ((r)[2]))




/*
* Multiplies xmm0,xmm1,xmm2 with a constant sse_float c
*/

#define _sse_vector_mul(c) \
__asm__ __volatile__ ("mulps %0, %%xmm0 \n\t" \
                      "mulps %0, %%xmm1 \n\t" \
                      "mulps %0, %%xmm2" \
                      : \
                      : \
                      "m" (c))


/*
* Multiplies xmm3,xmm4,xmm5 with a constant sse_float c
*/

#define _sse_vector_mul_up(c) \
__asm__ __volatile__ ("mulps %0, %%xmm3 \n\t" \
                      "mulps %0, %%xmm4 \n\t" \
                      "mulps %0, %%xmm5" \
                      : \
                      : \
                      "m" (c))

#define _sse_vector_i_mul_plusminus_up() \
__asm__ __volatile__ ("shufps $0xb1, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0xb1, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0xb1, %%xmm5, %%xmm5 \n\t" \
                      "mulps %0, %%xmm3 \n\t" \
                      "mulps %0, %%xmm4 \n\t" \
                      "mulps %0, %%xmm5 \n\t" \
                       : \
                       : \
                       "m" (_sse_internal_sgn23));

/*
* Adds xmm3,xmm4,xmm5 to xmm1,xmm2,xmm3
*/

#define _sse_vector_add() \
__asm__ __volatile__ ("addps %%xmm3, %%xmm0 \n\t" \
                      "addps %%xmm4, %%xmm1 \n\t" \
                      "addps %%xmm5, %%xmm2" \
                      : \
                      :)


/*
* Subtracts xmm3,xmm4,xmm5 from xmm1,xmm2,xmm3
*/

#define _sse_vector_sub() \
__asm__ __volatile__ ("subps %%xmm3, %%xmm0 \n\t" \
                      "subps %%xmm4, %%xmm1 \n\t" \
                      "subps %%xmm5, %%xmm2" \
                      : \
                      :)

/*
* Multiplies the high words xmm3,xmm4,xmm5 with -1 and adds these registers
* to xmm0,xmm1,xmm2
*/

#define _sse_vector_addsub() \
__asm__ __volatile__ ("mulps %0, %%xmm3 \n\t" \
                      "mulps %0, %%xmm4 \n\t" \
                      "mulps %0, %%xmm5 \n\t" \
                      "addps %%xmm3, %%xmm0 \n\t" \
                      "addps %%xmm4, %%xmm1 \n\t" \
                      "addps %%xmm5, %%xmm2" \
                      : \
                      : \
                      "m" (_sse_internal_sgn34))

/*
* Multiplies the low words xmm3,xmm4,xmm5 with -1 and adds these registers
* to xmm0,xmm1,xmm2
*/

#define _sse_vector_subadd() \
__asm__ __volatile__ ("mulps %0, %%xmm3 \n\t" \
                      "mulps %0, %%xmm4 \n\t" \
                      "mulps %0, %%xmm5 \n\t" \
                      "addps %%xmm3, %%xmm0 \n\t" \
                      "addps %%xmm4, %%xmm1 \n\t" \
                      "addps %%xmm5, %%xmm2" \
                      : \
                      : \
                      "m" (_sse_internal_sgn12))


/*
* Multiplies xmm3,xmm4,xmm5 with i and adds them to xmm1,xmm2,xmm3
*/

#define _sse_vector_i_add() \
__asm__ __volatile__ ("shufps $0xb1, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0xb1, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0xb1, %%xmm5, %%xmm5 \n\t" \
                      "mulps %0, %%xmm3 \n\t" \
                      "mulps %0, %%xmm4 \n\t" \
                      "mulps %0, %%xmm5 \n\t" \
                      "addps %%xmm3, %%xmm0 \n\t" \
                      "addps %%xmm4, %%xmm1 \n\t" \
                      "addps %%xmm5, %%xmm2" \
                      : \
                      : \
                      "m" (_sse_internal_sgn13))

/*
* Multiplies xmm3,xmm4,xmm5 with i and subtracts them from xmm1,xmm2,xmm3
*/

#define _sse_vector_i_sub() \
__asm__ __volatile__ ("shufps $0xb1, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0xb1, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0xb1, %%xmm5, %%xmm5 \n\t" \
                      "mulps %0, %%xmm3 \n\t" \
                      "mulps %0, %%xmm4 \n\t" \
                      "mulps %0, %%xmm5 \n\t" \
                      "addps %%xmm3, %%xmm0 \n\t" \
                      "addps %%xmm4, %%xmm1 \n\t" \
                      "addps %%xmm5, %%xmm2" \
                      : \
                      : \
                      "m" (_sse_internal_sgn24))

/*
* Exchanges the high and low words of xmm3,xmm4,xmm5, multiplies them with i
* and adds the result to xmm1,xmm2,xmm3
*/

#define _sse_vector_xch_i_add() \
__asm__ __volatile__ ("shufps $0x1b, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0x1b, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0x1b, %%xmm5, %%xmm5 \n\t" \
                      "mulps %0, %%xmm3 \n\t" \
                      "mulps %0, %%xmm4 \n\t" \
                      "mulps %0, %%xmm5 \n\t" \
                      "addps %%xmm3, %%xmm0 \n\t" \
                      "addps %%xmm4, %%xmm1 \n\t" \
                      "addps %%xmm5, %%xmm2" \
                      : \
                      : \
                      "m" (_sse_internal_sgn13))

#define _sse_vector_xch_i_mul_up() \
__asm__ __volatile__ ("shufps $0x1b, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0x1b, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0x1b, %%xmm5, %%xmm5 \n\t" \
                      "mulps %0, %%xmm3 \n\t" \
                      "mulps %0, %%xmm4 \n\t" \
                      "mulps %0, %%xmm5 \n\t" \
                      : \
                      : \
                      "m" (_sse_internal_sgn13))

#define _sse_vector_xch_i_mul_neg_up() \
__asm__ __volatile__ ("shufps $0x1b, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0x1b, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0x1b, %%xmm5, %%xmm5 \n\t" \
                      "mulps %0, %%xmm3 \n\t" \
                      "mulps %0, %%xmm4 \n\t" \
                      "mulps %0, %%xmm5 \n\t" \
                      : \
                      : \
                      "m" (_sse_internal_sgn24))


/*
* Exchanges the high and low words of xmm3,xmm4,xmm5, multiplies them with i
* and subtracts the result from xmm1,xmm2,xmm3
*/

#define _sse_vector_xch_i_sub() \
__asm__ __volatile__ ("shufps $0x1b, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0x1b, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0x1b, %%xmm5, %%xmm5 \n\t" \
                      "mulps %0, %%xmm3 \n\t" \
                      "mulps %0, %%xmm4 \n\t" \
                      "mulps %0, %%xmm5 \n\t" \
                      "addps %%xmm3, %%xmm0 \n\t" \
                      "addps %%xmm4, %%xmm1 \n\t" \
                      "addps %%xmm5, %%xmm2" \
                      : \
                      : \
                      "m" (_sse_internal_sgn24))

/*
* Multiplies the low and high words of xmm3,xmm4,xmm5 with i and -i
* respectively and adds these registers to xmm1,xmm2,xmm3
*/

#define _sse_vector_i_addsub() \
__asm__ __volatile__ ("shufps $0xb1, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0xb1, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0xb1, %%xmm5, %%xmm5 \n\t" \
                      "mulps %0, %%xmm3 \n\t" \
                      "mulps %0, %%xmm4 \n\t" \
                      "mulps %0, %%xmm5 \n\t" \
                      "addps %%xmm3, %%xmm0 \n\t" \
                      "addps %%xmm4, %%xmm1 \n\t" \
                      "addps %%xmm5, %%xmm2" \
                      : \
                      : \
                      "m" (_sse_internal_sgn14))

/*
* Multiplies the low and high words of xmm3,xmm4,xmm5 with -i and i
* respectively and adds these registers to xmm1,xmm2,xmm3
*/

#define _sse_vector_i_subadd() \
__asm__ __volatile__ ("shufps $0xb1, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0xb1, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0xb1, %%xmm5, %%xmm5 \n\t" \
                      "mulps %0, %%xmm3 \n\t" \
                      "mulps %0, %%xmm4 \n\t" \
                      "mulps %0, %%xmm5 \n\t" \
                      "addps %%xmm3, %%xmm0 \n\t" \
                      "addps %%xmm4, %%xmm1 \n\t" \
                      "addps %%xmm5, %%xmm2" \
                      : \
                      : \
                      "m" (_sse_internal_sgn23))

/*
* Exchanges the high and low words in xmm3,xmm4,xmm5 
*/

#define _sse_vector_xch() \
__asm__ __volatile__ ("shufps $0x4e, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0x4e, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0x4e, %%xmm5, %%xmm5" \
                      : \
                      :)

/*
* Multiplies a pair sl,sh of su3 vectors with an su3 matrix u,
* assuming sl and sh are in the low and high words of xmm0,xmm1,xmm2
*
* On output the result is in xmm3,xmm4,xmm5 and the registers 
* xmm0,xmm1,xmm2 are changed
*/


/* use array indices for u for SZIN */
/* note SZIN has the different ordering */




#define _sse_su3_multiply(u) \
__asm__ __volatile__ ("movss %0, %%xmm3 \n\t" \
                      "movss %1, %%xmm6 \n\t" \
                      "movss %2, %%xmm4 \n\t" \
                      "movss %3, %%xmm7 \n\t" \
                      "movss %4, %%xmm5 \n\t" \
                      "shufps $0x0, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm4, %%xmm4 \n\t" \
                      "mulps %%xmm0, %%xmm3 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm1, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm5, %%xmm5 \n\t" \
                      "mulps %%xmm0, %%xmm4 \n\t" \
                      "addps %%xmm6, %%xmm3 \n\t" \
                      "mulps %%xmm2, %%xmm7 \n\t" \
                      "mulps %%xmm0, %%xmm5 \n\t" \
                      "addps %%xmm7, %%xmm4 \n\t" \
                      "movss %5, %%xmm6 \n\t" \
                      "movss %6, %%xmm7 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm1, %%xmm6 \n\t" \
                      "mulps %%xmm2, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm5 \n\t" \
                      "addps %%xmm7, %%xmm3 \n\t" \
                      "movss %7, %%xmm6 \n\t" \
                      "movss %8, %%xmm7 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm1, %%xmm6 \n\t" \
                      "mulps %%xmm2, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm4 \n\t" \
                      "addps %%xmm7, %%xmm5" \
                      : \
                      : \
                      "m" ((u)[0][0][0]), \
                      "m" ((u)[1][0][0]), \
                      "m" ((u)[0][1][0]), \
                      "m" ((u)[2][1][0]), \
                      "m" ((u)[0][2][0]), \
                      "m" ((u)[1][2][0]), \
                      "m" ((u)[2][0][0]), \
                      "m" ((u)[1][1][0]), \
                      "m" ((u)[2][2][0])); \
__asm__ __volatile__ ("movss %0, %%xmm6 \n\t" \
                      "movss %1, %%xmm7 \n\t" \
                      "shufps $0xb1, %%xmm0, %%xmm0 \n\t" \
                      "shufps $0xb1, %%xmm1, %%xmm1 \n\t" \
                      "shufps $0xb1, %%xmm2, %%xmm2 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %9, %%xmm0 \n\t" \
                      "mulps %9, %%xmm1 \n\t" \
                      "mulps %9, %%xmm2 \n\t" \
                      "mulps %%xmm0, %%xmm6 \n\t" \
                      "mulps %%xmm1, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm3 \n\t" \
                      "addps %%xmm7, %%xmm4 \n\t" \
                      "movss %2, %%xmm6 \n\t" \
                      "movss %3, %%xmm7 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm2, %%xmm6 \n\t" \
                      "mulps %%xmm0, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm5 \n\t" \
                      "addps %%xmm7, %%xmm4 \n\t" \
                      "movss %4, %%xmm6 \n\t" \
                      "movss %5, %%xmm7 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm1, %%xmm6 \n\t" \
                      "mulps %%xmm0, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm3 \n\t" \
                      "addps %%xmm7, %%xmm5 \n\t" \
                      "movss %6, %%xmm0 \n\t" \
                      "movss %7, %%xmm6 \n\t" \
                      "movss %8, %%xmm7 \n\t" \
                      "shufps $0x0, %%xmm0, %%xmm0 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm2, %%xmm0 \n\t" \
                      "mulps %%xmm1, %%xmm6 \n\t" \
                      "mulps %%xmm2, %%xmm7 \n\t" \
                      "addps %%xmm0, %%xmm3 \n\t" \
                      "addps %%xmm6, %%xmm5 \n\t" \
                      "addps %%xmm7, %%xmm4" \
                      : \
                      : \
                      "m" ((u)[0][0][1]), \
                      "m" ((u)[1][1][1]), \
                      "m" ((u)[2][2][1]), \
                      "m" ((u)[0][1][1]), \
                      "m" ((u)[1][0][1]), \
                      "m" ((u)[0][2][1]), \
                      "m" ((u)[2][0][1]), \
                      "m" ((u)[1][2][1]), \
                      "m" ((u)[2][1][1]), \
                      "m" (_sse_internal_sgn13))

/*
* Multiplies a pair sl,sh of su3 vectors with an su3 matrix u^dagger, 
* assuming sl and sh are in the low and high words of xmm0,xmm1,xmm2
*
* On output the result is in xmm3,xmm4,xmm5 and the registers 
* xmm0,xmm1,xmm2 are changed
*/

#define _sse_su3_inverse_multiply(u) \
__asm__ __volatile__ ("movss %0, %%xmm3 \n\t" \
                      "movss %1, %%xmm6 \n\t" \
                      "movss %2, %%xmm4 \n\t" \
                      "movss %3, %%xmm7 \n\t" \
                      "movss %4, %%xmm5 \n\t" \
                      "shufps $0x0, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm4, %%xmm4 \n\t" \
                      "mulps %%xmm0, %%xmm3 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm1, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm5, %%xmm5 \n\t" \
                      "mulps %%xmm0, %%xmm4 \n\t" \
                      "addps %%xmm6, %%xmm3 \n\t" \
                      "mulps %%xmm2, %%xmm7 \n\t" \
                      "mulps %%xmm0, %%xmm5 \n\t" \
                      "addps %%xmm7, %%xmm4 \n\t" \
                      "movss %5, %%xmm6 \n\t" \
                      "movss %6, %%xmm7 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm1, %%xmm6 \n\t" \
                      "mulps %%xmm2, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm5 \n\t" \
                      "addps %%xmm7, %%xmm3 \n\t" \
                      "movss %7, %%xmm6 \n\t" \
                      "movss %8, %%xmm7 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm1, %%xmm6 \n\t" \
                      "mulps %%xmm2, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm4 \n\t" \
                      "addps %%xmm7, %%xmm5" \
                      : \
                      : \
                      "m" ((u)[0][0][0]), \
                      "m" ((u)[0][1][0]), \
                      "m" ((u)[1][0][0]), \
                      "m" ((u)[1][2][0]), \
                      "m" ((u)[2][0][0]), \
                      "m" ((u)[2][1][0]), \
                      "m" ((u)[0][2][0]), \
                      "m" ((u)[1][1][0]), \
                      "m" ((u)[2][2][0])); \
__asm__ __volatile__ ("movss %0, %%xmm6 \n\t" \
                      "movss %1, %%xmm7 \n\t" \
                      "shufps $0xb1, %%xmm0, %%xmm0 \n\t" \
                      "shufps $0xb1, %%xmm1, %%xmm1 \n\t" \
                      "shufps $0xb1, %%xmm2, %%xmm2 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %9, %%xmm0 \n\t" \
                      "mulps %9, %%xmm1 \n\t" \
                      "mulps %9, %%xmm2 \n\t" \
                      "mulps %%xmm0, %%xmm6 \n\t" \
                      "mulps %%xmm1, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm3 \n\t" \
                      "addps %%xmm7, %%xmm4 \n\t" \
                      "movss %2, %%xmm6 \n\t" \
                      "movss %3, %%xmm7 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm2, %%xmm6 \n\t" \
                      "mulps %%xmm0, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm5 \n\t" \
                      "addps %%xmm7, %%xmm4 \n\t" \
                      "movss %4, %%xmm6 \n\t" \
                      "movss %5, %%xmm7 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm1, %%xmm6 \n\t" \
                      "mulps %%xmm0, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm3 \n\t" \
                      "addps %%xmm7, %%xmm5 \n\t" \
                      "movss %6, %%xmm0 \n\t" \
                      "movss %7, %%xmm6 \n\t" \
                      "movss %8, %%xmm7 \n\t" \
                      "shufps $0x0, %%xmm0, %%xmm0 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm2, %%xmm0 \n\t" \
                      "mulps %%xmm1, %%xmm6 \n\t" \
                      "mulps %%xmm2, %%xmm7 \n\t" \
                      "addps %%xmm0, %%xmm3 \n\t" \
                      "addps %%xmm6, %%xmm5 \n\t" \
                      "addps %%xmm7, %%xmm4" \
                      : \
                      : \
                      "m" ((u)[0][0][1]), \
                      "m" ((u)[1][1][1]), \
                      "m" ((u)[2][2][1]), \
                      "m" ((u)[1][0][1]), \
                      "m" ((u)[0][1][1]), \
                      "m" ((u)[2][0][1]), \
                      "m" ((u)[0][2][1]), \
                      "m" ((u)[2][1][1]), \
                      "m" ((u)[1][2][1]), \
                      "m" (_sse_internal_sgn24));


#ifdef __cplusplus
}
#endif

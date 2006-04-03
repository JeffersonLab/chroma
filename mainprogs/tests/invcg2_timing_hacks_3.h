// -*- C++ -*-
// $Id: invcg2_timing_hacks_3.h,v 3.0 2006-04-03 04:59:14 edwards Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#ifndef __invcg2_timing_hacks_2_h__
#define __invcg2_timing_hacks_2_h__

#include "chromabase.h"
#include "linearop.h"
#include "actions/ferm/linop/dslash_w.h"

//! Highly optimised Conjugate-Gradient (CGNE) algorithm for a Even Odd Preconditioned
                                                                                
// Perversly theser are the types used in our axpys.
typedef OLattice< PSpinVector< PColorVector< RComplex< PScalar<REAL> >, Nc>, Ns> > LFerm;
                                                                                
typedef OScalar< PScalar < PScalar < RScalar< PScalar < REAL > > > > > LScal;
typedef OScalar< PScalar < PScalar < RScalar< PScalar < DOUBLE > > > > > LDble;                                                                                
// Get at the REAL embedded in an LScal
#define AT_REAL(a)  (a.elem().elem().elem().elem().elem())
                                                                                
// Get the first element of a vector over a subset
#define FIRST_ELEM(a,s) (&(a.elem(s.start()).elem(0).elem(0).real().elem()))
                                                                               
void InvCG2EvenOddPrecWilsLinOpTHack(const WilsonDslash &D,
                                const LFerm& chi,
                                LFerm& psi,
                                const LScal& mass,
                                const LScal& RsdCG,
                                int MaxCG,
                                int& n_count);
 

// GNUC vector type
typedef float v4sf __attribute__((mode(V4SF),aligned(16)));

// vaxpy3 and norm put together
inline
void vaxpy3_norm(REAL *Out,REAL *scalep,REAL *InScale, REAL *Add,int n_3vec,
		 REAL* dsum)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "SSE_TEST: vaxpy3_norm" << endl;
#endif

  int n_loops = n_3vec;

  v4sf vscalep = __builtin_ia32_loadss(scalep);
  asm("shufps\t$0,%0,%0" : "+x" (vscalep));

  REAL fzero = 0.0;
  register v4sf vsum = __builtin_ia32_loadss(&fzero);
  asm("shufps\t$0,%0,%0" : "+x" (vsum));

  for (; n_loops-- > 0; )
  {
    register v4sf vtmp;

    vtmp = __builtin_ia32_addps(__builtin_ia32_mulps(vscalep, __builtin_ia32_loadaps(InScale+ 0)), __builtin_ia32_loadaps(Add+ 0));
    vsum = __builtin_ia32_addps(vsum, __builtin_ia32_mulps(vtmp, vtmp));
    __builtin_ia32_storeaps(Out+ 0, vtmp);

    vtmp = __builtin_ia32_addps(__builtin_ia32_mulps(vscalep, __builtin_ia32_loadaps(InScale+ 4)), __builtin_ia32_loadaps(Add+ 4));
    vsum = __builtin_ia32_addps(vsum, __builtin_ia32_mulps(vtmp, vtmp));
    __builtin_ia32_storeaps(Out+ 4, vtmp);

    vtmp = __builtin_ia32_addps(__builtin_ia32_mulps(vscalep, __builtin_ia32_loadaps(InScale+ 8)), __builtin_ia32_loadaps(Add+ 8));
    vsum = __builtin_ia32_addps(vsum, __builtin_ia32_mulps(vtmp, vtmp));
    __builtin_ia32_storeaps(Out+ 8, vtmp);

    vtmp = __builtin_ia32_addps(__builtin_ia32_mulps(vscalep, __builtin_ia32_loadaps(InScale+12)), __builtin_ia32_loadaps(Add+12));
    vsum = __builtin_ia32_addps(vsum, __builtin_ia32_mulps(vtmp, vtmp));
    __builtin_ia32_storeaps(Out+12, vtmp);

    vtmp = __builtin_ia32_addps(__builtin_ia32_mulps(vscalep, __builtin_ia32_loadaps(InScale+16)), __builtin_ia32_loadaps(Add+16));
    vsum = __builtin_ia32_addps(vsum, __builtin_ia32_mulps(vtmp, vtmp));
    __builtin_ia32_storeaps(Out+16, vtmp);

    vtmp = __builtin_ia32_addps(__builtin_ia32_mulps(vscalep, __builtin_ia32_loadaps(InScale+20)), __builtin_ia32_loadaps(Add+20));
    vsum = __builtin_ia32_addps(vsum, __builtin_ia32_mulps(vtmp, vtmp));
    __builtin_ia32_storeaps(Out+20, vtmp);

    Out += 24; InScale += 24; Add += 24;
  }

  REAL fsum[4];
  __builtin_ia32_storeaps(fsum, vsum);
  *dsum = (REAL)(fsum[0] + fsum[1] + fsum[2] + fsum[3]);
}



#endif

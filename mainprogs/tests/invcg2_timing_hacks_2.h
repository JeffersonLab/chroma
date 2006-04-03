// -*- C++ -*-
// $Id: invcg2_timing_hacks_2.h,v 3.0 2006-04-03 04:59:14 edwards Exp $
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
 

#endif

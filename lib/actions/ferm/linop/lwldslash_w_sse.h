// -*- C++ -*-
// $Id: lwldslash_w_sse.h,v 1.1 2003-09-10 18:15:05 bjoo Exp $
/*! \file
 *  \brief Wilson Dslash linear operator
 */

#ifndef __lwldslash_sse_h__
#define __lwldslash_sse_h__

#include "linearop.h"

#include <sse_config.h>

using namespace QDP;

//! General Wilson-Dirac dslash
/*!
 * \ingroup linop
 *
 * DSLASH
 *
 * This routine is specific to Wilson fermions!
 *
 * Description:
 *
 * This routine applies the operator D' to Psi, putting the result in Chi.
 *
 *	       Nd-1
 *	       ---
 *	       \
 *   chi(x)  :=  >  U  (x) (1 - isign gamma  ) psi(x+mu)
 *	       /    mu			  mu
 *	       ---
 *	       mu=0
 *
 *	             Nd-1
 *	             ---
 *	             \    +
 *                +    >  U  (x-mu) (1 + isign gamma  ) psi(x-mu)
 *	             /    mu			   mu
 *	             ---
 *	             mu=0
 *
 */

#if SSE_PRECISION == 32
typedef float SSEREAL;
#warning "Building Packer for 32 Bits"
#elif SSE_PRECISION;
typedef double SSEREAL double;
#warning "Building Packer for 64 Bits"
#else
#error "Precision Not supported, define SSE_PRECISION"
#endif
                                                                                
typedef SSEREAL u_mat_array[3][3][2];

class SSEWilsonDslash : public DslashLinearOperator
{
public:
  //! Empty constructor. Must use create later
  SSEWilsonDslash() {}

  //! Full constructor
  SSEWilsonDslash(const multi1d<LatticeColorMatrix>& _u) {create(_u);}

  //! Creation routine
  void create(const multi1d<LatticeColorMatrix>& _u);

  //! No real need for cleanup here
  ~SSEWilsonDslash();

  /**
   * Apply a dslash
   *
   * \param psi     source                                      (Read)
   * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
   * \param cb	    Checkerboard of OUTPUT vector               (Read) 
   *
   * \return The output of applying dslash on dslash
   */
  LatticeFermion apply (const LatticeFermion& psi, enum LinOpSign isign, int cb) const;

  //! Subset is all here
  const OrderedSubset& subset() const {return all;}

private:
  multi3d<ColorMatrix> myu;
  
// Real CoeffWilsr_s;
};

#endif

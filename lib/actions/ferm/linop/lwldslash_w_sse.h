// -*- C++ -*-
// $Id: lwldslash_w_sse.h,v 1.5 2003-11-05 18:35:14 edwards Exp $
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

                                                                                
typedef PColorMatrix < RComplex <PScalar<REAL> >, Nc > PrimitiveSU3Matrix;

class SSEWilsonDslash : public DslashLinearOperator<LatticeFermion>
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
  multi1d<PrimitiveSU3Matrix> packed_gauge;
  
// Real CoeffWilsr_s;
};

#endif

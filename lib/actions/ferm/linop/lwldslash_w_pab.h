// -*- C++ -*-
// $Id: lwldslash_w_pab.h,v 1.3 2004-12-12 21:22:16 edwards Exp $
/*! \file
 *  \brief Wilson Dslash linear operator
 */

#ifndef __lwldslash_pab_h__
#define __lwldslash_pab_h__

#include "linearop.h"
#include <wfm.h>

using namespace QDP;

namespace Chroma 
{ 
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

                                                                                
typedef PColorMatrix < RComplex <REAL>, Nc > PrimitiveSU3Matrix;

class PABWilsonDslash : public DslashLinearOperator<LatticeFermion>
{
public:
  //! Empty constructor. Must use create later
  PABWilsonDslash() {}

  //! Full constructor
  PABWilsonDslash(const multi1d<LatticeColorMatrix>& _u) {create(_u);}

  //! Creation routine
  void create(const multi1d<LatticeColorMatrix>& _u);

  //! No real need for cleanup here
  ~PABWilsonDslash();

  /**
   * Apply a dslash
   *
   * \param chi     result                                      (Write)
   * \param psi     source                                      (Read)
   * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
   * \param cb	    Checkerboard of OUTPUT vector               (Read) 
   *
   * \return The output of applying dslash on psi
   */
  void apply (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign, int cb) const;

  //! Subset is all here
  const OrderedSubset& subset() const {return all;}

private:
  multi1d<PrimitiveSU3Matrix> packed_gauge;
  WilsonArg wil;
  unsigned long wil_cbsize;
// Real CoeffWilsr_s;
};


}; // End Namespace Chroma

using namespace Chroma;

#endif

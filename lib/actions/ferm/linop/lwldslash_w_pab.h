// -*- C++ -*-
// $Id: lwldslash_w_pab.h,v 1.4 2004-12-14 05:21:32 edwards Exp $
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

  //! Subset is all here
  const OrderedSubset& subset() const {return all;}

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

  //! Take deriv of D
  /*!
   * \param chi     left vector                                 (Read)
   * \param psi     right vector                                (Read)
   * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
   *
   * \return Computes   chi^dag * \dot(D} * psi  
   */
  void deriv(multi1d<LatticeColorMatrix>& ds_u, 
	     const LatticeFermion& chi, const LatticeFermion& psi, 
	     enum PlusMinus isign) const
  {
    ds_u.resize(Nd);

    multi1d<LatticeColorMatrix> ds_tmp(Nd);
    deriv(ds_u, chi, psi, isign, 0);
    deriv(ds_tmp, chi, psi, isign, 1);
    ds_u += ds_tmp;
  }


  //! Take deriv of D
  /*!
   * \param chi     left vector on cb                           (Read)
   * \param psi     right vector on 1-cb                        (Read)
   * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
   * \param cb	    Checkerboard of chi vector                  (Read)
   *
   * \return Computes   chi^dag * \dot(D} * psi  
   */
  void deriv(multi1d<LatticeColorMatrix>& ds_u, 
	     const LatticeFermion& chi, const LatticeFermion& psi, 
	     enum PlusMinus isign, int cb) const;


private:
  multi1d<PrimitiveSU3Matrix> packed_gauge;
  multi1d<LatticeColorMatrix> u;   // should not need this

  WilsonArg wil;
  unsigned long wil_cbsize;
// Real CoeffWilsr_s;
};


}; // End Namespace Chroma

using namespace Chroma;

#endif

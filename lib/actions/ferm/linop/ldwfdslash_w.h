// -*- C++ -*-
// $Id: ldwfdslash_w.h,v 1.3 2003-11-20 05:43:41 edwards Exp $
/*! \file
 *  \brief DW Dslash linear operator
 */

#ifndef __ldwfdslash_h__
#define __ldwfdslash_h__

#include "linearop.h"

using namespace QDP;

//! The even-odd preconditioned 5D DW dslash operator
/*!
 * \ingroup linop
 *
 * This routine is specific to Wilson fermions!
 *
 */

class DWDslash : public DslashLinearOperator<LatticeDWFermion>
{
public:
  //! Empty constructor. Must use create later
  DWDslash() {}

  //! Full constructor
  DWDslash(const multi1d<LatticeColorMatrix>& u_, const Real& WilsonMass_) 
    {create(u_, WilsonMass_);}

  //! Creation routine
  void create(const multi1d<LatticeColorMatrix>& u_, const Real& WilsonMass_);

  //! No real need for cleanup here
  ~DWDslash() {}

  /**
   * Apply a dslash
   *
   * \param chi     rsult                                       (Write)
   * \param psi     source                                      (Read)
   * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
   * \param cb	    Checkerboard of OUTPUT vector               (Read) 
   *
   * \return The output of applying dslash on psi
   */
  void apply (LatticeDWFermion& chi, const LatticeDWFermion& psi, 
	      enum PlusMinus isign, int cb) const;

  //! Subset is all here
  const OrderedSubset& subset() const {return all;}

private:
  multi1d<LatticeColorMatrix> u;
  Real  WilsonMass;
// Real CoeffWilsr_s;
};

#endif

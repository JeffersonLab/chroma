// -*- C++ -*-
// $Id: unprec_dwf_linop_w.h,v 1.2 2003-11-08 04:21:47 edwards Exp $
/*! \file
 *  \brief Unpreconditioned domain-wall fermion linear operator
 */

#ifndef __unprec_dwf_linop_w_h__
#define __unprec_dwf_linop_w_h__

#include "linearop.h"
#include "actions/ferm/linop/ldwfdslash_w.h"

using namespace QDP;

//! Unpreconditioned domain-wall Dirac operator
/*!
 * \ingroup linop
 *
 * This routine is specific to Wilson fermions!
 */

class UnprecDWLinOp : public LinearOperator<LatticeDWFermion>
{
public:
  //! Partial constructor
  UnprecDWLinOp() {}

  //! Full constructor
  UnprecDWLinOp(const multi1d<LatticeColorMatrix>& u, const Real& WilsonMass, const Real& m_q)
    {create(u,WilsonMass,m_q);}

  //! Creation routine
  void create(const multi1d<LatticeColorMatrix>& u, const Real& WilsonMass, const Real& m_q);

  //! Destructor is automatic
  ~UnprecDWLinOp() {}

  //! Only defined on the odd subset
  const OrderedSubset& subset() const {return all;}

  //! Apply the operator onto a source vector
  LatticeDWFermion operator() (const LatticeDWFermion& psi, enum LinOpSign isign) const;

private:
  Real WilsonMass;
  Real m_q;
  Real a5;
  multi1d<LatticeColorMatrix> u;
};

#endif

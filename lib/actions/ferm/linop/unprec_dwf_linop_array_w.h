// -*- C++ -*-
// $Id: unprec_dwf_linop_array_w.h,v 1.2 2003-11-13 04:13:06 edwards Exp $
/*! \file
 *  \brief Unpreconditioned domain-wall fermion linear operator
 */

#ifndef __unprec_dwf_linop_array_w_h__
#define __unprec_dwf_linop_array_w_h__

#include "linearop.h"
#include "actions/ferm/linop/dslash_w.h"

using namespace QDP;

//! Unpreconditioned domain-wall Dirac operator
/*!
 * \ingroup linop
 *
 * This routine is specific to Wilson fermions!
 */

class UnprecDWLinOpArray : public LinearOperator< multi1d<LatticeFermion> >
{
public:
  //! Partial constructor
  UnprecDWLinOpArray() {}

  //! Full constructor
  UnprecDWLinOpArray(const multi1d<LatticeColorMatrix>& u_, const Real& WilsonMass_, const Real& m_q, int N5_)
    {create(u_,WilsonMass_,m_q,N5_);}

  //! Creation routine
  void create(const multi1d<LatticeColorMatrix>& u_, const Real& WilsonMass_, const Real& m_q_, int N5_);

  //! Destructor is automatic
  ~UnprecDWLinOpArray() {}

  //! Only defined on the odd subset
  const OrderedSubset& subset() const {return all;}

  //! Apply the operator onto a source vector
  multi1d<LatticeFermion> operator() (const multi1d<LatticeFermion>& psi, enum PlusMinus isign) const;

private:
  Real WilsonMass;
  Real m_q;
  Real a5;
  int  N5;
  WilsonDslash  D;
};

#endif

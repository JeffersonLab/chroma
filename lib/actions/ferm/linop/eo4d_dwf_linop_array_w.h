// -*- C++ -*-
// $Id: eo4d_dwf_linop_array_w.h,v 1.1 2003-11-22 19:37:23 kostas Exp $
/*! \file
 *  \brief 4D Even Odd preconditioned domain-wall fermion linear operator
 */

#ifndef __eo4d_dwf_linop_array_w_h__
#define __eo4d_dwf_linop_array_w_h__

#include "linearop.h"
#include "actions/ferm/linop/dslash_w.h"

using namespace QDP;

//! 4D Even Odd preconditioned domain-wall Dirac operator
/*!
 * \ingroup linop
 *
 * This routine is specific to Wilson fermions!
 */

class EvenOdd4dDWLinOpArray : public LinearOperator< multi1d<LatticeFermion> >
{
public:
  //! Partial constructor
  EvenOdd4dDWLinOpArray() {}

  //! Full constructor
  EvenOdd4dDWLinOpArray(const multi1d<LatticeColorMatrix>& u_, const Real& WilsonMass_, const Real& m_q, int N5_)
    {create(u_,WilsonMass_,m_q,N5_);}

  //! Creation routine
  void create(const multi1d<LatticeColorMatrix>& u_, const Real& WilsonMass_, const Real& m_q_, int N5_);

  //! Destructor is automatic
  ~EvenOdd4dDWLinOpArray() {}

  //! Only defined on the odd sublattice
  const OrderedSubset& subset() const {return odd;}

  //! Apply the operator onto a source vector
  void operator()(multi1d<LatticeFermion>& chi, 
		  const multi1d<LatticeFermion>& psi, 
		  enum PlusMinus isign) const;

private:
  Real WilsonMass;
  Real m_q;
  Real a5;
  int  N5;
  EvenOdd4dDWDslashArray  Qslash;
  EvenOdd4dDWDiagArray    Qdiag ;
};

#endif

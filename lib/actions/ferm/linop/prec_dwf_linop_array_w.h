// -*- C++ -*-
// $Id: prec_dwf_linop_array_w.h,v 1.1 2003-11-22 21:34:01 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned domain-wall fermion linear operator
 */

#ifndef __unprec_dwf_linop_array_w_h__
#define __unprec_dwf_linop_array_w_h__

#include "linearop.h"
#include "actions/ferm/linop/dslash_w.h"

using namespace QDP;

//! Even-odd preconditioned domain-wall Dirac operator
/*!
 * \ingroup linop
 *
 * This routine is specific to Wilson fermions!
 */

class EvenOddPrecDWLinOpArray : public EvenOddPrecLinearOperator< multi1d<LatticeFermion> >
{
public:
  //! Partial constructor
  EvenOddPrecDWLinOpArray() {}

  //! Full constructor
  EvenOddPrecDWLinOpArray(const multi1d<LatticeColorMatrix>& u_, 
			  const Real& WilsonMass_, const Real& m_q, int N5_)
    {create(u_,WilsonMass_,m_q,N5_);}

  //! Creation routine
  void create(const multi1d<LatticeColorMatrix>& u_, 
	      const Real& WilsonMass_, const Real& m_q_, int N5_);

  //! Destructor is automatic
  ~EvenOddPrecDWLinOpArray() {}

  //! Apply the the even-odd block onto a source vector
  void evenEvenInvLinOp(multi1d<LatticeFermion>& chi, 
			const multi1d<LatticeFermion>& psi, 
			enum PlusMinus isign) const;

  //! Apply the the odd-even block onto a source vector
  void evenOddLinOp(multi1d<LatticeFermion>& chi, 
		    const multi1d<LatticeFermion>& psi, 
		    enum PlusMinus isign) const;

  //! Apply the the odd-odd block onto a source vector
  void oddEvenLinOp(multi1d<LatticeFermion>& chi, 
		    const multi1d<LatticeFermion>& psi, 
		    enum PlusMinus isign) const;

  //! Apply the operator onto a source vector
  void oddOddLinOp(multi1d<LatticeFermion>& chi, 
		   const multi1d<LatticeFermion>& psi, 
		   enum PlusMinus isign) const;

  //! Apply the operator onto a source vector
  void operator() (multi1d<LatticeFermion>& chi, 
		   const multi1d<LatticeFermion>& psi, 
		   enum PlusMinus isign) const;

private:
  Real WilsonMass;
  Real m_q;
  Real a5;
  int  N5;
  WilsonDslash  D;
};

#endif

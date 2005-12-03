// -*- C++ -*-
// $Id: unprec_hamberwu_linop_w.h,v 1.1 2005-12-03 04:20:20 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Hamber-Wu fermion linear operator
 */

#ifndef __unprec_hamberwu_linop_w_h__
#define __unprec_hamberwu_linop_w_h__

#include "linearop.h"
#include "actions/ferm/linop/dslash_w.h"


namespace Chroma 
{ 
  //! Unpreconditioned Hamber-Wu operator
  /*!
   * \ingroup linop
   */
  
  class UnprecHamberWuLinOp : public UnprecLinearOperator<LatticeFermion, multi1d<LatticeColorMatrix> >
  {
  public:
    //! Partial constructor
    UnprecHamberWuLinOp() {}

    //! Full constructor
    UnprecHamberWuLinOp(const multi1d<LatticeColorMatrix>& u_, 
			const Real& Mass_, const Real& u0_)
      {create(u_,Mass_,u0_);}

    //! Destructor is automatic
    ~UnprecHamberWuLinOp() {}

    //! Creation routine
    void create(const multi1d<LatticeColorMatrix>& u_, const Real& Mass_, const Real& u0_);

    //! Apply the operator onto a source vector
    void operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const;

    //! Derivative of unpreconditioned HamberWu dM/dU
    void deriv(multi1d<LatticeColorMatrix>& ds_u, 
	       const LatticeFermion& chi, const LatticeFermion& psi, 
	       enum PlusMinus isign) const;

    //! Return flops performed by the operator()
    unsigned long UnprecHamberWuLinOp::nFlops() const;

  private:
    Real Mass;
    Real u0;
    Real fact1;
    Real fact2;
    Real fact3;
    Real fact4;
    multi1d<LatticeColorMatrix> u_dble;
    WilsonDslash D;
  };

}; // End Namespace Chroma


#endif

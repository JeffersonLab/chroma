// -*- C++ -*-
// $Id: unprec_clover_linop_w.h,v 2.2 2005-12-18 23:53:26 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Clover fermion linear operator
 */

#ifndef __unprec_clover_linop_w_h__
#define __unprec_clover_linop_w_h__

#include "linearop.h"
#include "actions/ferm/linop/dslash_w.h"
#include "actions/ferm/linop/clover_term_w.h"


namespace Chroma 
{ 
  //! Unpreconditioned Clover-Dirac operator
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   */
  
  class UnprecCloverLinOp : public UnprecLinearOperator<LatticeFermion, multi1d<LatticeColorMatrix> >
  {
  public:
    //! Partial constructor
    UnprecCloverLinOp() {}

    //! Full constructor
    UnprecCloverLinOp(const multi1d<LatticeColorMatrix>& u_, const CloverFermActParams& param_)
      {create(u_,param_);}
    
    //! Destructor is automatic
    ~UnprecCloverLinOp() {}

    //! Creation routine
    void create(const multi1d<LatticeColorMatrix>& u_, 	
		const CloverFermActParams& param_);

    //! Apply the operator onto a source vector
    void operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const;

    //! Derivative of unpreconditioned Clover dM/dU
    void deriv(multi1d<LatticeColorMatrix>& ds_u, 
	       const LatticeFermion& chi, const LatticeFermion& psi, 
	       enum PlusMinus isign) const;

    //! Return flops performed by the operator()
    unsigned long nFlops() const;

  private:
    CloverFermActParams param;
    WilsonDslash        D;
    CloverTerm          A;
  };

}; // End Namespace Chroma


#endif

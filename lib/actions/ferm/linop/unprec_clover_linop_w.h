// -*- C++ -*-
// $Id: unprec_clover_linop_w.h,v 3.0 2006-04-03 04:58:51 edwards Exp $
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
  
  class UnprecCloverLinOp : public UnprecLinearOperator<LatticeFermion, 
			    multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Partial constructor
    UnprecCloverLinOp() {}

    //! Full constructor
    UnprecCloverLinOp(Handle< FermState<T,P,Q> > fs,
		      const CloverFermActParams& param_)
      {create(fs,param_);}
    
    //! Destructor is automatic
    ~UnprecCloverLinOp() {}

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return D.getFermBC();}

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > fs,
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

} // End Namespace Chroma


#endif

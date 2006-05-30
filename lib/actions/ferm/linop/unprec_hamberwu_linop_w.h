// -*- C++ -*-
// $Id: unprec_hamberwu_linop_w.h,v 3.1 2006-05-30 19:56:58 edwards Exp $
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
  
  class UnprecHamberWuLinOp : public UnprecLinearOperator<LatticeFermion, 
			      multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Partial constructor
    UnprecHamberWuLinOp() {}

    //! Full constructor
    UnprecHamberWuLinOp(Handle< FermState<T,P,Q> > fs,
			const Real& Mass_, const Real& u0_)
      {create(fs,Mass_,u0_);}

    //! Destructor is automatic
    ~UnprecHamberWuLinOp() {}

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return D.getFermBC();}

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > fs, const Real& Mass_, const Real& u0_);

    //! Apply the operator onto a source vector
    void operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const;

    //! Derivative of unpreconditioned HamberWu dM/dU
    void deriv(multi1d<LatticeColorMatrix>& ds_u, 
	       const LatticeFermion& chi, const LatticeFermion& psi, 
	       enum PlusMinus isign) const;

    //! Return flops performed by the operator()
    unsigned long nFlops() const;

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

} // End Namespace Chroma


#endif

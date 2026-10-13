// -*- C++ -*-
/*! \file
 *  \brief Unpreconditioned Clover fermion linear operator
 */

#ifndef __unprec_clover_linop_w_h__
#define __unprec_clover_linop_w_h__

#include "linearop.h"
#include "coarsenable_op.h"
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
  
  class UnprecCloverLinOp : public UnprecLinearOperator<LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >,
			    public MGCoarsenableOperator<LatticeFermion>
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

    //! Apply the operator onto a source std::vector
    void operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const;

    //! Derivative of unpreconditioned Clover dM/dU
    void deriv(multi1d<LatticeColorMatrix>& ds_u, 
	       const LatticeFermion& chi, const LatticeFermion& psi, 
	       enum PlusMinus isign) const;

    MGOperatorForm mgForm() const override
    { return MGOperatorForm::UNPRECONDITIONED; }

    const Subset& paritySubset(int p) const override { return rb[p]; }

    bool hasLocal    (int p, int q) const override { return p == q; }
    bool hasDirection(int p, int q) const override { return p != q; }

    void applyLocal(LatticeFermion& chi, const LatticeFermion& psi,
                    enum PlusMinus isign, int p, int q) const override
    {
      assert(p == q);
      A.apply(chi, psi, isign, p);             // clover incl. mass, cb = p
      getFermBC().modifyF(chi, rb[p]);
    }

    void applyDirection(LatticeFermion& chi, const LatticeFermion& psi,
                        enum PlusMinus isign, int dir, int p, int q) const override
    {
      assert(q == 1 - p);
      D.applyDirection(chi, psi, isign, dir, p); // output cb = p
      Real mhalf = -0.5;
      chi[rb[p]] *= mhalf;
      getFermBC().modifyF(chi, rb[p]);
    }

    //! Return flops performed by the operator()
    unsigned long nFlops() const;

  private:
    CloverFermActParams param;
    WilsonDslash        D;
    CloverTerm          A;
  };




  
} // End Namespace Chroma


#endif

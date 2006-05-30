// -*- C++ -*-
// $Id: unprec_w12_linop_w.h,v 3.1 2006-05-30 19:57:22 edwards Exp $
/*! \file
 *  \brief Unpreconditioned W12 action
 */

#ifndef __unprec_w12_linop_w_h__
#define __unprec_w12_linop_w_h__

#include "linearop.h"
#include "actions/ferm/linop/clover_term_w.h"


namespace Chroma 
{ 
  //! Unpreconditioned W12 operator
  /*!
   * \ingroup linop
   *
   * The W12 action does
   * Chi = (m0 - (2/3*((1/2)*(1/4))*sigma.F + W'  + (1/6)*W^2_mu) * Psi
   */
  class UnprecW12LinOp : public UnprecLinearOperator<LatticeFermion, 
			 multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Partial constructor
    UnprecW12LinOp() {}

    //! Full constructor
    UnprecW12LinOp(Handle< FermState<T,P,Q> > fs, const CloverFermActParams& param_)
      {create(fs,param_);}

    //! Destructor is automatic
    ~UnprecW12LinOp() {}

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return A.getFermBC();}

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > fs,
		const CloverFermActParams& param_);

    //! Apply the operator onto a source vector
    void operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const;

    //! Derivative of unpreconditioned W12 dM/dU
    void deriv(multi1d<LatticeColorMatrix>& ds_u, 
	       const LatticeFermion& chi, const LatticeFermion& psi, 
	       enum PlusMinus isign) const;

    //! Return flops performed by the operator()
    unsigned long nFlops() const;

  protected:
    //! GAMWM
    /*!
     * Description:
     *
     * This routine applies the operator W' to Psi, putting the result in Chi.
     *
     *   chi(x,mu)  :=  + (1 - isign gamma  ) U  (x) psi(x+mu)
     *                                    mu   mu
     *                                         +
     *                  + (1 + isign gamma  ) U  (x-mu) psi(x-mu)
     *                                    mu   mu
     */
    void gamW(multi1d<LatticeFermion>& chi, 
	      const LatticeFermion& psi,
	      int j_decay,
	      enum PlusMinus isign) const;


    //! GAMWMUM 
    /*! This routine applies the operator W' to Psi, putting the result in Chi.
     *
     *   chi(x,mu)  :=  + (1 - isign gamma  ) U  (x) psi  (x+mu)
     *                                    mu   mu       mu
     *                                         +
     *                  + (1 + isign gamma  ) U  (x-mu) psi  (x-mu)
     *                                    mu   mu          mu
     */
    void gamWmu(multi1d<LatticeFermion>& chi, 
		const multi1d<LatticeFermion>& psi,
		int j_decay,
		enum PlusMinus isign) const;

  private:
    CloverFermActParams param;
    Real aniso_fact;
    Real fact1;
    Real fact2;
    Real fact3;
    Real fact4;
    int  j_decay;
    multi1d<LatticeColorMatrix> u;  // fold in anisotropy
    CloverTerm          A;
  };

} // End Namespace Chroma


#endif

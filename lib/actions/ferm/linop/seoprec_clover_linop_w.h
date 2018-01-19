// -*- C++ -*-
/*! \file
 *  \brief Symmetric even-odd preconditioned Clover fermion linear operator
 */

#ifndef __seoprec_clover_linop_w_h__
#define __seoprec_clover_linop_w_h__

#include "state.h"
#include "fermbc.h"
#include "seoprec_logdet_linop.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "actions/ferm/linop/dslash_w.h"
#include "actions/ferm/linop/clover_term_w.h"


namespace Chroma 
{ 
  //! Symmetric even-odd preconditioned Clover-Dirac operator
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   *
   * The kernel for Clover fermions is
   *
   *      M  =  A + (d+M) - (1/2) D'
   */
  class SymEvenOddPrecCloverLinOp : public SymEvenOddPrecLogDetLinearOperator<LatticeFermion, 
									      multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Partial constructor
    SymEvenOddPrecCloverLinOp() {}

    //! Full constructor
    SymEvenOddPrecCloverLinOp(Handle< FermState<T,P,Q> > fs,
			      const CloverFermActParams& param_)
    {
      create(fs,param_);
    }

    //! Destructor is automatic
    ~SymEvenOddPrecCloverLinOp() {
      QDPIO::cout << "CLOV_LINOP: Time spent in clov deriv (total) = " << clov_deriv_time << std::endl;
      QDPIO::cout << "CLOV_LINOP: Time spent in clov apply/invapply (total) = " << clov_apply_time << std::endl;
    }

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return D.getFermBC();}

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > fs,
		const CloverFermActParams& param_);

    //! Apply the the even-even block onto a source std::vector
    void evenEvenLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
		       enum PlusMinus isign) const;

    //! Apply the inverse of the even-even block onto a source std::vector
    void evenEvenInvLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
			  enum PlusMinus isign) const;
  
    //! Apply the the even-odd block onto a source std::vector
    void evenOddLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
		      enum PlusMinus isign) const;

    //! Apply the the odd-even block onto a source std::vector
    void oddEvenLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
		      enum PlusMinus isign) const;

    //! Apply the the odd-odd block onto a source std::vector
    void oddOddLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
		     enum PlusMinus isign) const;

    //! Apply the inverse of the odd-odd block onto a source std::vector
    void oddOddInvLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
			enum PlusMinus isign) const;
  
    // Override inherited one with a few more funkies
    void operator()(LatticeFermion& chi, const LatticeFermion& psi, 
		    enum PlusMinus isign) const;

    //! Apply the even-even block onto a source std::vector
    void derivEvenEvenLinOp(multi1d<LatticeColorMatrix>& ds_u, 
			    const LatticeFermion& chi, const LatticeFermion& psi, 
			    enum PlusMinus isign) const;

    //! Apply the odd-odd block onto a source std::vector
    void derivOddOddLinOp(multi1d<LatticeColorMatrix>& ds_u, 
			  const LatticeFermion& chi, const LatticeFermion& psi, 
			  enum PlusMinus isign) const;

    //! Apply the even-even block onto a source std::vector
    void derivLogDetEvenEvenLinOp(multi1d<LatticeColorMatrix>& ds_u, 
				  enum PlusMinus isign) const;

    //! Apply the odd-odd block onto a source std::vector
    void derivLogDetOddOddLinOp(multi1d<LatticeColorMatrix>& ds_u, 
				enum PlusMinus isign) const;

    //! Apply the the even-odd block onto a source std::vector
    void derivEvenOddLinOp(multi1d<LatticeColorMatrix>& ds_u, 
			   const LatticeFermion& chi, const LatticeFermion& psi, 
			   enum PlusMinus isign) const;
 
    //! Apply the the odd-even block onto a source std::vector
    void derivOddEvenLinOp(multi1d<LatticeColorMatrix>& ds_u, 
			   const LatticeFermion& chi, const LatticeFermion& psi, 
			   enum PlusMinus isign) const;

    //! Return flops performed by the operator()
    unsigned long nFlops() const;

    //! Get the log det of the even even part
    Double logDetEvenEvenLinOp(void) const; 

    //! Get the log det of the even even part
    Double logDetOddOddLinOp(void) const; 

  private:
    CloverFermActParams param;
    WilsonDslash D;
    CloverTerm   clov;
    CloverTerm   invclov;  // needed for even-even and odd-odd
    mutable double clov_apply_time;
    mutable double clov_deriv_time;
    mutable StopWatch swatch;
  };



} // End Namespace Chroma


#endif

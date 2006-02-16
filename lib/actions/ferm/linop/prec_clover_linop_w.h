// -*- C++ -*-
// $Id: prec_clover_linop_w.h,v 2.5 2006-02-16 02:24:46 bjoo Exp $
/*! \file
 *  \brief Even-odd preconditioned Clover fermion linear operator
 */

#ifndef __prec_clover_linop_w_h__
#define __prec_clover_linop_w_h__

#include "prec_logdet_linop.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "actions/ferm/linop/dslash_w.h"
#include "actions/ferm/linop/clover_term_w.h"


namespace Chroma 
{ 
  //! Even-odd preconditioned Clover-Dirac operator
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   *
   * The kernel for Clover fermions is
   *
   *      M  =  A + (d+M) - (1/2) D'
   */
  class EvenOddPrecCloverLinOp : public EvenOddPrecLogDetLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> >
  {
  public:
    //! Partial constructor
    EvenOddPrecCloverLinOp() {}

    //! Full constructor
    EvenOddPrecCloverLinOp(const multi1d<LatticeColorMatrix>& u_, const CloverFermActParams& param_)
      {create(u_,param_);}

    //! Destructor is automatic
    ~EvenOddPrecCloverLinOp() {}

    //! Creation routine
    void create(const multi1d<LatticeColorMatrix>& u_, 	
		const CloverFermActParams& param_);

    //! Apply the the even-even block onto a source vector
    void evenEvenLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
		       enum PlusMinus isign) const;

    //! Apply the inverse of the even-even block onto a source vector
    void evenEvenInvLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
			  enum PlusMinus isign) const;
  
    //! Apply the the even-odd block onto a source vector
    void evenOddLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
		      enum PlusMinus isign) const;

    //! Apply the the odd-even block onto a source vector
    void oddEvenLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
		      enum PlusMinus isign) const;

    //! Apply the the odd-odd block onto a source vector
    void oddOddLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
		     enum PlusMinus isign) const;

    // Override inherited one with a few more funkies
    void operator()(LatticeFermion& chi, const LatticeFermion& psi, 
		    enum PlusMinus isign) const;

    //! Apply the even-even block onto a source vector
    void derivEvenEvenLinOp(multi1d<LatticeColorMatrix>& ds_u, 
			    const LatticeFermion& chi, const LatticeFermion& psi, 
			    enum PlusMinus isign) const;

    void derivEvenEvenLogDet(multi1d<LatticeColorMatrix>& ds_u, 
			     enum PlusMinus isign) const;

    //! Apply the the even-odd block onto a source vector
    void derivEvenOddLinOp(multi1d<LatticeColorMatrix>& ds_u, 
			   const LatticeFermion& chi, const LatticeFermion& psi, 
			   enum PlusMinus isign) const;
 
    //! Apply the the odd-even block onto a source vector
    void derivOddEvenLinOp(multi1d<LatticeColorMatrix>& ds_u, 
			   const LatticeFermion& chi, const LatticeFermion& psi, 
			   enum PlusMinus isign) const;

    //! Apply the the odd-odd block onto a source vector
    void derivOddOddLinOp(multi1d<LatticeColorMatrix>& ds_u, 
			  const LatticeFermion& chi, const LatticeFermion& psi, 
			  enum PlusMinus isign) const;

    //! Return flops performed by the operator()
    unsigned long nFlops() const;

    //! Get the log det of the even even part
    Double LogDetEvenEven(void) const ; 

  private:
    CloverFermActParams param;
    WilsonDslash D;
    CloverTerm   clov;
    CloverTerm   invclov;  // uggh, only needed for evenEvenLinOp
  };

}; // End Namespace Chroma


#endif

// -*- C++ -*-
// $Id: unprec_clover_linop_w.h,v 2.1 2005-12-03 21:19:38 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Clover fermion linear operator
 */

#ifndef __unprec_clover_linop_w_h__
#define __unprec_clover_linop_w_h__

#include "linearop.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "actions/ferm/linop/dslash_w.h"


namespace Chroma 
{ 
  //! Unpreconditioned Clover-Dirac operator
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   *
   * The kernel for Clover fermions is
   *
   *      M  =  (d+M) - (1/2) D'
   */
  class UnprecCloverLinOp : public UnprecLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> >
  {
  public:
    //! Full constructor
    UnprecCloverLinOp(const multi1d<LatticeColorMatrix>& u_, const CloverFermActParams& param_)
      {create(u_,param_);}

    //! Destructor is automatic
    ~UnprecCloverLinOp() {}

    //! Creation routine
    void create(const multi1d<LatticeColorMatrix>& u_, const CloverFermActParams& param_);

    //! Apply the the even-even block onto a source vector
    inline
    void evenEvenLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
		       enum PlusMinus isign) const
    {
      chi[rb[0]] = fact*psi;
    }

    //! Apply the inverse of the even-even block onto a source vector
    inline 
    void evenEvenInvLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
			  enum PlusMinus isign) const
    {
      chi[rb[0]] = invfact*psi;
    }
  
    //! Apply the the even-odd block onto a source vector
    void evenOddLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
		      enum PlusMinus isign) const;

    //! Apply the the odd-even block onto a source vector
    void oddEvenLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
		      enum PlusMinus isign) const;

    //! Apply the the odd-odd block onto a source vector
    inline 
    void oddOddLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
		     enum PlusMinus isign) const
    {
      chi[rb[1]] = fact*psi;
    }

    // Override inherited one with a few more funkies
    void operator()(LatticeFermion& chi, const LatticeFermion& psi, 
		    enum PlusMinus isign) const;


    //! Apply the even-even block onto a source vector
    void derivEvenEvenLinOp(multi1d<LatticeColorMatrix>& ds_u, 
			    const LatticeFermion& chi, const LatticeFermion& psi, 
			    enum PlusMinus isign) const
    {
      QDPIO::cerr << "Clover: not implemented" << endl;
      QDP_abort(1);
    }
  
    //! Apply the the even-odd block onto a source vector
    void derivEvenOddLinOp(multi1d<LatticeColorMatrix>& ds_u, 
			   const LatticeFermion& chi, const LatticeFermion& psi, 
			   enum PlusMinus isign) const
    {
      QDPIO::cerr << "Clover: not implemented" << endl;
      QDP_abort(1);
    }
 
    //! Apply the the odd-even block onto a source vector
    void derivOddEvenLinOp(multi1d<LatticeColorMatrix>& ds_u, 
			   const LatticeFermion& chi, const LatticeFermion& psi, 
			   enum PlusMinus isign) const
    {
      QDPIO::cerr << "Clover: not implemented" << endl;
      QDP_abort(1);
    }

    //! Apply the the odd-odd block onto a source vector
    void derivOddOddLinOp(multi1d<LatticeColorMatrix>& ds_u, 
			  const LatticeFermion& chi, const LatticeFermion& psi, 
			  enum PlusMinus isign) const
    {
      QDPIO::cerr << "Clover: not implemented" << endl;
      QDP_abort(1);
    }

  private:
    Real fact;
    Real invfact;
    CloverFermActParams param;
    WilsonDslash D;
  };

} // End Namespace Chroma


#endif

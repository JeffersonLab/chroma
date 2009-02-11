// -*- C++ -*-
// $Id: eoprec_clover_orbifold_linop_w.h,v 1.3 2009-02-11 06:17:16 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Clover fermion linear operator with orbifold
 *
 * 3D-Orbifold construction follows arXiv:0811.2127
 */

#ifndef __prec_clover_orbifold_linop_w_h__
#define __prec_clover_orbifold_linop_w_h__

#include "state.h"
#include "fermbc.h"
#include "eoprec_logdet_linop.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "actions/ferm/linop/dslash_w.h"
#include "actions/ferm/linop/clover_term_w.h"


namespace Chroma 
{ 
  //! Even-odd preconditioned Clover-Dirac operator with orbifold term
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   *
   * The kernel for Clover fermions is
   *
   *      M  =  A + (d+M) - (1/2) D' + 3d-orbifold
   */
  class EvenOddPrecCloverOrbifoldLinOp : public EvenOddPrecLogDetLinearOperator<LatticeFermion, 
				 multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Partial constructor
    EvenOddPrecCloverOrbifoldLinOp() {}

    //! Full constructor
    EvenOddPrecCloverOrbifoldLinOp(Handle< FermState<T,P,Q> > fs,
			   const CloverFermActParams& param_)
      {create(fs,param_);}

    //! Destructor is automatic
    ~EvenOddPrecCloverOrbifoldLinOp() {}

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return D.getFermBC();}

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > fs,
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

    //! Orbifold term
    void orbifold(LatticeFermion& chi, 
		  const LatticeFermion& psi, 
		  enum PlusMinus isign,
		  int z, int cb) const;

    //! Override inherited one with a few more funkies
    void operator()(LatticeFermion& chi, const LatticeFermion& psi, 
		    enum PlusMinus isign) const;

    //! Return flops performed by the operator()
    unsigned long nFlops() const;
 
    //! Get the log det of the even even part
    Double logDetEvenEvenLinOp(void) const; 

  private:
    CloverFermActParams param;
    WilsonDslash D;
    CloverTerm   clov;
    CloverTerm   invclov;  // uggh, only needed for evenEvenLinOp
  };

} // End Namespace Chroma


#endif

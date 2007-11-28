// -*- C++ -*-
// $Id: eoprec_slrc_linop_w.h,v 3.5 2007-11-28 22:09:15 bjoo Exp $
/*! \file
 *  \brief Even-odd preconditioned Clover linear operator (fat-relevant, thin-irrelevant terms)
 *
 * Here, the relevant terms are smeared and the irrelevant terms are not smeared.
 * Code provided by Thomas Kaltenbrunner.
 *
 */

#ifndef __eoprec_slrc_linop_w_h__
#define __eoprec_slrc_linop_w_h__

#include "state.h"
#include "fermbc.h"
#include "actions/ferm/fermstates/stout_fermstate_w.h"

#include "eoprec_logdet_linop.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "actions/ferm/linop/dslash_w.h"
#include "actions/ferm/linop/clover_term_w.h"

namespace Chroma 
{ 
  //! Even-odd preconditioned SLRC-Dirac operator
  /*!
   * \ingroup linop
   *
   * Here, the relevant terms are smeared and the irrelevant terms are not smeared.
   * The kernel for SLRC fermions is
   *
   *      M  =  A + (d+M) - (1/2) D'
   */
  class EvenOddPrecSLRCLinOp : public EvenOddPrecLogDetLinearOperator<LatticeFermion, 
				 multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Partial constructor
    EvenOddPrecSLRCLinOp() {}

    //! Full constructor
    EvenOddPrecSLRCLinOp(Handle< FermState<T,P,Q> > fs,
			 const CloverFermActParams& param_) : slrc_fs(fs)
      {create(fs, param_);}

    //! Destructor is automatic
    ~EvenOddPrecSLRCLinOp() {}

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return *(slrc_fs->getFermBC());}

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

    // Override inherited one with a few more funkies
    void operator()(LatticeFermion& chi, const LatticeFermion& psi, 
		    enum PlusMinus isign) const;

    //! Apply the even-even block onto a source vector
    void derivEvenEvenLinOp(multi1d<LatticeColorMatrix>& ds_u, 
			    const LatticeFermion& chi, const LatticeFermion& psi, 
			    enum PlusMinus isign) const;

    void derivLogDetEvenEvenLinOp(multi1d<LatticeColorMatrix>& ds_u, 
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
    Double logDetEvenEvenLinOp(void) const ; 

  private:
    Handle< FermState<T,P,Q> > slrc_fs;
    Handle< FermState<T,P,Q> > thin_fs;
    CloverFermActParams param;
    WilsonDslash D;
    CloverTerm   clov;
    CloverTerm   invclov;  // uggh, only needed for evenEvenLinOp
  };

} // End Namespace Chroma


#endif

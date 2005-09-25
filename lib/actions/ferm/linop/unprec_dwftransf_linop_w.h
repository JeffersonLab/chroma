// -*- C++ -*-
// $Id: unprec_dwftransf_linop_w.h,v 2.0 2005-09-25 21:04:30 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion linear operator
 */

#ifndef __unprec_dwftransf_linop_w_h__
#define __unprec_dwftransf_linop_w_h__

#include "linearop.h"
#include "handle.h"
#include "invtype.h"
#include "actions/ferm/linop/unprec_wilson_linop_w.h"


namespace Chroma 
{ 
  //! Apply   gamma_5 * H_T
  /*!
   * where
   *    H_T  =  (b5 + c5) H_w * ( 2 + (b5 + c5) D_w )^(-1)
   *         =  (b5 + c5) ( 2 + (b5 + c5) D_w^dag )^(-1) * H_w
   */
  class UnprecDWFTransfLinOp : public UnprecLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> >
  {
  public:
    //! Partial constructor
    UnprecDWFTransfLinOp() {}

    //! Full constructor
    UnprecDWFTransfLinOp(const multi1d<LatticeColorMatrix>& u_, 
			 const Real& Mass_,
			 const Real& b5_,
			 const Real& c5_,
			 const InvertParam_t& invParam_)
    {create(u_,Mass_, b5_, c5_, invParam_);}

    //! Destructor is automatic
    ~UnprecDWFTransfLinOp() {}

    //! Only defined on the entire lattice
    const OrderedSubset& subset() const {return all;}

    //! Creation routine
    void create(const multi1d<LatticeColorMatrix>& u_, 
		const Real& Mass_,
		const Real& b5_,
		const Real& c5_,
		const InvertParam_t& invParam_);

    //! Apply the operator onto a source vector
    void operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const;


  private:
    Real Mass;
    Real b5;
    Real c5;
    InvertParam_t invParam;
    Handle< LinearOperator<LatticeFermion> > D_w;
    Handle< LinearOperator<LatticeFermion> > D_denum;  
  };


  //! Apply   H_T * H_T
  /*!
   * where
   *    H_T^2 = (b5 + c5)^2 * H_w * [( 2 + (b5 + c5) D_w^dag )*( 2 + (b5 + c5) D_w )]^(-1) * H_w
   */
  class UnprecDWFTransfMdagMLinOp : public UnprecLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> >
  {
  public:
    //! Partial constructor
    UnprecDWFTransfMdagMLinOp() {}

    //! Full constructor
    UnprecDWFTransfMdagMLinOp(const multi1d<LatticeColorMatrix>& u_, 
			      const Real& Mass_,
			      const Real& b5_,
			      const Real& c5_,
			      const InvertParam_t& invParam_)
    {create(u_,Mass_, b5_, c5_, invParam_);}

    //! Destructor is automatic
    ~UnprecDWFTransfMdagMLinOp() {}

    //! Only defined on the entire lattice
    const OrderedSubset& subset() const {return all;}

    //! Creation routine
    void create(const multi1d<LatticeColorMatrix>& u_, 
		const Real& Mass_,
		const Real& b5_,
		const Real& c5_,
		const InvertParam_t& invParam_);

    //! Apply the operator onto a source vector
    void operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const;


  private:
    Real Mass;
    Real b5;
    Real c5;
    InvertParam_t invParam;
    Handle< LinearOperator<LatticeFermion> > D_w;
    Handle< LinearOperator<LatticeFermion> > D_denum;  
  };


  //! Operator to apply the denominator
  /*!
   * ( 2 + (b5 - c5) D )
   */
  class UnprecDWFTransfDenLinOp : public LinearOperator<LatticeFermion>
  {
  public:
    //! Partial constructor
    UnprecDWFTransfDenLinOp() {}

    //! Full constructor
    UnprecDWFTransfDenLinOp(const multi1d<LatticeColorMatrix>& u_, 
			    const Real& b5_minus_c5_,
			    const Handle<LinearOperator<LatticeFermion> > D_w_ ) 
      : b5_minus_c5(b5_minus_c5_), D_w(D_w_) {}


    //! Destructor is automatic
    ~UnprecDWFTransfDenLinOp() {}

    //! Only defined on the entire lattice
    const OrderedSubset& subset() const {return all;}

    //! Apply the operator onto a source vector
    void operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const 
    {

      // Apply chi ( 2 + (b5-c5) D_w ) psi
      (*D_w)(chi, psi, isign);
      chi *= b5_minus_c5;
      chi += Real(2)*psi;
    }


  private:
    Real b5_minus_c5;
    Handle< LinearOperator<LatticeFermion> > D_w;
  };

}; // End Namespace Chroma



#endif

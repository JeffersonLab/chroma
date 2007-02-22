// -*- C++ -*-
// $Id: unprec_dwftransf_linop_w.h,v 3.2 2007-02-22 21:11:47 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion linear operator
 */

#ifndef __unprec_dwftransf_linop_w_h__
#define __unprec_dwftransf_linop_w_h__

#include "linearop.h"
#include "handle.h"
#include "actions/ferm/linop/unprec_wilson_linop_w.h"
#include "actions/ferm/invert/syssolver_cg_params.h"


namespace Chroma 
{ 
  //! Apply   gamma_5 * H_T
  /*!
   * where
   *    H_T  =  (b5 + c5) H_w * ( 2 + (b5 + c5) D_w )^(-1)
   *         =  (b5 + c5) ( 2 + (b5 + c5) D_w^dag )^(-1) * H_w
   */
  class UnprecDWFTransfLinOp : public UnprecLinearOperator<LatticeFermion, 
			       multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Partial constructor
    UnprecDWFTransfLinOp() {}

    //! Full constructor
    UnprecDWFTransfLinOp(Handle< FermState<T,P,Q> > fs,
			 const Real& Mass_,
			 const Real& b5_,
			 const Real& c5_,
			 const SysSolverCGParams& invParam_)
    {create(fs,Mass_, b5_, c5_, invParam_);}

    //! Destructor is automatic
    ~UnprecDWFTransfLinOp() {}

    //! Only defined on the entire lattice
    const Subset& subset() const {return all;}

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > fs,
		const Real& Mass_,
		const Real& b5_,
		const Real& c5_,
		const SysSolverCGParams& invParam_);

    //! Apply the operator onto a source vector
    void operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const;

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return *fbc;}

  private:
    Real Mass;
    Real b5;
    Real c5;
    SysSolverCGParams invParam;
    Handle< LinearOperator<LatticeFermion> > D_w;
    Handle< LinearOperator<LatticeFermion> > D_denum;  
    Handle< FermBC<T,P,Q> > fbc;
  };


  //! Apply   H_T * H_T
  /*!
   * where
   *    H_T^2 = (b5 + c5)^2 * H_w * [( 2 + (b5 + c5) D_w^dag )*( 2 + (b5 + c5) D_w )]^(-1) * H_w
   */
  class UnprecDWFTransfMdagMLinOp : public UnprecLinearOperator<LatticeFermion, 
				    multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Partial constructor
    UnprecDWFTransfMdagMLinOp() {}

    //! Full constructor
    UnprecDWFTransfMdagMLinOp(Handle< FermState<T,P,Q> > fs,
			      const Real& Mass_,
			      const Real& b5_,
			      const Real& c5_,
			      const SysSolverCGParams& invParam_)
    {create(fs, Mass_, b5_, c5_, invParam_);}

    //! Destructor is automatic
    ~UnprecDWFTransfMdagMLinOp() {}

    //! Only defined on the entire lattice
    const Subset& subset() const {return all;}

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return *fbc;}

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > fs,
		const Real& Mass_,
		const Real& b5_,
		const Real& c5_,
		const SysSolverCGParams& invParam_);

    //! Apply the operator onto a source vector
    void operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const;


  private:
    Real Mass;
    Real b5;
    Real c5;
    SysSolverCGParams invParam;
    Handle< LinearOperator<LatticeFermion> > D_w;
    Handle< LinearOperator<LatticeFermion> > D_denum;  
    Handle< FermBC<T,P,Q> > fbc;
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
    UnprecDWFTransfDenLinOp(const Real& b5_minus_c5_,
			    const Handle<LinearOperator<LatticeFermion> > D_w_ ) 
      : b5_minus_c5(b5_minus_c5_), D_w(D_w_) {}


    //! Destructor is automatic
    ~UnprecDWFTransfDenLinOp() {}

    //! Only defined on the entire lattice
    const Subset& subset() const {return all;}

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

} // End Namespace Chroma



#endif

// -*- C++ -*-
// $Id: unprec_ht_contfrac5d_linop_array_w.h,v 3.1 2007-02-22 21:11:47 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned H_T kernel continued fraction (5D) operator
 */

#ifndef __unprec_ht_contfrac5d_linop_array_w_h__
#define __unprec_ht_contfrac5d_linop_array_w_h__

#include "linearop.h"
#include "state.h"


namespace Chroma 
{ 
  //! Unpreconditioned H_T kernel continued fraction (5D) operator
  /*!
   * \ingroup linop
   *
   * This operator applies the non-hermitian overlap operator
   *   Chi  =   ((1+m_q)/(1-m_q)*gamma_5 + B) . Psi
   *  where  B  is the continued fraction of the pole approx. to eps(H_T(m))
   *
   * Here, 
   *    H_T  =  (b5 + c5) H_w * ( 2 + (b5 - c5) D_w )^(-1)
   *         =  (b5 + c5) ( 2 + (b5 - c5) D_w^dag )^(-1) * H_w
   *
   * The denom is multiplied through the 5D system of equations to the RHS
   */

  class UnprecHTContFrac5DLinOpArray : public UnprecLinearOperatorArray<LatticeFermion, 
				       multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Full constructor
    UnprecHTContFrac5DLinOpArray(Handle< FermState<T,P,Q> > fs,
				 const Real& OverMass_, 
				 const Real& _m_q,
				 int _N5,
				 const Real& b5_,
				 const Real& c5_,
				 const multi1d<Real>& _alpha,
				 const multi1d<Real>& _beta,
				 const bool _isLastZeroP ) :
      OverMass(OverMass_), m_q(_m_q), N5(_N5), 
      alpha(_alpha), beta(_beta), 
      isLastZeroP(_isLastZeroP) 
    {
      QDPIO::cout << "LinOp isLastZeroP = " << isLastZeroP << endl;
      init(fs,b5_,c5_);
    }

    //! Length of DW flavor index/space
    int size() const {return N5;}

    //! Destructor is automatic
    ~UnprecHTContFrac5DLinOpArray() {}

    //! Only defined on the entire lattice
    const Subset& subset() const {return all;}

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return *fbc;}

    //! Apply the operator onto a source vector
    void operator() (multi1d<LatticeFermion>& chi, 
		     const multi1d<LatticeFermion>& psi, 
		     enum PlusMinus isign) const;

  protected:
    void init(Handle< FermState<T,P,Q> > fs,
	      const Real& b5, const Real& c5);

    //! Partial constructor
//    UnprecHTContFrac5DLinOpArray() {}
    //! Hide =
    void operator=(const UnprecHTContFrac5DLinOpArray&) {}

  private:
    Handle< LinearOperator<T> > D_w;
    Handle< LinearOperator<T> > D_denum;  
    Handle< FermBC<T,P,Q> > fbc;
    const Real OverMass;
    const Real m_q;
    const int  N5;    // Size of the 5th dimension
    Real scale_fac;
    Real a5;
    const multi1d<Real> alpha;
    const multi1d<Real> beta;
    const bool isLastZeroP;
  };

} // End Namespace Chroma


#endif

// -*- C++ -*-
// $Id: unprec_ht_contfrac5d_linop_array_w.h,v 1.2 2005-01-05 21:44:07 edwards Exp $
/*! \file
 *  \brief Unpreconditioned H_T kernel continued fraction (5D) operator
 */

#ifndef __unprec_ht_contfrac5d_linop_array_w_h__
#define __unprec_ht_contfrac5d_linop_array_w_h__

#include "linearop.h"
#include "state.h"

using namespace QDP;

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

  class UnprecHTContFrac5DLinOpArray : public UnprecLinearOperator< multi1d<LatticeFermion>, multi1d<LatticeColorMatrix> >
  {
  public:

    //! Full constructor
    UnprecHTContFrac5DLinOpArray(const multi1d<LatticeColorMatrix>& u_, 
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
      init(u_,b5_,c5_);
    }

    //! Length of DW flavor index/space
    int size() const {return N5;}

    //! Destructor is automatic
    ~UnprecHTContFrac5DLinOpArray() {}

    //! Only defined on the entire lattice
    const OrderedSubset& subset() const {return all;}

    //! Apply the operator onto a source vector
    void operator() (multi1d<LatticeFermion>& chi, 
		     const multi1d<LatticeFermion>& psi, 
		     enum PlusMinus isign) const;

  protected:
    void init(const multi1d<LatticeColorMatrix>& u,
	      const Real& b5, const Real& c5);

  private:
    Handle< LinearOperator<LatticeFermion> > D_w;
    Handle< LinearOperator<LatticeFermion> > D_denum;  
    const Real OverMass;
    const Real m_q;
    const int  N5;    // Size of the 5th dimension
    Real scale_fac;
    Real a5;
    const multi1d<Real> alpha;
    const multi1d<Real> beta;
    const bool isLastZeroP;
  };

}; // End Namespace Chroma

using namespace Chroma;

#endif

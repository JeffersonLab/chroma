// -*- C++ -*-
// $Id: eoprec_ovlap_contfrac5d_linop_array_opt_w.h,v 3.1 2006-10-19 16:01:30 edwards Exp $
/*! \file
 *  \brief Optimized Even-odd prec. 5D continued fraction linop
 */

#ifndef __prec_ovlap_contfrac5d_linop_array_opt_w_h__
#define __prec_ovlap_contfrac5d_linop_array_opt_w_h__

#include "actions/ferm/linop/eoprec_ovlap_contfrac5d_linop_base_array_w.h"

namespace Chroma 
{ 
  //! Optimized Even-odd prec. 5D continued fraction linop
  /*!
   * \ingroup linop
   *
   * This operator applies the extended version of the hermitian overlap operator
   *   Chi  =   ((1+m_q)/(1-m_q)*gamma_5 + B) . Psi
   *  where  B  is the continued fraction of the pole approx. to eps(H(m))
   *
   * The details of this operator are in some lattice proceedings 
   * by Joo,Kennedy,Wenger
   */

  class OptEvenOddPrecOvlapContFrac5DLinOpArray : public EvenOddPrecOvlapContFrac5DLinOpBaseArray
  {
  public:

    //! Full constructor
    /*! Pretty darn the same as for the unprec case
      except that the auxiliary linop M is no longer supplied, 
      but is created here 
    */
    OptEvenOddPrecOvlapContFrac5DLinOpArray(Handle< FermState<T,P,Q> > state,
					    const Real& _m_q,
					    const Real& _OverMass,
					    int _N5,
					    const Real& _scale_fac,
					    const multi1d<Real>& _alpha,
					    const multi1d<Real>& _beta,
					    const bool _isLastZeroP ) :
      EvenOddPrecOvlapContFrac5DLinOpBaseArray(state,
					       _m_q,
					       _OverMass,
					       _N5,
					       _scale_fac,
					       _alpha,
					       _beta,
					       _isLastZeroP) 
      {}
  
    //! Destructor is automatic
    ~OptEvenOddPrecOvlapContFrac5DLinOpArray() {}

  protected:

#if 0
    //! Apply the even-even (odd-odd) coupling piece of the domain-wall fermion operator
    /*!
     * \param chi     result     	                   (Modify)
     * \param psi     source     	                   (Read)
     * \param isign   Flag ( PLUS | MINUS )          (Read)
     * \param cb      checkerboard ( 0 | 1 )         (Read)
     *
     * Override with more optimized version for scalar-like machines
     */
    virtual
    void applyDiag(multi1d<LatticeFermion>& chi, 
		   const multi1d<LatticeFermion>& psi, 
		   enum PlusMinus isign,
		   const int cb) const;
#endif

    //! Apply the inverse even-even (odd-odd) coupling piece of the domain-wall fermion operator
    /*!
     * \param chi     result     	                   (Modify)
     * \param psi     source     	                   (Read)
     * \param isign   Flag ( PLUS | MINUS )   	   (Read)
     * \param cb      checkerboard ( 0 | 1 )         (Read)
     *
     * Override with more optimized version for scalar-like machines
     */
    void applyDiagInv(multi1d<LatticeFermion>& chi, 
		      const multi1d<LatticeFermion>& psi, 
		      enum PlusMinus isign,
		      const int cb) const;
  };

} // End Namespace Chroma


#endif

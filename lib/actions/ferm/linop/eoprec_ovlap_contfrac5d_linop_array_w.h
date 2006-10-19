// -*- C++ -*-
// $Id: eoprec_ovlap_contfrac5d_linop_array_w.h,v 3.1 2006-10-19 16:01:30 edwards Exp $
/*! \file
 *  \brief Even-odd prec. 5D continued fraction linop
 */

#ifndef __prec_ovlap_contfrac5d_linop_array_w_h__
#define __prec_ovlap_contfrac5d_linop_array_w_h__

#include "actions/ferm/linop/eoprec_ovlap_contfrac5d_linop_base_array_w.h"

namespace Chroma 
{ 
  //! Even-odd prec. 5D continued fraction linop
  /*!
   * \ingroup linop
   *
   * This operator applies the extended version of the hermitian overlap operator
   *   Chi  =   ((1+m_q)/(1-m_q)*gamma_5 + B) . Psi
   *  where  B  is the continued fraction of the pole approx. to eps(H(m))
   *
   * The details of this operator are in some lattice proceedings 
   * by Joo,Kennedy,Wenger
   *
   * Just uses default base class
   */

  class QDPEvenOddPrecOvlapContFrac5DLinOpArray : public EvenOddPrecOvlapContFrac5DLinOpBaseArray
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Full constructor
    /*! Pretty darn the same as for the unprec case
      except that the auxiliary linop M is no longer supplied, 
      but is created here 
    */
    QDPEvenOddPrecOvlapContFrac5DLinOpArray(Handle< FermState<T,P,Q> > state,
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
    ~QDPEvenOddPrecOvlapContFrac5DLinOpArray() {}
  };

} // End Namespace Chroma


#endif

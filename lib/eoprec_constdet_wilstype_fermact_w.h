// -*- C++ -*-
// $Id: eoprec_constdet_wilstype_fermact_w.h,v 3.1 2006-10-19 16:01:26 edwards Exp $

/*! @file
 * @brief Even-odd const determinant Wilson-like fermact
 */

#ifndef __eoprec_constdet_wilstype_fermact_w_h__
#define __eoprec_constdet_wilstype_fermact_w_h__

#include "eoprec_wilstype_fermact_w.h"
#include "eoprec_constdet_linop.h"

namespace Chroma
{
  //-------------------------------------------------------------------------------------------
  //! Even-odd preconditioned Wilson-like fermion actions specialised to Wilson Like (gauge independent diagonal term) actions.
  /*! @ingroup actions
   *
   * Even-odd preconditioned like Wilson-like fermion actions
   */
  template<typename T, typename P, typename Q>
  class EvenOddPrecConstDetWilsonTypeFermAct : public EvenOddPrecWilsonTypeFermAct<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddPrecConstDetWilsonTypeFermAct() {}

    //! Override to produce an even-odd prec. linear operator for this action
    /*! Covariant return rule - override base class function */
    virtual EvenOddPrecConstDetLinearOperator<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const = 0;

  };


  //! Even-odd preconditioned Wilson-like fermion actions including derivatives
  /*! @ingroup actions
   *
   * Even-odd preconditioned like Wilson-like fermion actions
   * Here, use arrays of matter fields.
   */
  template<typename T, typename P, typename Q>
  class EvenOddPrecConstDetWilsonTypeFermAct5D : public EvenOddPrecWilsonTypeFermAct5D<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddPrecConstDetWilsonTypeFermAct5D() {}

    //! Override to produce an even-odd prec. linear operator for this action
    /*! Covariant return rule - override base class function */
    virtual EvenOddPrecConstDetLinearOperatorArray<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const = 0;

    //! Override to produce an even-odd prec. Pauli-Villars linear operator for this action
    /*! Covariant return rule - override base class function */
    virtual EvenOddPrecConstDetLinearOperatorArray<T,P,Q>* linOpPV(Handle< FermState<T,P,Q> > state) const = 0;

  };

}


#endif

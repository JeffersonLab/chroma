// -*- C++ -*-
// $Id: eoprec_logdet_wilstype_fermact_w.h,v 3.1 2006-10-19 16:01:26 edwards Exp $

/*! @file
 * @brief Even-odd const determinant Wilson-like fermact
 */

#ifndef __eoprec_logdet_wilstype_fermact_w_h__
#define __eoprec_logdet_wilstype_fermact_w_h__

#include "eoprec_wilstype_fermact_w.h"
#include "eoprec_logdet_linop.h"

namespace Chroma
{
  //-------------------------------------------------------------------------------------------
  //! Even-odd preconditioned Wilson-like fermion action, specialised to clover like (gauge dependent diagonal term with exactly known derivative) structure
  /*! @ingroup actions
   *
   * Even-odd preconditioned like Wilson-like fermion actions
   */
  template<typename T, typename P, typename Q>
  class EvenOddPrecLogDetWilsonTypeFermAct : public EvenOddPrecWilsonTypeFermAct<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddPrecLogDetWilsonTypeFermAct() {}

    //! Override to produce an even-odd prec. linear operator for this action
    /*! Covariant return rule - override base class function */
    virtual EvenOddPrecLogDetLinearOperator<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const = 0;

  };

}


#endif

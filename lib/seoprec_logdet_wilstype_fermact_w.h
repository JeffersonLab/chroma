// -*- C++ -*-
/*! @file
 * @brief Symmetric even-odd const determinant Wilson-like fermact
 */

#ifndef __seoprec_logdet_wilstype_fermact_w_h__
#define __seoprec_logdet_wilstype_fermact_w_h__

#include "seoprec_wilstype_fermact_w.h"
#include "seoprec_logdet_linop.h"

namespace Chroma
{
  //-------------------------------------------------------------------------------------------
  //! Symmetric even-odd preconditioned Wilson-like fermion action, specialised to clover like (gauge dependent diagonal term with exactly known derivative) structure
  /*! @ingroup actions
   *
   * Even-odd preconditioned like Wilson-like fermion actions
   */
  template<typename T, typename P, typename Q>
  class SymEvenOddPrecLogDetWilsonTypeFermAct : public SymEvenOddPrecWilsonTypeFermAct<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~SymEvenOddPrecLogDetWilsonTypeFermAct() {}

    //! Override to produce an even-odd prec. linear operator for this action
    /*! Covariant return rule - override base class function */
    virtual SymEvenOddPrecLogDetLinearOperator<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const = 0;

  };

}


#endif

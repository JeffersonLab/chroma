// -*- C++ -*-
/*! @file
 * @brief Symmetric even-odd preconditioned Wilson-like fermion actions
 */

#ifndef __seoprec_wilstype_fermact_w_h__
#define __seoprec_wilstype_fermact_w_h__

#include "wilstype_fermact_w.h"
#include "seoprec_linop.h"

namespace Chroma
{
  //-------------------------------------------------------------------------------------------
  //! Symmetric even-odd preconditioned Wilson-like fermion actions including derivatives
  /*! @ingroup actions
   *
   * Even-odd preconditioned like Wilson-like fermion actions
   */
  template<typename T, typename P, typename Q>
  class SymEvenOddPrecWilsonTypeFermAct : public WilsonTypeFermAct<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~SymEvenOddPrecWilsonTypeFermAct() {}

    //! Override to produce an even-odd prec. linear operator for this action
    /*! Covariant return rule - override base class function */
    virtual SymEvenOddPrecLinearOperator<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const = 0;

    //! Return quark prop solver, solution of unpreconditioned system
    /*! Default implementation provided */
    virtual SystemSolver<T>* qprop(Handle< FermState<T,P,Q> > state,
				   const GroupXML_t& invParam) const;
  };




}


#endif

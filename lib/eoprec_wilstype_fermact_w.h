// -*- C++ -*-
// $Id: eoprec_wilstype_fermact_w.h,v 3.1 2006-10-19 16:01:26 edwards Exp $

/*! @file
 * @brief Even-odd preconditioned Wilson-like fermion actions
 */

#ifndef __eoprec_wilstype_fermact_w_h__
#define __eoprec_wilstype_fermact_w_h__

#include "wilstype_fermact_w.h"
#include "eoprec_linop.h"

namespace Chroma
{
  //-------------------------------------------------------------------------------------------
  //! Even-odd preconditioned Wilson-like fermion actions including derivatives
  /*! @ingroup actions
   *
   * Even-odd preconditioned like Wilson-like fermion actions
   */
  template<typename T, typename P, typename Q>
  class EvenOddPrecWilsonTypeFermAct : public WilsonTypeFermAct<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddPrecWilsonTypeFermAct() {}

    //! Override to produce an even-odd prec. linear operator for this action
    /*! Covariant return rule - override base class function */
    virtual EvenOddPrecLinearOperator<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const = 0;

    //! Return quark prop solver, solution of unpreconditioned system
    /*! Default implementation provided */
    virtual SystemSolver<T>* qprop(Handle< FermState<T,P,Q> > state,
				   const GroupXML_t& invParam) const;
  };



  //-------------------------------------------------------------------------------------------
  //! Even-odd preconditioned Wilson-like fermion actions including derivatives
  /*! @ingroup actions
   *
   * Even-odd preconditioned like Wilson-like fermion actions
   * Here, use arrays of matter fields.
   */
  template<typename T, typename P, typename Q>
  class EvenOddPrecWilsonTypeFermAct5D : public WilsonTypeFermAct5D<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddPrecWilsonTypeFermAct5D() {}

    //! Override to produce an even-odd prec. linear operator for this action
    /*! Covariant return rule - override base class function */
    virtual EvenOddPrecLinearOperatorArray<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const = 0;

    //! Override to produce an even-odd prec. Pauli-Villars linear operator for this action
    /*! Covariant return rule - override base class function */
    virtual EvenOddPrecLinearOperatorArray<T,P,Q>* linOpPV(Handle< FermState<T,P,Q> > state) const = 0;

    //! Return quark prop solver, solution of unpreconditioned system
    /*! Default implementation provided */
    virtual SystemSolverArray<T>* qpropT(Handle< FermState<T,P,Q> > state,
					 const GroupXML_t& invParam) const;
  };

}


#endif

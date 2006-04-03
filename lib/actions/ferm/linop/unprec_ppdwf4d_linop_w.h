// -*- C++ -*-
// $Id: unprec_ppdwf4d_linop_w.h,v 3.0 2006-04-03 04:58:52 edwards Exp $
/*! \file
 *  \brief Unpreconditioned projected DWF operator to 4D using prec 5D bits
 */

#ifndef __unprec_ppdwf4d_linop_w_h__
#define __unprec_ppdwf4d_linop_w_h__

#include "linearop.h"
#include "eo_prec_linop.h"
#include "handle.h"
#include "invtype.h"



namespace Chroma
{
  //! Unpreconditioned projected DWF operator to 4D, using prec. 5D pieces
  /*!
   * \ingroup linop
   */
  template<typename T, typename P, typename Q>
  class UnprecPPDWF4DLinOp : public LinearOperator<T>
  {
  public:
    //! Initialize pointer with existing pointer
    /*! Requires that the pointer p is a return value of new */
    UnprecPPDWF4DLinOp(EvenOddPrecLinearOperatorArray<T,P,Q>* D_, 
		       EvenOddPrecLinearOperatorArray<T,P,Q>* PV_,
		       const InvertParam_t& invParam_) : 
      D(D_), PV(PV_), invParam(invParam_) {}

    //! Copy pointer (one more owner)
    UnprecPPDWF4DLinOp(Handle< EvenOddPrecLinearOperatorArray<T,P,Q> > D_, 
		       Handle< EvenOddPrecLinearOperatorArray<T,P,Q> > PV_,
		       const InvertParam_t& invParam_) : 
      D(D_), PV(PV_), invParam(invParam_) {}

    //! Destructor
    ~UnprecPPDWF4DLinOp() {}

    //! Length of internal 5D
    int size() const {return D->size();}

    //! Operator lives on the entire lattice
    inline const OrderedSubset& subset() const {return all;}

    //! Apply the operator onto a source vector
    /*! For this operator, the sign is ignored */
    void operator() (T& chi, const T& psi, enum PlusMinus isign) const;

  private:
    Handle< EvenOddPrecLinearOperatorArray<T,P,Q> > D;
    Handle< EvenOddPrecLinearOperatorArray<T,P,Q> > PV;
    const InvertParam_t invParam;
  };

}


#endif

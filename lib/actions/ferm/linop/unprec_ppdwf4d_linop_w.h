// -*- C++ -*-
// $Id: unprec_ppdwf4d_linop_w.h,v 1.5 2005-01-20 03:07:34 edwards Exp $
/*! \file
 *  \brief Unpreconditioned projected DWF operator to 4D using prec 5D bits
 */

#ifndef __unprec_ppdwf4d_linop_w_h__
#define __unprec_ppdwf4d_linop_w_h__

#include "linearop.h"
#include "handle.h"
#include "invtype.h"



namespace Chroma
{
  //! Unpreconditioned projected DWF operator to 4D, using prec. 5D pieces
  /*!
   * \ingroup linop
   */
  template<typename T, typename P>
  class UnprecPPDWF4DLinOp : public LinearOperator<T>
  {
  public:
    //! Initialize pointer with existing pointer
    /*! Requires that the pointer p is a return value of new */
    UnprecPPDWF4DLinOp(const EvenOddPrecLinearOperator< multi1d<T>, P >* D_, 
		       const EvenOddPrecLinearOperator< multi1d<T>, P >* PV_,
		       const InvertParam_t& invParam_) : 
      D(D_), PV(PV_), invParam(invParam_) {}

    //! Copy pointer (one more owner)
    UnprecPPDWF4DLinOp(Handle<const EvenOddPrecLinearOperator< multi1d<T>, P > > D_, 
		       Handle<const EvenOddPrecLinearOperator< multi1d<T>, P > > PV_,
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
    const Handle< const EvenOddPrecLinearOperator< multi1d<T>, P > > D;
    const Handle< const EvenOddPrecLinearOperator< multi1d<T>, P > > PV;
    const InvertParam_t& invParam;
  };

}


#endif

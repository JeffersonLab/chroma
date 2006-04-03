// -*- C++ -*-
// $Id: unprec_dwf4d_linop_w.h,v 3.0 2006-04-03 04:58:51 edwards Exp $
/*! \file
 *  \brief Unpreconditioned projected DWF operator to 4D
 */

#ifndef __unprec_dwf4d_linop_w_h__
#define __unprec_dwf4d_linop_w_h__

#include "linearop.h"
#include "handle.h"
#include "invtype.h"



namespace Chroma
{
  //! Unpreconditioned projected DWF operator to 4D
  /*!
   * \ingroup linop
   */
  template<typename T>
  class UnprecDWF4DLinOp : public LinearOperator<T>
  {
  public:
    //! Initialize pointer with existing pointer
    /*! Requires that the pointer p is a return value of new */
    UnprecDWF4DLinOp(LinearOperatorArray<T>* D_, 
		     LinearOperatorArray<T>* PV_,
		     const InvertParam_t& invParam_) : 
      D(D_), PV(PV_), invParam(invParam_) {}

    //! Copy pointer (one more owner)
    UnprecDWF4DLinOp(Handle< LinearOperatorArray<T> > D_, 
		     Handle< LinearOperatorArray<T> > PV_,
		     const InvertParam_t& invParam_) : 
      D(D_), PV(PV_), invParam(invParam_) {}

    //! Destructor
    ~UnprecDWF4DLinOp() {}

    //! Length of internal 5D
    int size() const {return D->size();}

    //! Operator lives on the entire lattice
    inline const OrderedSubset& subset() const {return all;}

    //! Apply the operator onto a source vector
    /*! For this operator, the sign is ignored */
    void operator() (T& chi, const T& psi, enum PlusMinus isign) const;

  private:
    Handle< LinearOperatorArray<T> > D;
    Handle< LinearOperatorArray<T> > PV;
    const InvertParam_t& invParam;
  };

}


#endif

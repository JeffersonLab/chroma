// -*- C++ -*-
// $Id: unprec_dwf4d_linop_w.h,v 1.1 2004-11-08 05:41:29 edwards Exp $
/*! \file
 *  \brief Unpreconditioned projected DWF operator to 4D
 */

#ifndef __unprec_dwf4d_linop_w_h__
#define __unprec_dwf4d_linop_w_h__

#include "linearop.h"
#include "handle.h"
#include "invtype.h"

using namespace QDP;
using namespace Chroma;


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
    UnprecDWF4DLinOp(const LinearOperator< multi1d<T> >* D_, 
		     const LinearOperator< multi1d<T> >* PV_,
		     const InvertParam_t& invParam_) : 
      D(D_), PV(PV_), invParam(invParam_) {}

    //! Copy pointer (one more owner)
    UnprecDWF4DLinOp(Handle<const LinearOperator< multi1d<T> > > D_, 
		     Handle<const LinearOperator< multi1d<T> > > PV_,
		     const InvertParam_t& invParam) : 
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
    const Handle< const LinearOperator< multi1d<T> > > D;
    const Handle< const LinearOperator< multi1d<T> > > PV;
    const InvertParam_t& invParam;
  };

}

using namespace Chroma;

#endif

// -*- C++ -*-
// $Id: unprec_pdwf4d_linop_w.h,v 3.3 2007-02-22 21:11:47 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned projected DWF operator to 4D using prec 5D bits
 */

#ifndef __unprec_pdwf4d_linop_w_h__
#define __unprec_pdwf4d_linop_w_h__

#include "linearop.h"
#include "eoprec_linop.h"
#include "handle.h"
#include "actions/ferm/invert/syssolver_cg_params.h"


namespace Chroma
{
  //! Unpreconditioned projected DWF operator to 4D, using prec. 5D pieces
  /*!
   * \ingroup linop
   */
  template<typename T, typename P, typename Q>
  class UnprecPDWF4DLinOp : public LinearOperator<T>
  {
  public:
    //! Initialize pointer with existing pointer
    /*! Requires that the pointer p is a return value of new */
    UnprecPDWF4DLinOp(EvenOddPrecLinearOperatorArray<T,P,Q>* D_, 
		      UnprecLinearOperatorArray<T,P,Q>* PV_,
		      const SysSolverCGParams& invParam_) : 
      D(D_), PV(PV_), invParam(invParam_) {}

    //! Copy pointer (one more owner)
    UnprecPDWF4DLinOp(Handle< EvenOddPrecLinearOperatorArray<T,P,Q> > D_, 
		      Handle< UnprecLinearOperatorArray<T,P,Q> > PV_,
		      const SysSolverCGParams& invParam_) : 
      D(D_), PV(PV_), invParam(invParam_) {}

    //! Destructor
    ~UnprecPDWF4DLinOp() {}

    //! Length of internal 5D
    int size() const {return D->size();}

    //! Operator lives on the entire lattice
    inline const Subset& subset() const {return all;}

    //! Apply the operator onto a source vector
    /*! For this operator, the sign is ignored */
    void operator() (T& chi, const T& psi, enum PlusMinus isign) const;

  protected:
    //! Hide default constructor
    UnprecPDWF4DLinOp() {}

  private:
    Handle< EvenOddPrecLinearOperatorArray<T,P,Q> > D;
    Handle< UnprecLinearOperatorArray<T,P,Q> > PV;
    SysSolverCGParams invParam;
  };

}


#endif

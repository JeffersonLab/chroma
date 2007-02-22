// -*- C++ -*-
// $Id: tprec_linop.h,v 3.2 2007-02-22 21:11:45 bjoo Exp $
/*! @file
 * @brief Time-preconditioned Linear Operators
 */

#ifndef __tprec_linop_h__
#define __tprec_linop_h__

#include "linearop.h"

namespace Chroma
{

  //-----------------------------------------------------------------------------------
  //! Time preconditioned linear operator
  /*! @ingroup linop
   *
   * Support for time preconditioned linear operators
   * Given a matrix M written in block form:
   *
   *  M = D_t  +  D_s
   *
   * The preconditioning consists of multiplying by the inverse
   * of the time operator
   *
   * This class is used to implement the resulting linear operator
   *
   *      M'     =  1 +  D_t^(-1)*D_s
   *
   *      M'^dag =  1 +  D_s^dag * (D_t^(-1))^dag
   *
   * The non-symmetrical nature of the daggered version means the two
   * cases (no-dagger and dagger) must be handled separately. This is
   * in contrast to the standard (4D) even-odd precond. case where the
   * the daggered version has the same structure, except the dagger
   * is pushed down into the individual pieces.
   */

  template<typename T, typename P, typename Q>
  class TimePrecLinearOperator : public DiffLinearOperator<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~TimePrecLinearOperator() {}

    //! Defined on the entire lattice
    const Subset& subset() const {return all;}

    //! Return the fermion BC object for this linear operator
    virtual const FermBC<T,P,Q>& getFermBC() const = 0;

    //! The time direction
    int tDir() const = 0;

    //! Apply the time block onto a source vector
    /*! This does not need to be optimized */
    virtual void timeLinOp(T& chi, const T& psi, 
			   enum PlusMinus isign) const = 0;
  
    //! Apply the inverse of the time block onto a source vector
    virtual void timeInvLinOp(T& chi, const T& psi, 
			      enum PlusMinus isign) const = 0;
  
    //! Apply the the space block onto a source vector
    virtual void spaceLinOp(T& chi, const T& psi, 
			    enum PlusMinus isign) const = 0;

    //! Apply the operator onto a source vector
    virtual void operator() (T& chi, const T& psi, 
			     enum PlusMinus isign) const
    {
      T   tmp1, tmp2; moveToFastMemoryHint(tmp1); moveToFastMemoryHint(tmp2);

      switch (isign)
      {
      case PLUS:
	//  chi   =  psi  +  D_t^(-1)*D_s*psi
	spaceLinOp(tmp1, psi, isign);
	timeInvLinOp(tmp2, tmp1, isign);
	chi = psi + tmp2;
	break;

      case MINUS:
	//  chi   =  psi  +  D_s^dag*D_t^(-1)^dag*psi
	timeInvLinOp(tmp1, psi, isign);
	spaceLinOp(tmp2, tmp1, isign);
	chi = psi + tmp2;
	break;

      default:
	QDPIO::cerr << "unknown sign" << endl;
	QDP_abort(1);
      }
    }

    //! Apply the UNPRECONDITIONED operator onto a source vector
    /*! Mainly intended for debugging */
    virtual void unprecLinOp(T& chi, const T& psi, 
			     enum PlusMinus isign) const
    {
      T   tmp1, tmp2;  moveToFastMemoryHint(tmp1); moveToFastMemoryHint(tmp2);

      //  chi   =  D_t*psi  +  D_s*psi
      timeLinOp(tmp1, psi, isign);
      spaceLinOp(tmp2, psi, isign);
      chi = tmp1 + tmp2;
    }

    //! Apply the even-even block onto a source vector
    virtual void derivTimeLinOp(P& ds_u, const T& chi, const T& psi, 
				enum PlusMinus isign) const = 0;
  
    //! Apply the space block onto a source vector
    virtual void derivSpaceLinOp(P& ds_u, const T& chi, const T& psi, 
				 enum PlusMinus isign) const = 0;
 
    //! Apply the derivative of the operator onto a source vector
    /*! User should make sure deriv routines do a resize  */
    virtual void deriv(P& ds_u, const T& chi, const T& psi, 
		       enum PlusMinus isign) const = 0;

    //! Apply the derivative of the UNPRECONDITIONED operator onto a source vector
    /*! Mainly intended for debugging */
    virtual void derivUnprecLinOp(P& ds_u, const T& chi, const T& psi, 
				  enum PlusMinus isign) const = 0;
  };

}

#endif

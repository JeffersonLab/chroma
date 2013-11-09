// -*- C++ -*-
// $Id: tprec_logdet_linop.h,v 3.2 2007-02-22 21:11:45 bjoo Exp $
/*! @file
 * @brief Time-preconditioned Linear Operators
 */

#ifndef __tprec_logdet_linop_h__
#define __tprec_logdet_linop_h__

#include "tprec_linop.h"

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
   *
   * Here we assume D_t does depend on the gauge fields, and that we can 
   * exactly simulate the determinant without pseudofermion fields.
   * We know det D_t and can write it in the action as 
   * exp( log det D_t ) = exp( Tr Ln D_t )
   * 
   * Since we can explicitly evaluate Tr Ln D_t for the action, we
   * can also evaluate the force contribution
   * 
   * d/dt ( Tr Ln D_t ) = Tr ( (d/dt)D_t D_t^{-1} )
   *
   * and  d/dt ( D_t^{-1} ) = D_t^{-1} [ (d/dt)D_t ] D_t^{-1}
   * 
   * hence we have functions  logDetTimeLinOp()
   *  and                     derivTimeLinOp()
   *  and                     derivLogDetTimeLinOp()
   */

  template<typename T, typename P, typename Q>
  class TimePrecLogDetLinearOperator : public TimePrecLinearOperator<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~TimePrecLogDetLinearOperator() {}

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

    //! Apply the even-even block onto a source vector
    virtual void derivTimeLinOp(P& ds_u, const T& chi, const T& psi, 
				enum PlusMinus isign) const
    {
      QDPIO::cerr << "derivTime: not implemented" << endl;
      QDP_abort(1);
    }
  
    //! Apply the space block onto a source vector
    virtual void derivSpaceLinOp(P& ds_u, const T& chi, const T& psi, 
				   enum PlusMinus isign) const
    {
      QDPIO::cerr << "derivSpace: not implemented" << endl;
      QDP_abort(1);
    }
 
    //! Apply the derivative of the operator onto a source vector
    /*! User should make sure deriv routines do a resize  */
    virtual void deriv(P& ds_u, const T& chi, const T& psi, 
		       enum PlusMinus isign) const
    {
      // Need deriv of  chi^dag*[psi  +  D_t^(-1)*D_s*psi]
      enum PlusMinus msign = (isign == PLUS) ? MINUS : PLUS;

      //
      // Make sure the deriv routines do a resize !!!
      //
      T   tmp1, tmp2, tmp3;  moveToFastMemoryHint(tmp1); moveToFastMemoryHint(tmp2);
      P   ds_1;    moveToFastMemoryHint(ds_1);

      switch (isign)
      {
      case PLUS:
	//  ds_u   =  chi^dag * D_t^(-1) * D'_s * psi
	timeInvLinOp(tmp1, chi, msign);
	derivSpaceLinOp(ds_u, tmp1, psi, isign);

	//  ds_u  -=  chi^dag * D_t^(-1) * D'_t * D_t^(-1) * D_s * psi
	spaceLinOp(tmp2, psi, isign);
	timeInvLinOp(tmp3, tmp2, isign);
	derivTimeLinOp(ds_1, tmp1, tmp3, isign);
	ds_u -= ds_1;
	break;

      case MINUS:
	//  ds_u   =  chi^dag * D'_s^dag * D_t^(-1)^dag * psi
	timeInvLinOp(tmp1, psi, isign);
	derivSpaceLinOp(ds_u, chi, tmp1, isign);

	//  ds_u  -=  chi^dag * D_s^dag * D_t^(-1)^dag * D'_t^dag * D_t^(-1)^dag * psi
	spaceLinOp(tmp2, chi, msign);
	timeInvLinOp(tmp3, tmp2, msign);
	derivTimeLinOp(ds_1, tmp3, tmp1, isign);
	ds_u -= ds_1;
	break;

      default:
	QDPIO::cerr << "unknown sign" << endl;
	QDP_abort(1);
      }
    }

    //! Apply the derivative of the UNPRECONDITIONED operator onto a source vector
    /*! Mainly intended for debugging */
    virtual void derivUnprecLinOp(P& ds_u, const T& chi, const T& psi, 
				  enum PlusMinus isign) const
    {
      P   ds_tmp;  // deriv routines should resize

      //  ds_u = chi^dag * D'_t*psi  +  chi^dag * D'_s * psi

      //  ds_u  =  chi^dag * D'_t * psi
      derivTimeLinOp(ds_u, chi, psi, isign);

      //  ds_u +=  chi^dag * D'_s * psi
      derivSpaceLinOp(ds_tmp, chi, psi, isign);
      ds_u += ds_tmp;
    }

    //! Get the force from the Time Trace Log
    virtual void derivLogDetTimeLinOp(P& ds_u, enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }

    //! Get the log det of the even even part
    virtual Double logDetTimeLinOp(void) const = 0;
  };

}

#endif

// -*- C++ -*-
// $Id: teoprec_logdet_linop.h,v 3.2 2007-02-22 21:11:45 bjoo Exp $
/*! @file
 * @brief Time-preconditioned Linear Operators
 */

#ifndef __teoprec_logdet_linop_h__
#define __teoprec_logdet_linop_h__

#include "teoprec_linop.h"

namespace Chroma
{

  //-----------------------------------------------------------------------------------
  //! Even-odd and time preconditioned linear operator
  /*! @ingroup linop
   *
   * Given a matrix M written in separate time and space operators
   *
   *  M = D_t  +  D_s
   *
   * The time preconditioning consists of multiplying by the inverse
   * of the time operator to give
   *
   *     M'  =  1 +  D_t^(-1)*D_s
   *
   * For spatial even-odd precond., this interface requires the D_t to
   * be diagonal in spatial components - have no space-space or space-time
   * components. Hence, this means no even-odd component. 
   * Also, the D_s must no time coupling, and no even-even or odd-odd
   * component - e.g. all diagonal terms have been pushed into D_t.
   *
   * Rewrite M' into block form:
   *
   *       [      A             D        ]
   *       [       E,E           E,O     ]
   *  M' = [                             ]
   *       [      D             A        ]
   *       [       O,E           O,O     ]
   *
   * where the E,O refer to only the spatial dimensions and
   *
   *  A(e,e) =  1 + D_t^(-1)(e,e) * D_s(e,e) -> 1   # D_s(e,e)=0
   *  D(e,o) =  D_t^(-1)(e,e) * D_s(e,o)            # is a matrix in time-time
   *  D(o,e) =  D_t^(-1)(o,o) * D_s(o,e)            # is a matrix in time-time
   *  A(o,o) =  1 + D_t^(-1)(o,o) * D_s(o,o) -> 1   # D_s(o,o)=0
   *
   *  A(e,e)^(-1) = <diagonal in space, matrix in time-time components>
   *
   * Spatial even-odd preconditioning consists of using the triangular matrices
   *
   *      [      1              0        ]
   *      [       E,E            E,O     ]
   *  L = [                              ]
   *      [     D     A^(-1)    1        ]
   *      [      O,E   E,E        O,O    ]
   *
   * and
   *
   *      [      A              D       ]
   *      [       E,E            E,O    ]
   *  U = [                             ]
   *      [      0              1       ]
   *      [       O,E            O,O    ]
   *
   * The preconditioned matrix is formed from
   *
   *  ~
   *  M   =  L^-1 * M' * U^-1
   *
   * where
   *
   *           [      1              0        ]
   *           [       E,E            E,O     ]
   *  L^(-1) = [                              ]
   *           [   - D     A^(-1)    1        ]
   *           [      O,E   E,E        O,O    ]
   *
   * and
   *
   *           [      A^(-1)       - A^(-1) D       ]
   *           [       E,E            E,E    E,O    ]
   *  U^(-1) = [                                    ]
   *           [      0                1            ]
   *           [       O,E              O,O         ]
   *
   * Resulting in a new  M
   *
   *      [      1                    0                      ]
   *  ~   [       E,E                  E,O                   ]
   *  M = [                                                  ]
   *      [      0                A     -  D    A^(-1)  D    ]
   *      [       O,E              O,O      O,E   E,E    E,O ]
   *
   *
   * This class is used to implement the resulting linear operator
   *
   *  ~
   *  M  =  A(o,o) - D(o,e) * A^-1(e,e) * D(e,o)
   *     =  1   -  D_t^(-1)(o,o) * D_s(o,e) * D_t^(-1)(e,e) * D_s(e,o)
   *
   *  ~
   *  M^dag  = 1  -  D_s(o,e)^dag * D_t^(-1)(e,e)^dag * D_s(e,o)^dag * D_t^(-1)(o,o)^dag
   *
   * By construction, the linear operator is ONLY defined on the odd subset
   *
   * The non-symmetrical nature of the daggered version means the two
   * cases (no-dagger and dagger) must be handled separately. This is
   * in contrast to the standard (4D) even-odd precond. case where the
   * the daggered version has the same structure, except the dagger
   * is pushed down into the individual pieces.
   *
   * Here A^{-1}_{ee} does depend on the gauge fields but we can 
   * exactly simulate the determinant without pseudofermion fields
   * since we know det A_{ee} and can write it in the action as 
   * exp( log det A ) = exp( Tr Ln A )
   * 
   * Since we can explicitly evaluate Tr Ln A for the action, we
   * can also evaluate the force contribution
   * 
   * d/dt ( Tr Ln A ) = Tr ( (d/dt)A A^{-1} )
   *
   * and  d/dt ( A^{-1} ) = A^{-1} [ (d/dt)A ] A^{-1}
   * 
   * hence we have functions  logDetEvenEvenLinOp()
   *  and                     derivEvenEvenLinOp()
   *  and                     derivLogDetEvenEvenLinOp()
   */

  template<typename T, typename P, typename Q>
  class EvenOddTimePrecLogDetLinearOperator : public EvenOddTimePrecLinearOperator<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddTimePrecLogDetLinearOperator() {}

    //! Only defined on the odd lattice
    const Subset& subset() const {return rb[1];}   // not correct, need space-rb

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
      T   tmp1; moveToFastMemoryHint(tmp1);
      T   tmp2; moveToFastMemoryHint(tmp2);
      T   tmp3; moveToFastMemoryHint(tmp3);
      P   ds_1; moveToFastMemoryHint(ds_1);

      switch (isign)
       {
      case PLUS:
	// Need deriv of  - chi^dag*(D_t^(-1)(o,o) * D_s(o,e) * D_t^(-1)(e,e) * D_s(e,o))*psi
	//
	// NOTE: even with even-odd decomposition, the ds_u will still have contributions
	// on all cb. So, no adding of ds_1 onto ds_u under a subset
	//
	//  ds_u  =  + chi^dag * Dtinv(o,o) * Dt'(o,o) * Dtinv(o,o)*D_s(o,e)*Dtinv(e,e)*D_s(e,o) * psi
	evenOddSpaceLinOp(tmp1, psi, isign);
	evenEvenTimeInvLinOp(tmp2, tmp1, isign);
	oddEvenSpaceLinOp(tmp1, tmp2, isign);
	oddOddTimeInvLinOp(tmp3, chi, msign);
	derivOddOddTimeLinOp(ds_u, tmp3, tmp1, isign);

	//  ds_u  -=  chi^dag * Dtinv(o,o) * D'_s(o,e) *Dtinv(e,e)*D_s(e,o) * psi
	evenOddSpaceLinOp(tmp1, psi, isign);
	evenEvenTimeInvLinOp(tmp2, tmp1, isign);
	oddOddTimeInvLinOp(tmp3, chi, msign);
	derivOddEvenTimeLinOp(ds_1, tmp3, tmp1, isign);
	ds_u -= ds_1;

	//  ds_u  +=  chi^dag * Dtinv(o,o) * D_s(o,e) *Dtinv(e,e) * Dt'(e,e) * Dtinv(e,e)*D_s(e,o) * psi
	evenOddSpaceLinOp(tmp2, psi, isign);
	evenEvenTimeInvLinOp(tmp1, tmp2, msign);
	oddOddTimeInvLinOp(tmp3, chi, msign);
	evenOddSpaceLinOp(tmp2, tmp3, msign);
	evenEvenTimeInvLinOp(tmp3, tmp2, msign);
	derivEvenEvenTimeLinOp(ds_1, tmp3, tmp1, isign);
	ds_u += ds_1;

	//  ds_u  -=  chi^dag * Dtinv(o,o) * D_s(o,e) *Dtinv(e,e) * D'_s(e,o) * psi
	oddOddTimeInvLinOp(tmp3, chi, msign);
	evenOddSpaceLinOp(tmp2, tmp3, msign);
	evenEvenTimeInvLinOp(tmp3, tmp2, msign);
	derivEvenEvenTimeLinOp(ds_1, tmp3, psi, isign);
	ds_u -= ds_1;
	break;

      case MINUS:
	// Need deriv of  
	//   - chi^dag*(D_s(o,e)^dag * D_t^(-1)(e,e)^dag * D_s(e,o)^dag * D_t^(-1)(o,o)^dag)*psi
	//
	//  ds_u   =  chi^dag*D_s(o,e)^dag* Dtinv(e,e)^dag * Dt'(e,e)^dag *Dtinv(e,e)^dag*D_s(e,o)^dag*Dtinv(o,o)^dag*psi
	oddOddTimeInvLinOp(tmp1, psi, isign);
	evenOddSpaceLinOp(tmp2, tmp1, isign);
	evenEvenTimeInvLinOp(tmp1, tmp2, isign);
	evenOddSpaceLinOp(tmp2, chi, msign);
	evenEvenTimeInvLinOp(tmp3, tmp2, msign);
	derivEvenEvenTimeLinOp(ds_u, tmp3, tmp1, isign);

	//  ds_u  -=  chi^dag * D'_s(o,e)^dag * Dtinv(e,e)^dag*D_s(e,o)^dag*Dtinv(o,o)^dag*psi
	oddOddTimeInvLinOp(tmp1, psi, isign);
	evenOddSpaceLinOp(tmp2, tmp1, isign);
	evenEvenTimeInvLinOp(tmp1, tmp2, isign);
	derivOddEvenTimeLinOp(ds_1, chi, tmp1, isign);
	ds_u -= ds_1;

	//  ds_u  +=  chi^dag*D_s(o,e)^dag*Dtinv(e,e)^dag*D_s(e,o)^dag*Dtinv(o,o)^dag * Dt'(o,o)^dag *Dtinv(o,o)^dag*psi
	oddOddTimeInvLinOp(tmp1, psi, isign);
	evenOddSpaceLinOp(tmp2, chi, msign);
	evenEvenTimeInvLinOp(tmp3, tmp2, msign);
	oddEvenSpaceLinOp(tmp2, tmp3, msign);
	oddOddTimeInvLinOp(tmp3, tmp2, msign);
	derivOddOddTimeLinOp(ds_1, tmp3, tmp1, isign);
	ds_u += ds_1;

	//  ds_u  -=  chi^dag*D_s(o,e)^dag* Dtinv(e,e)^dag *D'_s(e,o)^dag*Dtinv(o,o)^dag*psi
	oddOddTimeInvLinOp(tmp1, psi, isign);
	evenOddSpaceLinOp(tmp2, chi, msign);
	evenEvenTimeInvLinOp(tmp3, tmp2, msign);
	derivEvenOddTimeLinOp(ds_1, tmp3, tmp1, isign);
	ds_u -= ds_1;
	break;

      default:
	QDPIO::cerr << "unknown sign" << endl;
	QDP_abort(1);
      }
    }

    //! Apply the even-even block onto a source vector
    virtual void derivEvenEvenTimeLinOp(P& ds_u, const T& chi, const T& psi, 
					enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenEvenTime: not implemented" << endl;
      QDP_abort(1);
    }
  
    //! Apply the the even-odd block onto a source vector
    virtual void derivEvenOddSpaceLinOp(P& ds_u, const T& chi, const T& psi, 
					enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOddSpace: not implemented" << endl;
      QDP_abort(1);
    }
 
    //! Apply the the odd-even block onto a source vector
    virtual void derivOddEvenSpaceLinOp(P& ds_u, const T& chi, const T& psi, 
					enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOddSpace: not implemented" << endl;
      QDP_abort(1);
    }

    //! Apply the the odd-odd block onto a source vector
    virtual void derivOddOddTimeLinOp(P& ds_u, const T& chi, const T& psi, 
				      enum PlusMinus isign) const
    {
      QDPIO::cerr << "OddOddTime: not implemented" << endl;
      QDP_abort(1);
    }
  };


}



#endif

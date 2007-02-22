// -*- C++ -*-
// $Id: teoprec_linop.h,v 3.2 2007-02-22 21:11:45 bjoo Exp $
/*! @file
 * @brief Even-odd Time-preconditioned Linear Operators
 */

#ifndef __teoprec_linop_h__
#define __teoprec_linop_h__

#include "linearop.h"

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
   */

  template<typename T, typename P, typename Q>
  class EvenOddTimePrecLinearOperator : public DiffLinearOperator<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddTimePrecLinearOperator() {}

    //! Only defined on the odd lattice
    const Subset& subset() const {return rb[1];}   // not correct, need space-rb

    //! Apply the even-even block onto a source vector
    /*! This does not need to be optimized */
    virtual void evenEvenTimeLinOp(T& chi, const T& psi, 
				   enum PlusMinus isign) const = 0;
  
    //! Apply the inverse of the even-even block onto a source vector
    virtual void evenEvenTimeInvLinOp(T& chi, const T& psi, 
				      enum PlusMinus isign) const = 0;
  
    //! Apply the the even-odd block onto a source vector
    virtual void evenOddSpaceLinOp(T& chi, const T& psi, 
				   enum PlusMinus isign) const = 0;

    //! Apply the the odd-even block onto a source vector
    virtual void oddEvenSpaceLinOp(T& chi, const T& psi, 
				   enum PlusMinus isign) const = 0;

    //! Apply the the odd-odd block onto a source vector
    virtual void oddOddTimeLinOp(T& chi, const T& psi, 
				 enum PlusMinus isign) const = 0;

    //! Apply the operator onto a source vector
    virtual void operator() (T& chi, const T& psi, 
			     enum PlusMinus isign) const
    {
      T   tmp1, tmp2; moveToFastMemoryHint(tmp1); moveToFastMemoryHint(tmp2);

      switch (isign)
      {
      case PLUS:
	//  chi   =  psi  -  D_t^(-1)(o,o)*D_s(o,e)*D_t^(-1)(e,e)*D_s(e,o)*psi
	evenOddSpaceLinOp(tmp1, psi, isign);
	evenEvenTimeInvLinOp(tmp2, tmp1, isign);
	OddEvenSpaceLinOp(tmp1, tmp2, isign);
	oddOddTimeInvLinOp(tmp2, tmp1, isign);
	chi[rb[1]] = psi - tmp2;
	break;

      case MINUS:
	//  chi   =  psi  -  D_s(o,e)^dag * D_t^(-1)(e,e)^dag * D_s(e,o)^dag * D_t^(-1)(o,o)^dag*psi
	oddOddTimeInvLinOp(tmp1, psi, isign);
	evenOddSpaceLinOp(tmp2, tmp1, isign);
	evenEvenTimeInvLinOp(tmp1, tmp2, isign);
	OddEvenSpaceLinOp(tmp2, tmp1, isign);
	chi[rb[1]] = psi - tmp2;
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
      T   tmp1, tmp2; moveToFastMemoryHint(tmp1); moveToFastMemoryHint(tmp2);

      //  chi   =  D_t*psi  +  D_s*psi
      evenEvenTimeLinOp(tmp1, psi, isign);
      evenOddSpaceLinOp(tmp2, psi, isign);
      chi[rb[0]] = tmp1 + tmp2;

      oddOddTimeLinOp(tmp1, psi, isign);
      oddEvenSpaceLinOp(tmp2, psi, isign);
      chi[rb[1]] = tmp1 + tmp2;
    }

    //! Apply the even-even block onto a source vector
    virtual void derivEvenEvenTimeLinOp(P& ds_u, const T& chi, const T& psi, 
					enum PlusMinus isign) const = 0;
  
    //! Apply the the even-odd block onto a source vector
    virtual void derivEvenOddSpaceLinOp(P& ds_u, const T& chi, const T& psi, 
					enum PlusMinus isign) const = 0;
 
    //! Apply the the odd-even block onto a source vector
    virtual void derivOddEvenSpaceLinOp(P& ds_u, const T& chi, const T& psi, 
					enum PlusMinus isign) const = 0;

    //! Apply the the odd-odd block onto a source vector
    virtual void derivOddOddTimeLinOp(P& ds_u, const T& chi, const T& psi, 
				      enum PlusMinus isign) const = 0;

    //! Apply the derivative of the operator onto a source vector
    /*! User should make sure deriv routines do a resize  */
    virtual void deriv(P& ds_u, const T& chi, const T& psi, 
		       enum PlusMinus isign) const = 0;

    //! Apply the derivative of the UNPRECONDITIONED operator onto a source vector
    /*! Mainly intended for debugging */
    virtual void derivUnprecLinOp(P& ds_u, const T& chi, const T& psi, 
				  enum PlusMinus isign) const
    {
      T   tmp1, tmp2; moveToFastMemoryHint(tmp1); moveToFastMemoryHint(tmp2);
      P   ds_tmp; moveToFastMemoryHint(ds_tmp);

      //  chi   =  D_t*psi  +  D_s*psi

      //  ds_u  =  chi^dag * (D_t)'_ee * psi
      derivEvenEvenTimeLinOp(ds_u, chi, psi, isign);

      //  ds_u +=  chi^dag * (D_s)'_eo * psi
      derivEvenOddSpaceLinOp(ds_tmp, chi, psi, isign);
      ds_u += ds_tmp;

      //  ds_u +=  chi^dag * (D_t)'_oo * psi
      derivOddOddTimeLinOp(ds_tmp, chi, psi, isign);
      ds_u += ds_tmp;

      //  ds_u +=  chi^dag * (D_s)'_oe * psi
      derivOddEvenSpaceLinOp(ds_tmp, chi, psi, isign);
      ds_u += ds_tmp;
    }

  };


}



#endif

// -*- C++ -*-
// $Id: tpreclinop.h,v 1.1 2005-02-01 21:24:05 edwards Exp $
/*! @file
 * @brief Time-preconditioned Linear Operators
 */

#ifndef __tpreclinop_h__
#define __tpreclinop_h__

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
   *      M'  =  1 +  D_t^(-1)*D_s
   *
   */

  template<typename T, typename P>
  class TimePrecLinearOperator : public DiffLinearOperator<T,P>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~TimePrecLinearOperator() {}

    //! Only defined on the odd lattice
    const OrderedSubset& subset() const {return all;}

    //! The time direction
    int tDir() const = 0;

    //! Apply the time block onto a source vector
    /*! This does not need to be optimized */
    virtual void timeLinOp(T& chi, const T& psi, 
			   enum PlusMinus isign) const = 0;
  
    //! Apply the inverse of the time block onto a source vector
    virtual void timeInvLinOp(T& chi, const T& psi, 
			      enum PlusMinus isign) const = 0;
  
    //! Apply the the spatial block onto a source vector
    virtual void spatialLinOp(T& chi, const T& psi, 
			      enum PlusMinus isign) const = 0;

    //! Apply the operator onto a source vector
    virtual void operator() (T& chi, const T& psi, 
			     enum PlusMinus isign) const
    {
      T   tmp1, tmp2;  // if an array is used here, the space is not reserved

      switch (isign)
      {
      case PLUS:
	//  chi   =  psi  +  D_t^(-1)*D_s*psi
	spatialLinOp(tmp1, psi, isign);
	timeInvLinOp(tmp2, tmp1, isign);
	chi = psi + tmp2;
	break;

      case MINUS:
	//  chi   =  psi  +  D_s^dag*D_t^(-1)^dag*psi
	timeInvLinOp(tmp1, psi, isign);
	spatialLinOp(tmp2, tmp1, isign);
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
      T   tmp1, tmp2;  // if an array is used here, the space is not reserved

      //  chi   =  D_t*psi  +  D_s*psi
      timeLinOp(tmp1, psi, isign);
      spatialLinOp(tmp2, psi, isign);
      chi = tmp1 + tmp2;
    }

    //! Apply the even-even block onto a source vector
    virtual void derivtimeLinOp(P& ds_u, const T& chi, const T& psi, 
				enum PlusMinus isign) const
    {
      QDPIO::cerr << "derivTime: not implemented" << endl;
      QDP_abort(1);
    }
  
    //! Apply the spatial block onto a source vector
    virtual void derivSpatialLinOp(P& ds_u, const T& chi, const T& psi, 
				   enum PlusMinus isign) const
    {
      QDPIO::cerr << "derivSpatial: not implemented" << endl;
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
      T   tmp1, tmp2, tmp3;  // if an array is used here, the space is not reserved
      P   ds_1;  // deriv routines should resize

      switch (isign)
      {
      case PLUS:
	//  ds_u   =  chi^dag * D_t^(-1) * D'_s * psi
	timeInvLinOp(tmp1, chi, msign);
	derivSpatialLinOp(ds_u, tmp1, psi, isign);

	//  ds_u  -=  chi^dag * D_t^(-1) * D'_t * D_t^(-1) * D_s * psi
	spatialLinOp(tmp2, psi, isign);
	timeInvLinOp(tmp3, tmp2, isign);
	derivTimeLinOp(ds_1, tmp1, tmp3, isign);
	ds_u -= ds_1;
	break;

      case MINUS:
	//  ds_u   =  chi^dag * D'_s^dag * D_t^(-1)^dag * psi
	timeInvLinOp(tmp1, psi, isign);
	derivSpatialLinOp(ds_u, chi, tmp1, isign);

	//  ds_u  -=  chi^dag * D_s^dag * D_t^(-1)^dag * D'_t^dag * D_t^(-1)^dag * psi
	spatialLinOp(tmp2, chi, msign);
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
      derivSpatialLinOp(ds_tmp, chi, psi, isign);
      ds_u += ds_tmp;
    }
  };


  //! Even-odd and time preconditioned linear operator
  /*! @ingroup linop
   *
   * Support for even-odd preconditioned linear operators
   * Given a matrix M written in block form:
   *
   *      [      A             D        ]
   *      [       E,E           E,O     ]
   *  M = [                             ]
   *      [      D             A        ]
   *      [       O,E           O,O     ]
   *
   * The preconditioning consists of using the triangular matrices
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
   *  M'   =  L^-1 * M * U^-1
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
   *      [       E,E                  E,O                   ]
   *      [                                                  ]
   *      [      0                A     -  D    A^(-1)  D    ]
   *      [       O,E              O,O      O,E   E,E    E,O ]
   *
   *
   * This class is used to implement the resulting linear operator
   *
   *      ~
   *      M  =  A(o,o) - D(o,e) . A^-1(e,e) . D(e,o)
   *
   * By construction, the linear operator is ONLY defined on the odd subset
   *
   */

  template<typename T, typename P>
  class EvenOddTimePrecLinearOperator : public DiffLinearOperator<T,P>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddTimePrecLinearOperator() {}

    //! Only defined on the odd lattice
    const OrderedSubset& subset() const {return rb[1];}

    //! Apply the even-even block onto a source vector
    /*! This does not need to be optimized */
    virtual void evenEvenLinOp(T& chi, const T& psi, 
			       enum PlusMinus isign) const = 0;
  
    //! Apply the inverse of the even-even block onto a source vector
    virtual void evenEvenInvLinOp(T& chi, const T& psi, 
				  enum PlusMinus isign) const = 0;
  
    //! Apply the the even-odd block onto a source vector
    virtual void evenOddLinOp(T& chi, const T& psi, 
			      enum PlusMinus isign) const = 0;

    //! Apply the the odd-even block onto a source vector
    virtual void oddEvenLinOp(T& chi, const T& psi, 
			      enum PlusMinus isign) const = 0;

    //! Apply the the odd-odd block onto a source vector
    virtual void oddOddLinOp(T& chi, const T& psi, 
			     enum PlusMinus isign) const = 0;

    //! Apply the operator onto a source vector
    virtual void operator() (T& chi, const T& psi, 
			     enum PlusMinus isign) const
    {
      T   tmp1, tmp2;  // if an array is used here, the space is not reserved

      /*  Tmp1   =  D     A^(-1)     D    Psi  */
      /*      O      O,E        E,E   E,O    O */
      evenOddLinOp(tmp1, psi, isign);
      evenEvenInvLinOp(tmp2, tmp1, isign);
      oddEvenLinOp(tmp1, tmp2, isign);

      /*  Chi   =  A    Psi  -  Tmp1  */
      /*     O      O,O    O        O */
      oddOddLinOp(chi, psi, isign);
      chi[rb[1]] -= tmp1;
    }

    //! Apply the UNPRECONDITIONED operator onto a source vector
    /*! Mainly intended for debugging */
    virtual void unprecLinOp(T& chi, const T& psi, 
			     enum PlusMinus isign) const
    {
      T   tmp1, tmp2;  // if an array is used here, the space is not reserved

      /*  Chi   =   A    Psi   +    D    Psi   */
      /*     E       E,E    O        E,O    O  */
      evenEvenLinOp(tmp1, psi, isign);
      evenOddLinOp(tmp2, psi, isign);
      chi[rb[0]] = tmp1 + tmp2;

      /*  Chi   =  A    Psi  -  Tmp1  */
      /*     O      O,O    O        O */
      oddEvenLinOp(tmp1, psi, isign);
      oddOddLinOp(tmp2, psi, isign);
      chi[rb[1]] = tmp1 + tmp2;
    }

    //! Apply the even-even block onto a source vector
    virtual void derivEvenEvenLinOp(P& ds_u, const T& chi, const T& psi, 
				    enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }
  
    //! Apply the the even-odd block onto a source vector
    virtual void derivEvenOddLinOp(P& ds_u, const T& chi, const T& psi, 
				   enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }
 
    //! Apply the the odd-even block onto a source vector
    virtual void derivOddEvenLinOp(P& ds_u, const T& chi, const T& psi, 
				   enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }

    //! Apply the the odd-odd block onto a source vector
    virtual void derivOddOddLinOp(P& ds_u, const T& chi, const T& psi, 
				  enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }

    //! Apply the derivative of the operator onto a source vector
    /*! User should make sure deriv routines do a resize  */
    virtual void deriv(P& ds_u, const T& chi, const T& psi, 
		       enum PlusMinus isign) const
    {
      // Need deriv of  (A_oo - D_oe*Ainv_ee*D_eo*psi_e)
      enum PlusMinus msign = (isign == PLUS) ? MINUS : PLUS;

      //
      // Make sure the deriv routines do a resize !!!
      //
      T   tmp1, tmp2, tmp3;  // if an array is used here, the space is not reserved
      P   ds_1;  // deriv routines should resize

      //
      // NOTE: even with even-odd decomposition, the ds_u will still have contributions
      // on all cb. So, no adding of ds_1 onto ds_u under a subset
      //
      //  ds_u  =  chi^dag * A'_oo * psi
      derivOddOddLinOp(ds_u, chi, psi, isign);

      //  ds_u  -=  chi^dag * D'_oe * Ainv_ee * D_eo * psi_o
      evenOddLinOp(tmp1, psi, isign);
      evenEvenInvLinOp(tmp2, tmp1, isign);
      derivOddEvenLinOp(ds_1, chi, tmp2, isign);
      ds_u -= ds_1;

      //  ds_u  +=  chi^dag * D_oe * Ainv_ee * A'_ee * Ainv_ee * D_eo * psi_o
      evenOddLinOp(tmp1, psi, isign);
      evenEvenInvLinOp(tmp2, tmp1, isign);
      evenOddLinOp(tmp1, chi, msign);
      evenEvenInvLinOp(tmp3, tmp1, msign);
      derivEvenEvenLinOp(ds_1, tmp3, tmp2, isign);
      ds_u += ds_1;

      //  ds_u  -=  chi^dag * D_oe * Ainv_ee * D'_eo * psi_o
      evenOddLinOp(tmp1, chi, msign);
      evenEvenInvLinOp(tmp3, tmp1, msign);
      derivEvenOddLinOp(ds_1, tmp3, psi, isign);
      ds_u -= ds_1;
    }

    //! Apply the derivative of the UNPRECONDITIONED operator onto a source vector
    /*! Mainly intended for debugging */
    virtual void derivUnprecLinOp(P& ds_u, const T& chi, const T& psi, 
				  enum PlusMinus isign) const
    {
      T   tmp1, tmp2;  // if an array is used here, the space is not reserved
      P   ds_tmp;  // deriv routines should resize

      //  ds_u  =  chi^dag * A'_ee * psi
      derivEvenEvenLinOp(ds_u, chi, psi, isign);

      //  ds_u +=  chi^dag * D'_eo * psi
      derivEvenOddLinOp(ds_tmp, chi, psi, isign);
      ds_u += ds_tmp;

      //  ds_u +=  chi^dag * D'_oe * psi
      derivOddEvenLinOp(ds_tmp, chi, psi, isign);
      ds_u += ds_tmp;

      //  ds_u +=  chi^dag * A'_oo * psi
      derivOddOddLinOp(ds_tmp, chi, psi, isign);
      ds_u += ds_tmp;
    }

  };


}



#endif

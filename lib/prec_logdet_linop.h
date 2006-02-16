
// -*- C++ -*-
// $Id: prec_logdet_linop.h,v 2.2 2006-02-16 02:29:44 bjoo Exp $

/*! @file
 * @brief Preconditioned  Linear Operators  where the Even Even part depends on the Gauge Field and we can evaluate Log Det E, where E is the Even Even part. Essentially this is for things like clover.
 */
#ifndef PREC_LOGDET_LINOP_H
#define PREC_LOGDET_LINOP_H

#include "eo_prec_linop.h"

using namespace QDP::Hints;

namespace Chroma {


  //! Even-odd preconditioned linear operator
  /*! @ingroup linop
   *
   * Support for even-odd preconditioned linear operators
   * Given a matrix M written in block form:
   *
   *      [      A             D        ]
   *      [       E,E           E,O     ]
   *      [                             ]
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
   *      [      1            A^{-1}  D       ]
   *      [       E,E          EE      E,O    ]
   *  U = [                                   ]
   *      [      0                    1       ]
   *      [       O,E                  O,O    ]
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
   *           [      1            - A^(-1) D       ]
   *           [       E,E            E,E    E,O    ]
   *  U^(-1) = [                                    ]
   *           [      0                1            ]
   *           [       O,E              O,O         ]
   *
   * Resulting in a new  M
   *
   *      [      A                    0                      ]
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
   * where A^{-1}_{ee} does depend on the gauge fields but we can 
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
   * hence we have functions  lnDetEvenEven()
   *  and                     derivEvenEven()
   *  and                     derivEvenEvenInv()
   *  --- this needs work, may not even need more than lnDetEvenEven
   *
   *
   */

  template<typename T, typename P>
  class EvenOddPrecLogDetLinearOperator : public EvenOddPrecLinearOperator<T,P>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddPrecLogDetLinearOperator() {}

    //! Only defined on the odd lattice
    const OrderedSubset& subset() const {return rb[1];}

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
      moveToFastMemoryHint(tmp1);
      moveToFastMemoryHint(tmp2);
      moveToFastMemoryHint(tmp3);

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

    //! Get the force from the EvenEven Trace Log
    virtual void derivEvenEvenLogDet(P& ds_u, enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }

    //! Get the log det of the even even part
    virtual Double LogDetEvenEven(void) const = 0;
  };


  //! Even-odd preconditioned 5D linear operator
  /*! @ingroup linop
   *
   * Support for even-odd preconditioned linear operators
   * Given a matrix M written in block form:
   *
   *      [      A             D        ]
   *      [       E,E           E,O     ]
   *      [                             ]
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
   *      [      1            A^{-1}  D       ]
   *      [       E,E          EE      E,O    ]
   *  U = [                                   ]
   *      [      0                    1       ]
   *      [       O,E                  O,O    ]
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
   *           [      1            - A^(-1) D       ]
   *           [       E,E            E,E    E,O    ]
   *  U^(-1) = [                                    ]
   *           [      0                1            ]
   *           [       O,E              O,O         ]
   *
   * Resulting in a new  M
   *
   *      [      A                    0                      ]
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
   * where A^{-1}_{ee} does depend on the gauge fields but we can 
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
   * hence we have functions  lnDetEvenEven()
   *  and                     derivEvenEven()
   *  and                     derivEvenEvenInv()
   *  --- this needs work, may not even need more than lnDetEvenEven
   *
   *
   */


  template<typename T, typename P>
  class EvenOddPrecLogDetLinearOperator< multi1d<T>, P > : public EvenOddPrecLinearOperator< multi1d<T>, P >
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddPrecLogDetLinearOperator() {}

    //! Only defined on the odd lattice
    const OrderedSubset& subset() const {return rb[1];}

    virtual int size() const = 0;

    //! Apply the operator onto a source vector
    /*! User should make sure deriv routines do a resize  */
    virtual void deriv(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
		       enum PlusMinus isign) const
    {
      // Need deriv of  (A_oo - D_oe*Ainv_ee*D_eo*psi_e)
      enum PlusMinus msign = (isign == PLUS) ? MINUS : PLUS;

      //
      // Make sure the deriv routines do a resize !!!
      //
      multi1d<T>   tmp1(size()), tmp2(size()), tmp3(size());  // if an array is used here, the space is not reserved
      moveToFastMemoryHint(tmp1);
      moveToFastMemoryHint(tmp2);
      moveToFastMemoryHint(tmp3);

      P            ds_1;  // deriv routines should resize

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

    //! Apply the even-even block onto a source vector
    virtual void derivEvenEvenLinOp(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
				    enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }
  
    //! Apply the the even-odd block onto a source vector
    virtual void derivEvenOddLinOp(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
				   enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }
 
    //! Apply the the odd-even block onto a source vector
    virtual void derivOddEvenLinOp(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
				   enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }

    //! Apply the the odd-odd block onto a source vector
    virtual void derivOddOddLinOp(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi, 
				  enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }

    //! Get the force from the EvenEven Trace Log
    virtual void derivEvenEvenLogDet(P& ds_u, enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }

    //! Get the log det of the even even part
    virtual Double LogDetEvenEven(void) const = 0;
  };



}
#endif

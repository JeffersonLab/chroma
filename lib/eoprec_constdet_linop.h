// -*- C++ -*-
// $Id: eoprec_constdet_linop.h,v 3.2 2007-02-22 21:11:44 bjoo Exp $
/*! @file
 * @brief Preconditioned 4D Linop with Gauge Independent Even-Even part
 */

#ifndef __eoprec_constdet_linop_h__
#define __eoprec_constdet_linop_h__

#include "eoprec_linop.h"

using namespace QDP::Hints;

namespace Chroma 
{

  //----------------------------------------------------------------
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
   * where A^{-1}_{ee} is independent of the gauge fields. This
   * means that the det A_{ee} is an irrelevant constant and that
   * the force term due to the A_{ee} part is zero.
   *
   * This structure suits most of the linear operators we use, and 
   * It simplifies the force term.
   * By construction, the linear operator is ONLY defined on the odd subset
   *
   */

  template<typename T, typename P, typename Q>
  class EvenOddPrecConstDetLinearOperator : public EvenOddPrecLinearOperator<T,P,Q>  
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddPrecConstDetLinearOperator() {}

    //! Return the fermion BC object for this linear operator
    virtual const FermBC<T,P,Q>& getFermBC() const = 0;

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
      T   tmp1, tmp2;  // if an array is used here, the space is not reserved
      moveToFastMemoryHint(tmp1);
      moveToFastMemoryHint(tmp2);

      P   ds_1;  // deriv routines should resize

      //
      // NOTE: even with even-odd decomposition, the ds_u will still have contributions
      // on all cb. So, no adding of ds_1 onto ds_u under a subset
      //
            
      //  ds_u  -=  chi^dag * D'_oe * Ainv_ee * D_eo * psi_o
      evenOddLinOp(tmp1, psi, isign);
      evenEvenInvLinOp(tmp2, tmp1, isign);
      derivOddEvenLinOp(ds_u, chi, tmp2, isign);
      

      //  ds_u  -=  chi^dag * D_oe * Ainv_ee * D'_eo * psi_o
      evenOddLinOp(tmp1, chi, msign);
      evenEvenInvLinOp(tmp2, tmp1, msign);
      derivEvenOddLinOp(ds_1, tmp2, psi, isign);
      ds_u += ds_1;

      // This is yucky and ws should have a *= function 
      // so that I wouldn't have to expose the fact that I expect
      // ds_u to be an Nd dimensional array
      for(int mu=0; mu < Nd; mu++) { 
	ds_u[mu] *= Real(-1);
      }

      getFermBC().zero(ds_u);
    }


    //! Apply the derivative of the operator onto a source vector
    /*! User should make sure deriv routines do a resize  */
    virtual void deriv(TowerArray<typename PQTraits<Q>::Base_t>& ds_u,
	       const Tower<T>& chi,
	       const Tower<T>& psi,
               const P& p,
	       enum PlusMinus isign)
    {
      // Need deriv of  (A_oo - D_oe*Ainv_ee*D_eo*psi_e)
      enum PlusMinus msign = (isign == PLUS) ? MINUS : PLUS;

      //
      // Make sure the deriv routines do a resize !!!
      //
      Tower<T>   tmp1(chi.size());
      Tower<T>   tmp2(chi.size());  // if an array is used here, the space is not reserved
    
      TowerArray<typename PQTraits<Q>::Base_t>  ds_1(Nd,chi.size());  // deriv routines should resize

      //
      // NOTE: even with even-odd decomposition, the ds_u will still have contributions
      // on all cb. So, no adding of ds_1 onto ds_u under a subset
      //
            
      //  ds_u  -=  chi^dag * D'_oe * Ainv_ee * D_eo * psi_o
      evenOddLinOp(tmp1, psi, p, isign);
      evenEvenInvLinOp(tmp2, tmp1, p, isign);
      derivOddEvenLinOp(ds_u, chi, tmp2, p, isign);
      

      //  ds_u  -=  chi^dag * D_oe * Ainv_ee * D'_eo * psi_o
      evenOddLinOp(tmp1, chi, p, msign);
      evenEvenInvLinOp(tmp2, tmp1, p, msign);
      derivEvenOddLinOp(ds_1, tmp2, psi, p, isign);
      ds_u += ds_1;

      // This is yucky and ws should have a *= function 
      // so that I wouldn't have to expose the fact that I expect
      // ds_u to be an Nd dimensional array
      for(int mu=0; mu < Nd; mu++) { 
	ds_u[mu] *= Real(-1);
      }

      getFermBC().zero(ds_u);
    }

    //! Apply the even-even block onto a source vector
    virtual void derivEvenEvenLinOp(P& ds_u, const T& chi, const T& psi, 
				    enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }

    virtual void derivEvenEvenLinOp(TowerArray<typename PQTraits<Q>::Base_t>& ds_u,
				    const Tower<T>& chi,
				    const Tower<T>& psi,
				    const P& p,
				    enum PlusMinus isign)
    {
	QDPIO::cerr << "Not Implemented" << endl;
	QDP_abort(1);
    }

    //! Apply the the even-odd block onto a source vector
    virtual void derivEvenOddLinOp(P& ds_u, const T& chi, const T& psi, 
				   enum PlusMinus isign) const
    {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }

    virtual void derivEvenOddLinOp(TowerArray<typename PQTraits<Q>::Base_t>& ds_u,
				   const Tower<T>& chi,
				   const Tower<T>& psi,
				   const P& p,
				   enum PlusMinus isign)
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

 
    virtual void derivOddEvenLinOp(TowerArray<typename PQTraits<Q>::Base_t>& ds_u,
				   const Tower<T>& chi,
				   const Tower<T>& psi,
				   const P& p,
				   enum PlusMinus isign)
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

    virtual void derivOddOddLinOp(TowerArray<typename PQTraits<Q>::Base_t>& ds_u,
				  const Tower<T>& chi,
				  const Tower<T>& psi,
				  const P& p,
				  enum PlusMinus isign)
   {
      QDPIO::cerr << "EvenOdd: not implemented" << endl;
      QDP_abort(1);
    }
  };


  //----------------------------------------------------------------
  //! Even-odd preconditioned linear operator including derivatives for arrays
  /*! @ingroup linop
   *
   * Support for even-odd preconditioned linear operators with derivatives
   *
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

  template<typename T, typename P, typename Q>
  class EvenOddPrecConstDetLinearOperatorArray : public EvenOddPrecLinearOperatorArray<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddPrecConstDetLinearOperatorArray() {}

    //! Only defined on the odd lattice
    const Subset& subset() const {return rb[1];}

    //! Expected length of array index
    virtual int size(void) const = 0;

    //! Return the fermion BC object for this linear operator
    virtual const FermBC<T,P,Q>& getFermBC() const = 0;

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
      //      moveToFastMemoryHint(tmp3);

      P            ds_1;  // deriv routines should resize

      //
      // NOTE: even with even-odd decomposition, the ds_u will still have contributions
      // on all cb. So, no adding of ds_1 onto ds_u under a subset
      //
      //  ds_u  =  chi^dag * A'_oo * psi
      //      derivOddOddLinOp(ds_u, chi, psi, isign);

      


      //  ds_u  -=  chi^dag * D'_oe * Ainv_ee * D_eo * psi_o
      evenOddLinOp(tmp1, psi, isign);
      evenEvenInvLinOp(tmp2, tmp1, isign);
      derivOddEvenLinOp(ds_u, chi, tmp2, isign);

      //  ds_u  -=  chi^dag * D_oe * Ainv_ee * D'_eo * psi_o
      evenOddLinOp(tmp1, chi, msign);
      evenEvenInvLinOp(tmp2, tmp1, msign);
      derivEvenOddLinOp(ds_1, tmp2, psi, isign);
      ds_u += ds_1;

      // This is yucky and ws should have a *= function 
      // so that I wouldn't have to expose the fact that I expect
      // ds_u to be an Nd dimensional array
      for(int mu=0; mu < Nd; mu++) { 
	ds_u[mu] *= Real(-1);
      }

      getFermBC().zero(ds_u);
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
  };




}
#endif

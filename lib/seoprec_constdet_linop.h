// -*- C++ -*-
/*! @file
 * @brief Preconditioned 4D Linop with Gauge Independent Even-Even part
 */

#ifndef __seoprec_constdet_linop_h__
#define __seoprec_constdet_linop_h__

#include "seoprec_linop.h"

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
 * The preconditioning consists of factorizing as
 *
 *
 *      M_unprec = M_diag x M'
 *
 *
 *     M_diag = [   A_E,E    0   ]
 *              [   0      A_O,O ]
 *
 *              [ 1        M_E,O ]
 *     M' =     [                ]
 *              [ M_O,E      1   ]
 *
 *
 *     M_E,O = A^{-1}_E,E D_E,O
 *     M_O,E = A^{-1}_O,O D_E,O
 *
 *  We then do a schur precondition of M'
 *
 *  as M' =  L  M  R
 *
 *
 *     L  =   [  1      0  ]
 *            [            ]
 *            [  M_O,E  1  ]
 *
 *
 *
 *     R  =   [  1    M_E,O ]
 *            [             ]
 *            [  0      1   ]
 *
 *
 *     M  =   [   1                0     ]
 *            [                          ]
 *            [   0                S     ]
 *
 * where S = 1 - ( M_O,E )( M_E,O )
 *
 *
 *    L^{-1}  =   [   1     0  ]
 *                [            ]
 *                [- M_O,E  1  ]
 *
 *
 *
 *   R^{-1}  =   [  1    - M_E,O ]
 *               [               ]
 *               [  0       1    ]
 *
 *
 *
 * For props we need to solve:
 *
 *       M  x' = L^{-1} (M_diag)^{-1} b
 *
 * and then x = R^{-1} x'  (backsubstitution)
 *
 * for HMC we have det(L)=1 det(R)=1
 * det(M_diag) is handled by TraceLog/LogDet calls to clover
 * and we pseudofermionize only the 'S' operators
 *
 *
 * Structure: capture D_oe, D_eo, A_ee, A_oo, A^{-1}_oo and A^{-1}_ee
 * as
 *    unprecEvenOddLinOp
 *    unprecOddEvenLinOp
 *    scaleEvenEvenLinOp
 *    scaleOddOddLinOp
 *    scaleOddOddInvLinOp
 *    scaleEvenEvenInvLinOp
 *
 *    then we can write category defaults for:
 *
 *    then
 *      evenOddLinOp is M_eo = A^{-1}_ee D_eo or  h.c.
 *      oddEvenLinOp is M_oe = A^{-1}_oo D_oe or  h.c.
 *
 *    the so called Jacobi Operator:
 *         [ 1             M_eo ]
 *         [ M_oe           1   ]
 *
 *    the Schur Operator:  S = 1 - M_oe M_eo
 *
 *    the Unprec Operator: diag( A_ee A_oo ) M_jacobi
 *
 *    We can also write category default for the forces in terms
 *    of
 *         derivEvenOddLinOp
 *    and  derivOddEvenLinOp
 *
 *    Const-detness or log-det ness will depend on how the
 *    derivEvenOdd and derivOddEven operators are coded

 *
 */
  template<typename T, typename P, typename Q>
  class SymEvenOddPrecConstDetLinearOperator : public SymEvenOddPrecLinearOperator<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~SymEvenOddPrecConstDetLinearOperator() {}

    //! Return the fermion BC object for this linear operator
    virtual const FermBC<T,P,Q>& getFermBC() const = 0;

    //! Apply the inverse of the even-even block onto a source std::vector
    virtual void scaleEvenEvenInvLinOp(T& chi, const T& psi,
    		enum PlusMinus isign) const override = 0;

    //! Apply the inverse of the odd-odd block onto a source std::vector
    virtual void scaleOddOddInvLinOp(T& chi, const T& psi,
    		enum PlusMinus isign) const override = 0;

    //! Apply the even-even block onto a source std::vector
    /*! This does not need to be optimized */
    virtual void scaleEvenEvenLinOp(T& chi, const T& psi,
    		enum PlusMinus isign) const override = 0;

    //! Apply the odd-odd block onto a source std::vector
    /*! This does not need to be optimized */
    virtual void scaleOddOddLinOp(T& chi, const T& psi,
    		enum PlusMinus isign) const override = 0;

    //! Apply the even-odd block onto a source std::vector
    /*! This does not need to be optimized */
    virtual void unprecEvenOddLinOp(T& chi, const T& psi,
    		enum PlusMinus isign) const override = 0;

    //! Apply the odd-even block onto a source std::vector
    /*! This does not need to be optimized */
    virtual void unprecOddEvenLinOp(T& chi, const T& psi,
    		enum PlusMinus isign) const override = 0;

    virtual void derivUnprecEvenOddLinOp(P& ds_u, const T& chi, const T& psi,
    		enum PlusMinus isign) const = 0;

    virtual void derivUnprecOddEvenLinOp(P& ds_u, const T& chi, const T& psi,
    		enum PlusMinus isign) const = 0;

    //! Apply the the even-odd block onto a source std::vector
    virtual void derivEvenOddLinOp(P& ds_u, const T& chi, const T& psi, 
				   enum PlusMinus isign) const override
    {
    	ds_u.resize(Nd);
    	ds_u = zero;

    	T tmp; moveToFastMemoryHint(tmp);
    	if( isign == PLUS ) {
    		scaleEvenEvenInvLinOp(tmp,chi,MINUS);
    		derivUnprecEvenOddLinOp(ds_u, tmp, psi, PLUS);
    	}
    	else {
    		scaleOddOddInvLinOp(tmp,psi,MINUS);
    		derivUnprecEvenOddLinOp(ds_u,chi, tmp, MINUS);

    	}
    	getFermBC().zero(ds_u);
    }
 
    //! Apply the the odd-even block onto a source std::vector
    virtual void derivOddEvenLinOp(P& ds_u, const T& chi, const T& psi, 
				   enum PlusMinus isign) const override
    {
    	ds_u.resize(Nd);
    	ds_u = zero;

    	T tmp; moveToFastMemoryHint(tmp);
    	if( isign == PLUS ) {
    		scaleOddOddInvLinOp(tmp,chi,MINUS);
    		derivUnprecOddEvenLinOp(ds_u, tmp, psi, PLUS);
    	}
    	else {
    		scaleEvenEvenLinOp(tmp,psi,MINUS);
    		derivUnprecEvenOddLinOp(ds_u,chi, tmp, MINUS);

    	}
    	getFermBC().zero(ds_u);
    }



  };



}
#endif

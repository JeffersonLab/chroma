// -*- C++ -*-
/*! @file
 * @brief Symmetric preconditioned linear pperators where the even-even and odd-odd parts depends on the gauge field.
 *
 * We assume we can evaluate Log Det E and Log Det O, where E/E is the (Even Even)[Odd Odd] part. 
 * Essentially this is for things like clover.
 */

#ifndef __seoprec_logdet_linop_h__
#define __seoprec_logdet_linop_h__

#include "seoprec_constdet_linop.h"

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
 *    unprecEvenEvenLinOp
 *    unprecOddOddLinOp
 *    unprecOddOddInvLinOp
 *    unprecEvenEvenInvLinOp
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


  // Inherit from ConstDet and extend structure -- Then in Monomials we can pass this as a constdet
  template<typename T, typename P, typename Q>
  class SymEvenOddPrecLogDetLinearOperator : public SymEvenOddPrecConstDetLinearOperator<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~SymEvenOddPrecLogDetLinearOperator() {}

    //! Only defined on the odd lattice
    const Subset& subset() const {return rb[1];}

    //! Return the fermion BC object for this linear operator
    virtual const FermBC<T,P,Q>& getFermBC() const = 0;

	//! Apply the even-even block onto a source std::vector
	/*! This does not need to be optimized */
	virtual void unprecEvenEvenLinOp(T& chi, const T& psi,
			enum PlusMinus isign) const override = 0;

	//! Apply the inverse of the even-even block onto a source std::vector
	virtual void unprecEvenEvenInvLinOp(T& chi, const T& psi,
			enum PlusMinus isign) const override = 0;

	//! Apply the the odd-odd block onto a source std::vector
	virtual void unprecOddOddLinOp(T& chi, const T& psi,
			enum PlusMinus isign)  const override  = 0;

	//! Apply the inverse of the odd-odd block onto a source std::vector
	virtual void unprecOddOddInvLinOp(T& chi, const T& psi,
			enum PlusMinus isign) const override = 0;

	//! Apply the the even-odd block onto a source std::vector
	virtual void unprecEvenOddLinOp(T& chi, const T& psi,
			enum PlusMinus isign) const override = 0;

	//! Apply the the odd-even block onto a source std::vector
	virtual void unprecOddEvenLinOp(T& chi, const T& psi,
			enum PlusMinus isign) const override = 0;

    // Deriv of A_ee
    virtual void derivUnprecEvenEvenLinOp(P& ds_u, const T& chi, const T& psi,
    		enum PlusMinus isign) const = 0;

    // Deriv of A_oo
    virtual void derivUnprecOddOddLinOp(P& ds_u, const T& chi, const T& psi,
    		enum PlusMinus isign) const = 0;

    // Deriv of  D_eo
    virtual void derivUnprecEvenOddLinOp(P& ds_u, const T& chi, const T& psi,
    		enum PlusMinus isign) const override = 0;

    // Deriv of  D_eo
    virtual void derivUnprecOddEvenLinOp(P& ds_u, const T& chi, const T& psi,
    		enum PlusMinus isign) const override = 0;

    //! deriv of  A^{-1} = - A^{-1} deriv(A) A^{-1}
    virtual void derivUnprecEvenEvenInvLinOp(P& ds_u, const T& chi, const T& psi,
				    enum PlusMinus isign) const
    {
    	enum PlusMinus msign = ( isign == PLUS ) ? MINUS : PLUS;
    	T Achi = zero;
    	T Apsi = zero;
    	unprecEvenEvenInvLinOp(Achi,chi, msign);
    	unprecEvenEvenInvLinOp(Apsi,psi, isign);
    	derivUnprecEvenEvenLinOp(ds_u, Achi,Apsi, isign);
    	for(int mu=0; mu < Nd; ++mu) {
    		ds_u[mu] *= Real(-1);
    	}
    	getFermBC().zero(ds_u);

    }

    //! deriv of  A^{-1} = - A^{-1} deriv(A) A^{-1}
    virtual void derivUnprecOddOddInvLinOp(P& ds_u, const T& chi, const T& psi,
    		enum PlusMinus isign) const
    {
    	enum PlusMinus msign = ( isign == PLUS ) ? MINUS : PLUS;
    	T Achi = zero;
    	T Apsi = zero;
    	unprecOddOddInvLinOp(Achi,chi, msign);
    	unprecOddOddInvLinOp(Apsi,psi, isign);
    	derivUnprecOddOddLinOp(ds_u, Achi,Apsi, isign);
    	for(int mu=0; mu < Nd; ++mu) ds_u[mu] *= Real(-1);
    	getFermBC().zero(ds_u);

    }

    //! Apply the the even-odd block onto a source std::vector
    //
    //  Connect chi_even through a derivative with psi_odd
    virtual void derivEvenOddLinOp(P& ds_u, const T& chi, const T& psi, 
				   enum PlusMinus isign) const override
    {
    	T tmp; moveToFastMemoryHint(tmp);

    	P ds_tmp;
    	if( isign == PLUS ) {
    		unprecEvenOddLinOp(tmp,psi, PLUS);
    		derivUnprecEvenEvenInvLinOp(ds_u, chi, tmp, PLUS);

    		unprecEvenEvenInvLinOp(tmp, chi, MINUS);
    		derivUnprecEvenOddLinOp(ds_tmp, tmp,psi, PLUS);
    		ds_u += ds_tmp;
    	}
    	else {
    		// W_o = A^{-dag}_oo Y_o
    		unprecOddOddInvLinOp(tmp, psi, MINUS);
    		// X_e^\dagger d [ D^\dagger ]_eo W_o
    		derivUnprecEvenOddLinOp(ds_u, chi, tmp, MINUS);

    		unprecOddEvenLinOp(tmp,chi,PLUS);
    		derivUnprecOddOddInvLinOp(ds_tmp, tmp,psi, MINUS);
    		ds_u += ds_tmp;
    	}
    	getFermBC().zero(ds_u);

    }
 
    //! Apply the the odd-even block onto a source std::vector
    //
    //  Connect chi_odd with psi_even
    virtual void derivOddEvenLinOp(P& ds_u, const T& chi, const T& psi, 
    		enum PlusMinus isign) const override
    {
    	T tmp; moveToFastMemoryHint(tmp);

    	P ds_tmp;
    	if( isign == PLUS ) {
    		unprecOddEvenLinOp(tmp,psi, PLUS);
    		derivUnprecOddOddInvLinOp(ds_u, chi, tmp, PLUS);

    		unprecOddOddInvLinOp(tmp, chi, MINUS);
    		derivUnprecOddEvenLinOp(ds_tmp, tmp,psi, PLUS);
    		ds_u += ds_tmp;
    	}
    	else {
    		unprecEvenEvenInvLinOp(tmp, psi, MINUS);
    		derivUnprecOddEvenLinOp(ds_u, chi, tmp, MINUS);

    		unprecEvenOddLinOp(tmp,chi,PLUS);
    		derivUnprecEvenEvenInvLinOp(ds_tmp, tmp,psi, MINUS);
    		ds_u += ds_tmp;
    	}
    	getFermBC().zero(ds_u);

    }


    //! Get the force from the EvenEven Trace Log
    virtual void derivLogDetEvenEvenLinOp(P& ds_u, enum PlusMinus isign) const = 0;

    //! Get the log det of the even even part
    virtual Double logDetEvenEvenLinOp(void) const = 0;

    //! Get the force from the OddOdd Trace Log
    virtual void derivLogDetOddOddLinOp(P& ds_u, enum PlusMinus isign) const = 0;

    //! Get the log det of the odd odd part
    virtual Double logDetOddOddLinOp(void) const = 0;
  };


}
#endif

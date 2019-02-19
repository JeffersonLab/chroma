// -*- C++ -*-
/*! @file
 * @brief Base class for symmetric even-odd preconditioned 4D and 5D Linop
 */

#ifndef __seoprec_linop_h__
#define __seoprec_linop_h__

#include "chromabase.h"
#include "linearop.h"

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
 *    ConstDet will treat A as a constant (muliplicative factor)
 *    LogDet will evaluate chain rule in derivEvenOddLinOp
 *       and also provide lodDet terms and TraceLog forces
 *
 */

template<typename T, typename P, typename Q>
class SymEvenOddPrecLinearOperator : public DiffLinearOperator<T,P,Q>
{
public:
	//! Virtual destructor to help with cleanup;
	virtual ~SymEvenOddPrecLinearOperator() {}

	//! Only defined on the odd lattice
	const Subset& subset() const {return rb[1];}

	//! Return the fermion BC object for this linear operator
	virtual const FermBC<T,P,Q>& getFermBC() const = 0;

	//! Apply the even-even block onto a source std::vector
	/*! This does not need to be optimized */
	virtual void scaleEvenEvenLinOp(T& chi, const T& psi,
			enum PlusMinus isign) const = 0;

	//! Apply the inverse of the even-even block onto a source std::vector
	virtual void scaleEvenEvenInvLinOp(T& chi, const T& psi,
			enum PlusMinus isign) const = 0;

	//! Apply the the odd-odd block onto a source std::vector
	virtual void scaleOddOddLinOp(T& chi, const T& psi,
			enum PlusMinus isign) const = 0;

	//! Apply the inverse of the odd-odd block onto a source std::vector
	virtual void scaleOddOddInvLinOp(T& chi, const T& psi,
			enum PlusMinus isign) const = 0;

	//! Apply the the even-odd block onto a source std::vector
	virtual void unprecEvenOddLinOp(T& chi, const T& psi,
			enum PlusMinus isign) const = 0;

	//! Apply the the odd-even block onto a source std::vector
	virtual void unprecOddEvenLinOp(T& chi, const T& psi,
			enum PlusMinus isign) const = 0;


	//! Apply the EvenOdd Linop for which we have a category default)
	virtual void evenOddLinOp(T& chi, const T& psi,
			enum PlusMinus isign) const {

		T tmp; moveToFastMemoryHint(tmp);

		if( isign == PLUS ) {
			unprecEvenOddLinOp(tmp, psi, PLUS);
			scaleEvenEvenInvLinOp(chi,tmp,PLUS);
		}
		else {
			scaleOddOddInvLinOp(tmp, psi, MINUS);
			unprecEvenOddLinOp(chi, tmp, MINUS);
		}

		// Do I want to apply BC's here?
		getFermBC().modifyF(chi);
	}

	virtual void oddEvenLinOp(T& chi, const T& psi,
			enum PlusMinus isign) const {

		T tmp; moveToFastMemoryHint(tmp);

		if( isign == PLUS ) {
			unprecOddEvenLinOp(tmp, psi, PLUS);
			scaleOddOddInvLinOp(chi,tmp,PLUS);
		}
		else {
			scaleEvenEvenInvLinOp(tmp, psi, MINUS);
			unprecOddEvenLinOp(chi, tmp, MINUS);
		}

		// Do I want to apply BCs here?
		getFermBC().modifyF(chi);
	}

	// The odd odd block. Typically this is the identity
	// but may include something like a twisted clover term.
	virtual void oddOddLinOp(T& chi, const T& psi, enum PlusMinus isign) const =0;

	//! Apply the operator onto a source std::vector
	virtual void operator() (T& chi, const T& psi,
			enum PlusMinus isign) const
	{
		T   tmp1; moveToFastMemoryHint(tmp1);
		T   tmp2; moveToFastMemoryHint(tmp2);

		/*  t1 =  M_eo psi */
		evenOddLinOp(tmp1,psi,isign);

		/*  t2 = M_oe t1 */
		oddEvenLinOp(tmp2,tmp1,isign);

		/* chi = M_oo psi - M_oe M_eo psi = (M_oo - M_oe M_eo) psi */
		/* NB: the construction takes care of daggering
		 *       as the M_eo^\dagger will reverse operator order etc */
		oddOddLinOp(chi,psi,isign);
		chi[rb[1]] -=  tmp2;
		getFermBC().modifyF(chi);
	}

	//! Apply the Jacobi Operator
	virtual void jacobiOp(T& chi, const T& psi,
			enum PlusMinus isign) const
	{
		T tmp; moveToFastMemoryHint(tmp);

		//         [   1     A^{-1}_ee D_eo ]
		//         [                        ]
		//         [ A^{-1}_oo D_oe    M_oo    ]
		//
		//   Or hermitian conjugate
		oddOddLinOp(chi,psi, isign);
		chi[rb[0]] = psi;

		evenOddLinOp(tmp, psi, isign);
		chi[ rb[0]] += tmp;

		oddEvenLinOp(tmp, psi, isign);
		chi[ rb[1]] += tmp;

		getFermBC().modifyF(chi);
	}

	//! Apply the UNPRECONDITIONED operator onto a source std::vector
	/*! Mainly intended for debugging */
	virtual void unprecLinOp(T& chi, const T& psi,
			enum PlusMinus isign) const
	{
		T   tmp1;  moveToFastMemoryHint(tmp1);
		QDPIO::cout << "M_Unprec : " << std::endl;
		if ( isign == PLUS ) {
			// Apply Jacobi Operator and scale for M_unprec:
			//
			//   [ A_ee    0   ]  [   1     A^{-1}_ee D_eo    ]
			//   [             ]  [                           ]
			//   [   0   A_oo  ]  [ A^{-1}_oo D_oe    M_oo    ]
			//
			jacobiOp(tmp1, psi, isign);

			scaleEvenEvenLinOp(chi,tmp1,isign);
			scaleOddOddLinOp(chi,tmp1, isign);
		}
		else {
			// For Hermitian conjugate scale first then apply Jacobi operator
			//
			// [ 1     (D^\dag)_eo A^{-dag}_oo ] [ A^dag_ee     0  ]
			// [                               ] [                 ]
			// [ (D^\dag)_oe A^{-dag}_ee     1 ] [  0     A^dag_oo ]
			scaleEvenEvenLinOp(tmp1,psi,isign);
			scaleOddOddLinOp(tmp1,psi,isign);
			jacobiOp(chi,tmp1,isign);

		}
		getFermBC().modifyF(chi);
	}


	//! Compute force coming from EvenOdd (A^{-1}_ee D_eo) term
	virtual void derivEvenOddLinOp(P& ds_u, const T& chi, const T& psi,
			enum PlusMinus isign) const = 0;

	//! Compute force coming from OddEven (A^{-1}_oo D_oe ) term
	virtual void derivOddEvenLinOp(P& ds_u, const T& chi, const T& psi,
			enum PlusMinus isign) const = 0;


	//! Compute the force coming from the M_oo part. Typically this would be zero
    virtual void derivOddOddLinOp(P& ds_u, const T& chi, const T& psi,
    		enum PlusMinus isign) const = 0;

	//! Deriv
	virtual void deriv(P& ds_u, const T& chi, const T& psi,
			enum PlusMinus isign) const
	{
		T M_eo_psi; moveToFastMemoryHint(M_eo_psi);
		T M_oe_dag_chi; moveToFastMemoryHint(M_oe_dag_chi);

		enum PlusMinus msign = ( isign == PLUS ) ? MINUS : PLUS;

		evenOddLinOp( M_eo_psi, psi, isign);
		evenOddLinOp( M_oe_dag_chi, chi, msign);

		ds_u.resize(Nd);
		ds_u = zero;
		P ds_tmp;
		ds_tmp.resize(Nd);

		derivOddOddLinOp(ds_u, chi,psi,isign);

		ds_tmp = zero;
		derivOddEvenLinOp(ds_tmp,chi,M_eo_psi,isign);
		ds_u -= ds_tmp;

		ds_tmp = zero;
		derivEvenOddLinOp(ds_tmp, M_oe_dag_chi,psi,isign);
		ds_u -= ds_tmp;
	

		getFermBC().zero(ds_u);
	}

	//! Apply the derivative of the operator onto a source std::vector
	/*! User should make sure deriv routines do a resize  */
	virtual void derivMultipole(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi,
			enum PlusMinus isign) const
	{
		T M_eo_psi; moveToFastMemoryHint(M_eo_psi);
		T M_oe_dag_chi; moveToFastMemoryHint(M_oe_dag_chi);

		enum PlusMinus msign = ( isign == PLUS ) ? MINUS : PLUS;


		ds_u.resize(Nd);
		ds_u = zero;
		P ds_tmp;
		ds_tmp.resize(Nd);
		for(int i=0; i < chi.size(); ++i) {
			ds_tmp = zero;
			derivOddOddLinOp(ds_tmp, chi[i],psi[i], isign);
			ds_u += ds_tmp;

			evenOddLinOp( M_eo_psi, psi[i], isign);
			evenOddLinOp( M_oe_dag_chi, chi[i], msign);

			derivOddEvenLinOp(ds_tmp,chi[i],M_eo_psi,isign);
			ds_u -= ds_tmp;
			derivEvenOddLinOp(ds_tmp, M_oe_dag_chi,psi[i],isign);
			ds_u -= ds_tmp;
		}
		getFermBC().zero(ds_u);

	}

};


} // End namespace Chroma

#endif

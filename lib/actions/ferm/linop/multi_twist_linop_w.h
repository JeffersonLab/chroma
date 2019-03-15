/*
 * multi_twist_linop_w.h
 *
 *  Created on: Mar 14, 2019
 *      Author: bjoo
 */

#ifndef LIB_ACTIONS_FERM_LINOP_MULTI_TWIST_LINOP_W_H_
#define LIB_ACTIONS_FERM_LINOP_MULTI_TWIST_LINOP_W_H_

#include "chromabase.h"
#include "actions/ferm/linop/shifted_linop_w.h"

namespace Chroma  {

template<typename T, typename P, typename Q, template<typename, typename,
		typename > class LinOp>
class MultiTwistLinOp: public DiffLinearOperatorArray<T, P, Q> {
public:
	MultiTwistLinOp(const Handle< LinOp<T,P,Q> >& base_op, const multi1d<Real>& Twists) :
		_theLinOp(base_op), _Twists(Twists), _N(Twists.size()){}

	~MultiTwistLinOp() {}

	void operator()(multi1d<T>& psi, const multi1d<T>& chi, enum PlusMinus isign) const {
		AssertCompatible(chi,psi);
		for(int i=0; i < _N; ++i) {
			TwistedShiftedLinOp<T,P,Q,LinOp> theShiftedOp((*_theLinOp), _Twists[i]);
			theShiftedOp(psi[i],chi[i],isign);
		}
	}

	//! Apply the derivative of the operator onto a source std::vector
	/*! Default implementation */
	void deriv(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi,
			enum PlusMinus isign) const override
	{
		assertArraySizes(chi,psi,_N);
		ds_u.resize(Nd);
		ds_u = zero;

		P ds_tmp;
		for(int i=0; i < _N; ++i) {
			TwistedShiftedLinOp<T,P,Q,LinOp> theShiftedOp((*_theLinOp), _Twists[i]);
			theShiftedOp.deriv(ds_tmp, chi[i], psi[i], isign);
			ds_u += ds_tmp;
		}
	}


	int size() const {
		return _N;
	}

	//! Return the subset on which the operator acts
	const Subset& subset() const override {
		return _theLinOp->subset();
	}

	//! Return the fermion BC object for this linear operator
	const FermBC<T,P,Q>& getFermBC() const override {
		return _theLinOp->getFermBC();
	}

	const multi1d<Real>& getTwists() const {
		return _Twists;
	}
private:
	inline
	void AssertCompatible(const multi1d<T>& psi, const multi1d<T>& chi) const {
		if (psi.size() != chi.size() ) {
			QDPIO::cout << "MultiTwistLinOp: chi.size() != psi.size() " << std::endl;
			QDP_abort(1);
		}
		if( psi.size() != _N ) {
			QDPIO::cout << "MultiTwistLinOp: psi.size() != N " << std::endl;
			QDP_abort(1);
		}
	}
	const Handle< LinOp<T,P,Q> > _theLinOp;
	const multi1d<Real> _Twists;
	const int _N;
};

template<typename T, typename P, typename Q>
using SymEvenOddPrecLogDetMultiTwistLinOp = MultiTwistLinOp<T,P,Q,SymEvenOddPrecLogDetLinearOperator>;

template<typename T, typename P, typename Q>
using EvenOddPrecMultiTwistLinOp = MultiTwistLinOp<T,P,Q,EvenOddPrecLinearOperator>;


} // Namespace

#endif /* LIB_ACTIONS_FERM_LINOP_MULTI_TWIST_LINOP_W_H_ */

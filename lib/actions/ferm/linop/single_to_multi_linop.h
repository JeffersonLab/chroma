/*
 * single_to_multi_linop.h
 *
 *  Created on: Mar 11, 2019
 *      Author: bjoo
 */

#ifndef LIB_ACTIONS_FERM_LINOP_SINGLE_TO_MULTI_LINOP_H_
#define LIB_ACTIONS_FERM_LINOP_SINGLE_TO_MULTI_LINOP_H_

#include "chromabase.h"
#include "linearop.h"
#include <iostream>
namespace Chroma {

namespace {
	template<typename T>
	void assertArraySizes(const multi1d<T>& chi, const multi1d<T>& psi, const int N)
	{
		if ( chi.size() != psi.size() ) {
			QDPIO::cout << "Array sizes dont match in SingleRHSToMultiRHSDiffLinOp" << std::endl;
			QDP_abort(1);
	    }
		if ( chi.size() != N ) {
			QDPIO::cout << "Array sizes do not match N" << std::endl;
			QDP_abort(1);
		}
	}
}

template<typename T>
class SingleRHSToMultiRHSProxyLinOp : public LinearOperatorArray<T>
{
public:
	SingleRHSToMultiRHSProxyLinOp(const Handle< LinearOperator<T> >& base_op, const int N)
		: _theLinOp(base_op), _N(N) {}
	~SingleRHSToMultiRHSProxyLinOp() {}

	int size() const {
		return _N;
	}
	 //! Return the subset on which the operator acts
	const Subset& subset() const override {
		return _theLinOp->subset();
	}

	//! Apply the operator onto a source std::vector to some precision
	void operator() (multi1d<T>& chi, const multi1d<T>& psi,
				     enum PlusMinus isign) const override
	{
		assertArraySizes(chi,psi,_N);
		for(int i=0; i < _N; ++i) {
			(*_theLinOp)(chi[i],psi[i],isign);
		}
	}

private:

	const Handle< LinearOperator<T> >& _theLinOp;
	const int _N;
};


template<typename T, typename P, typename Q>
class SingleRHSToMultiRHSProxyDiffLinOp : public DiffLinearOperatorArray<T,P,Q>
{
public:
	SingleRHSToMultiRHSProxyDiffLinOp(const Handle< DiffLinearOperator<T,P,Q> >& base_op,
				const int N) : _theLinOp(base_op), _N(N) {}
	~SingleRHSToMultiRHSProxyDiffLinOp() {}

	//! Return the subset on which the operator acts
	const Subset& subset() const override {
		return _theLinOp->subset();
	}

	int size() const {
		return _N;
	}

	//! Return the fermion BC object for this linear operator
	const FermBC<T,P,Q>& getFermBC() const override {
		return _theLinOp->getFermBC();
	}

	//! Apply the operator onto a source std::vector to some precision
	void operator() (multi1d<T>& chi, const multi1d<T>& psi,
				     enum PlusMinus isign) const override
	{
		assertArraySizes(chi,psi,_N);
		for(int i=0; i < _N; ++i) {
			(*_theLinOp)(chi[i],psi[i],isign);
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
			_theLinOp->deriv(ds_tmp, chi[i], psi[i], isign);
			ds_u += ds_tmp;
		}
	}


private:
	const Handle< DiffLinearOperator<T,P,Q> > _theLinOp;
	const int _N;
};



}
#endif /* LIB_ACTIONS_FERM_LINOP_SINGLE_TO_MULTI_LINOP_H_ */

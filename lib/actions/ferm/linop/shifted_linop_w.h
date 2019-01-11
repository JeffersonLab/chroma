/* even-odd preconditioned fermion linear operator with
 * shifted mass term i*\mu*\gamma5*A_oo
 * based on template, which can be used both to 
 * symmetric and asymmetric eo precondition
 */

#ifndef __SHIFTED_LINOP_W_H__
#define __SHIFTED_LINOP_W_H__

#include "seoprec_constdet_linop.h"
#include "eoprec_constdet_linop.h"
#include "seoprec_logdet_linop.h"
#include "eoprec_logdet_linop.h"
#include "eoprec_linop.h"
#include "actions/ferm/linop/dslash_w.h"

namespace Chroma
{
	using namespace QDP::Hints;

	template<typename T, typename P, typename Q,
		template<typename, typename, typename> class LinOp>
			class TwistedShiftedLinOp : public DiffLinearOperator<T,P,Q> {};

	// Symmetric partial specialization
	template<typename T, typename P, typename Q>
		class TwistedShiftedLinOp<T, P, Q, SymEvenOddPrecLogDetLinearOperator> : public DiffLinearOperator<T,P,Q>  {
			public:
				TwistedShiftedLinOp(const SymEvenOddPrecLogDetLinearOperator<T,P,Q>& base_op_, const Real& mu_):
					base_op(base_op_), mu(mu_){}

				const Subset& subset() const override { return base_op.subset(); }

				const FermBC<T,P,Q>& getFermBC() const override  { return base_op.getFermBC() ; }

				void operator()(T& out, const T& in, enum PlusMinus isign) const override
				{
					START_CODE();
					T tmp;
					moveToFastMemoryHint(tmp);

					if(isign == PLUS){
						// shift with i\gamma_5 A_oo
						base_op(out, in, isign);
						base_op.unprecOddOddLinOp(tmp, in, isign);
						out[rb[1]] += mu*(Gamma(15)*timesI(tmp));
					}else{
						// shift with -i\gamma_5 A_oo
						base_op(out, in, isign);
						base_op.unprecOddOddLinOp(tmp, in, isign);
						out[rb[1]] -= mu*(Gamma(15)*timesI(tmp));
					}
					END_CODE();
				}

				void deriv(P& ds_u, const T& Y, const T& X, enum PlusMinus isign)const override
				{
					START_CODE();
					//const SymEvenOddPrecConstDetLinearOperator<T,P,Q>& constdet_op = 
					//	static_cast<const SymEvenOddPrecConstDetLinearOperator<T,P,Q>&>(base_op);
					
					P ds_extra;
					ds_extra.resize(Nd);
					ds_u.resize(Nd);
					for(int i=0; i<Nd; ++i){
						ds_extra[i] = zero;
						ds_u[i] = zero;
					}
					//constdet_op.deriv(ds_u, Y, X, isign);
					base_op.deriv(ds_u, Y, X, isign);
					// add deriv of shifted mass term
					T Y_prime = zero;
					Y_prime = Gamma(15)*Y;
					base_op.derivUnprecOddOddLinOp(ds_extra, Y_prime, X, isign);
					for(int i=0; i<Nd; ++i)
					ds_extra[i] *= mu;

					if(isign == PLUS){
						for(int i=0; i<Nd; ++i)
						ds_u[i] += timesI(ds_extra[i]);
					}else{
						for(int i=0; i<Nd; ++i)
						ds_u[i] -= timesI(ds_extra[i]);
					}
					END_CODE();
				}

			private:
				const Real mu;
				//SymEvenOddPrecConstDetLinearOperator<T,P,Q> base_op;
				const SymEvenOddPrecLogDetLinearOperator<T,P,Q>& base_op;
		};

	// Asymmetric partial specialization
	// EvenOddPrecLogDetLinearOperator and EvenOddPrecConstDetLinearOperator both
	// inherit from EvenOddPrecLinearOperator
	// different from SymEvenOddPrecLogDetLinearOperator inherit from
	// SymEvenOddPrecConstDetLinearOperator
	// therefore slightly different specialization parameter for LinOp
	template<typename T, typename P, typename Q>
		class TwistedShiftedLinOp<T, P, Q, EvenOddPrecLinearOperator> :
			public DiffLinearOperator<T,P,Q>  {
				public:
				TwistedShiftedLinOp(const EvenOddPrecLinearOperator<T,P,Q>& base_op_, const Real& mu_):
					base_op(base_op_), mu(mu_){}

				const Subset& subset() const override { return base_op.subset(); }
				const FermBC<T,P,Q>& getFermBC() const override  { return base_op.getFermBC(); }
				void operator()(T& out, const T& in, enum PlusMinus isign)const override
				{
					START_CODE();
					if(isign == PLUS){
						// shift with i\gamma_5
						base_op(out, in, isign);
						out[rb[1]] += mu*(Gamma(15)*timesI(in));
					}else{
						// shift with -i\gamma_5
						base_op(out, in, isign);
						out[rb[1]] -= mu*(Gamma(15)*timesI(in));
					}
					END_CODE();
				}

				void deriv(P& ds_u, const T& Y, const T& X, enum PlusMinus isign)const override
				{
					// its +/- i gamma_5 so indep of gauge fields
					base_op.deriv(ds_u, Y, X, isign);
				}

			private:
				const Real mu;
				const EvenOddPrecLinearOperator<T,P,Q>& base_op;
		};
}
#endif

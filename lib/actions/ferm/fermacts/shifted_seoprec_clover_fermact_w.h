// -*- C++ -*-
/*! \file
 *  \brief Symmetric even-odd preconditioned Clover fermion action
 *  with shifted mass term
 */

#ifndef __SHIFTED_SEOPREC_CLOVER_FERMACT_W_H__
#define __SHIFTED_SEOPREC_CLOVER_FERMACT_W_H__

#include "seoprec_logdet_wilstype_fermact_w.h"
#include "actions/ferm/linop/lgherm_w.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"

namespace Chroma
{
	namespace ShiftSymEvenOddPrecCloverFermActEnv
	{
		extern const std::string name;
		bool registerAll();
	}

	//! Symmetric even-odd preconditioned Clover fermion action
	/*! \ingroup fermacts
	 *
	 * Symmetric Even-odd preconditioned clover fermion action. 
	 * Only defined on odd subset.
	 */

	class ShiftSymEvenOddPrecCloverFermAct: public SymEvenOddPrecLogDetWilsonTypeFermAct<LatticeFermion,
	multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
	{

		public:
			// Typedefs to save typing
			typedef LatticeFermion               T;
			typedef multi1d<LatticeColorMatrix>  P;
			typedef multi1d<LatticeColorMatrix>  Q;

			//! Partial constructor
			ShiftSymEvenOddPrecCloverFermAct() {}

			//! General FermState
		//	ShiftSymEvenOddPrecCloverFermAct(Handle< CreateFermState<T,P,Q> > cfs_, 
		//			const CloverFermActParams& param_, const Real& mu_) : 
		//		cfs(cfs_), param(param_), mu(mu_) {}

			ShiftSymEvenOddPrecCloverFermAct(Handle< CreateFermState<T,P,Q> > cfs_, 
					const CloverFermActParams& param_) : 
				cfs(cfs_), param(param_){}
			
			//! Copy constructor
			ShiftSymEvenOddPrecCloverFermAct(const ShiftSymEvenOddPrecCloverFermAct& a) : 
				cfs(a.cfs), param(a.param), mu(a.mu) {}

			//! Produce a linear operator for this action
			SymEvenOddPrecLogDetLinearOperator<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const;

			//! Produce the gamma_5 hermitian operator H_w
			LinearOperator<LatticeFermion>* hermitianLinOp(Handle< FermState<T,P,Q> > state) const 
			{ 
				return new lgherm<LatticeFermion>(linOp(state));
			}

			//! Set shifted mass mu
			void setMu(const Real& mu_){
				mu = mu_;
			}
			//! Destructor is automatic
			~ShiftSymEvenOddPrecCloverFermAct() {}

		protected:
			//! Return the fermion BC object for this action
			const CreateFermState<T,P,Q>& getCreateState() const {return *cfs;}

			//! Assignment
			void operator=(const ShiftSymEvenOddPrecCloverFermAct& a) {}

		private:
			Handle< CreateFermState<T,P,Q> >  cfs;
			CloverFermActParams param;
			Real mu;
	};

} // End Namespace Chroma
#endif

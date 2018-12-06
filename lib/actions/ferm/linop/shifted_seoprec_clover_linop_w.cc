// -*- C++ -*-
/*! \file
 *  \brief Symmetric even-odd preconditioned Clover fermion linear operator
 *  \with shifted mass term i*\mu*\gamma5*A_oo
 */
#include "actions/ferm/linop/shifted_seoprec_clover_linop_w.h"

namespace Chroma
{
	using namespace QDP::Hints;

	//! Creation routine with Anisotropy
	/*!
	 * \param u_ 	    gauge field     	       (Read)
	 * \param param_  fermion kappa   	       (Read)
	 */
	void ShiftSymEvenOddPrecCloverLinOp::create(Handle<FermState<T,P,Q> >fs,
			const CloverFermActParams& param_)
	{
		START_CODE();
		param = param_;
		// TODO:: may remove twisted mass term in this shifted mass linop
		QDPIO::cout << "Using Twisted Mass: " << param.twisted_m_usedP << std::endl;
		if(param.twisted_m_usedP){
			QDPIO::cout << "Twisted Mass is " << param.twisted_m << std::endl;
		}
		clov.create(fs, param);

		invclov.create(fs, param, clov);  // make a copy
		invclov.choles(0);  // invert the cb=0 part
		invclov.choles(1);  // invert the cb=1 part

		D.create(fs, param.anisoParam);

		clov_deriv_time = 0;
		clov_apply_time = 0;
		END_CODE();
	}

	//! Apply the odd-odd block onto a source std::vector
	void ShiftSymEvenOddPrecCloverLinOp::
		unprecOddOddLinOp(LatticeFermion& chi, const LatticeFermion& psi,
				enum PlusMinus isign) const
		{
			START_CODE();
			clov.apply(chi, psi, isign, 1);
			END_CODE();
		}

	//! Apply the inverse of the odd-odd block onto a source std::vector
	void ShiftSymEvenOddPrecCloverLinOp::
		unprecOddOddInvLinOp(LatticeFermion& chi, const LatticeFermion& psi,
				enum PlusMinus isign) const
		{
			START_CODE();
			invclov.apply(chi, psi, isign, 1);
			END_CODE();
		}

	//! Apply the even-even block onto a source std::vector
	void ShiftSymEvenOddPrecCloverLinOp::
		unprecEvenEvenLinOp(LatticeFermion& chi, const LatticeFermion& psi,
				enum PlusMinus isign) const
		{
			START_CODE();
			clov.apply(chi, psi, isign, 0);
			END_CODE();
		}

	//! Apply the inverse of the even-even block onto a source std::vector
	void ShiftSymEvenOddPrecCloverLinOp::
		unprecEvenEvenInvLinOp(LatticeFermion& chi, const LatticeFermion& psi,
				enum PlusMinus isign) const
		{
			START_CODE();
			invclov.apply(chi, psi, isign, 0);
			END_CODE();
		}


	//! Apply even-odd linop component
	/*!
	 * The operator acts on the entire even sublattice
	 *
	 * \param chi 	  Pseudofermion field     	       (Write)
	 * \param psi 	  Pseudofermion field     	       (Read)
	 * \param isign   Flag ( PLUS | MINUS )   	       (Read)
	 */
	void ShiftSymEvenOddPrecCloverLinOp::
		unprecEvenOddLinOp(LatticeFermion& chi, const LatticeFermion& psi,
				enum PlusMinus isign) const
		{
			START_CODE();

			Real mhalf = -0.5;
			D.apply(chi, psi, isign, 0);
			chi[rb[0]] *= mhalf;

			END_CODE();
		}

	//! Apply odd-even linop component
	/*!
	 * The operator acts on the entire odd sublattice
	 *
	 * \param chi 	  Pseudofermion field     	       (Write)
	 * \param psi 	  Pseudofermion field     	       (Read)
	 * \param isign   Flag ( PLUS | MINUS )   	       (Read)
	 */
	void ShiftSymEvenOddPrecCloverLinOp::
		unprecOddEvenLinOp(LatticeFermion& chi, const LatticeFermion& psi,
				enum PlusMinus isign) const
		{
			START_CODE();

			Real mhalf = -0.5;
			D.apply(chi, psi, isign, 1);
			chi[rb[1]] *= mhalf;

			END_CODE();
		}


	//! Apply even-odd preconditioned Clover fermion linear operator
	/*!
	 * \param chi 	  Pseudofermion field     	       (Write)
	 * \param psi 	  Pseudofermion field     	       (Read)
	 * \param isign   Flag ( PLUS | MINUS )   	       (Read)
	 */
	void ShiftSymEvenOddPrecCloverLinOp::
		operator()(T& chi, const T& psi,
				enum PlusMinus isign) const
		{
			START_CODE();

			T tmp1; moveToFastMemoryHint(tmp1);
			T tmp2; moveToFastMemoryHint(tmp2);
			Real mquarter = -0.25;

			if(isign == PLUS){
				// tmp2_o = A^(-1)_oo D_oe A^(-1)_ee D_eo psi_o
				D.apply(tmp1, psi, isign, 0);
				invclov.apply(tmp2, tmp1, isign, 0);
				D.apply(tmp1, tmp2, isign, 1);
				invclov.apply(tmp2, tmp1, isign, 1);

				// tmp1_o for shited mass term
				clov.apply(tmp1, psi, isign, 1);

				// add shifted mass term i*mu*gamma_5*A_oo
				chi[rb[1]] = psi + mquarter*tmp2 + Gamma(15)*timesI(tmp1);
			}else{
				invclov.apply(tmp1, psi, isign, 1);
				D.apply(tmp2, tmp1, isign, 0);
				invclov.apply(tmp1, tmp2, isign, 0);
				D.apply(tmp2, tmp1, isign, 1);

				// tmp1_o for shifted mass term
				clov.apply(tmp1, psi, isign, 1);

				// add shifted mass term -i*mu*gamma_5*A_oo
				chi[rb[1]] = psi + mquarter*tmp2 - Gamma(15)*timesI(tmp1);
			}
			getFermBC().modifyF(chi);

			END_CODE();
		}

	//! Deriv of even-even block
	void ShiftSymEvenOddPrecCloverLinOp::
		derivUnprecEvenEvenLinOp(multi1d<LatticeColorMatrix>& ds_u,
				const LatticeFermion& chi, const LatticeFermion& psi,
				enum PlusMinus isign) const
		{
			START_CODE();
			clov.deriv(ds_u, chi, psi, isign, 0);
			END_CODE();
		}

	//! Deriv of odd-odd block
	void ShiftSymEvenOddPrecCloverLinOp::
		derivUnprecOddOddLinOp(multi1d<LatticeColorMatrix>& ds_u,
				const LatticeFermion& chi, const LatticeFermion& psi,
				enum PlusMinus isign) const
		{
			START_CODE();
			clov.deriv(ds_u, chi, psi, isign, 1);
			END_CODE();
		}

	//! Deriv of Trace Log of even-even block
	void ShiftSymEvenOddPrecCloverLinOp::
		derivLogDetEvenEvenLinOp(multi1d<LatticeColorMatrix>& ds_u,
				enum PlusMinus isign) const
		{
			START_CODE();
			invclov.derivTrLn(ds_u, isign, 0);
			END_CODE();
		}

	//! Deriv of Trace Log of odd-odd block
	void ShiftSymEvenOddPrecCloverLinOp::
		derivLogDetOddOddLinOp(multi1d<LatticeColorMatrix>& ds_u,
				enum PlusMinus isign) const
		{
			START_CODE();
			invclov.derivTrLn(ds_u, isign, 1);
			END_CODE();
		}

	//! Deriv of even-odd block
	void ShiftSymEvenOddPrecCloverLinOp::
		derivUnprecEvenOddLinOp(multi1d<LatticeColorMatrix>& ds_u,
				const LatticeFermion& chi, const LatticeFermion& psi,
				enum PlusMinus isign) const
		{
			START_CODE();
			ds_u.resize(Nd);
			D.deriv(ds_u, chi, psi, isign, 0);
			for(int id=0; id<Nd; id++){
				ds_u[id] *= Real(-0.5);
			}
			END_CODE();
		}

	//! Deriv of odd-even block
	void ShiftSymEvenOddPrecCloverLinOp::
		derivUnprecOddEvenLinOp(multi1d<LatticeColorMatrix>& ds_u,
				const LatticeFermion& chi, const LatticeFermion& psi,
				enum PlusMinus isign) const
		{
			START_CODE();
			ds_u.resize(Nd);
			D.deriv(ds_u, chi, psi, isign, 1);
			for(int id=0; id<Nd; id++){
				ds_u[id] *= Real(-0.5);
			}
			END_CODE();
		}

	//! Deriv
	void ShiftSymEvenOddPrecCloverLinOp::
		deriv(P& ds_u, const T& chi, const T& psi,
				enum PlusMinus isign) const
		{
			T M_eo_psi; moveToFastMemoryHint(M_eo_psi);
			T M_oe_dag_chi; moveToFastMemoryHint(M_oe_dag_chi);
			T tmp;
			T Ainv_left;
			T Ainv_right;

			if(isign == PLUS){
				D.apply(tmp, psi, PLUS, 0);
				invclov.apply(M_eo_psi, tmp, PLUS, 0);

				invclov.apply(tmp, chi, MINUS, 1);
				D.apply(M_oe_dag_chi, tmp, MINUS, 0);


				M_eo_psi[rb[0]] *= Real(-0.5);
				M_oe_dag_chi[rb[0]] *= Real(-0.5);

				ds_u.resize(Nd);
				ds_u = zero;
				P ds_tmp;
				ds_tmp.resize(Nd);

				// total deriv (chain rule)
				D.apply(tmp, M_eo_psi, PLUS, 1);
				invclov.apply(Ainv_right, tmp, PLUS, 1);
				invclov.apply(Ainv_left, chi, MINUS, 1);
				clov.deriv(ds_u, Ainv_left, Ainv_right, PLUS, 1);

				invclov.apply(tmp, chi, MINUS, 1);
				D.deriv(ds_tmp, tmp, M_eo_psi, PLUS, 1);
				ds_u -= ds_tmp;

				invclov.apply(Ainv_left, M_oe_dag_chi, PLUS, 0);
				clov.deriv(ds_tmp, Ainv_left, M_eo_psi, PLUS, 0);
				// match the total coeffient
				for(int id=0; id<Nd; ++id){
					ds_tmp[id] *= Real(-2.0);
				}
				ds_u += ds_tmp;

				D.deriv(ds_tmp, Ainv_left, psi, PLUS, 0);
				ds_u -= ds_tmp;

				// shifted mass term
				clov.deriv(ds_tmp, chi, Gamma(15)*timesI(psi), PLUS, 1);
				// match the total coeffient
				for(int id=0; id<Nd; ++id){
					ds_tmp[id] *= Real(-2.0);
					//ds_u[id] += mu* Gamma(15) * ds_tmp[id];
					ds_u[id] += mu * ds_tmp[id];
				}
			}else{
				invclov.apply(tmp, psi, MINUS, 1);
				D.apply(M_eo_psi, tmp, MINUS, 0);

				D.apply(tmp, chi, PLUS, 0);
				invclov.apply(M_oe_dag_chi, tmp, PLUS, 0);

				M_eo_psi[rb[0]] *= Real(-0.5);
				M_oe_dag_chi[rb[0]] *= Real(-0.5);

				ds_u.resize(Nd);
				ds_u = zero;
				P ds_tmp;
				ds_tmp.resize(Nd);

				// total deriv (chain rule)
				invclov.apply(Ainv_right, M_eo_psi, MINUS, 0);
				clov.deriv(ds_u, M_oe_dag_chi, Ainv_right, MINUS, 0);
				// match the total coeffient
				for(int id=0; id<Nd; ++id){
					ds_u[id] *= Real(-2.0);
				}

				invclov.apply(tmp, M_eo_psi, MINUS, 0);
				D.deriv(ds_tmp, chi, tmp, MINUS, 1);
				ds_u -= ds_tmp;

				invclov.apply(tmp, psi, MINUS, 1);
				D.deriv(ds_tmp, M_oe_dag_chi, tmp, MINUS, 0);
				ds_u -= ds_tmp;

				D.apply(tmp, M_oe_dag_chi, PLUS, 1);
				invclov.apply(Ainv_left, tmp, PLUS, 1);
				invclov.apply(Ainv_right, psi, MINUS, 1);
				clov.deriv(ds_tmp, Ainv_left, Ainv_right, MINUS, 1);

				// shifted mass term
				clov.deriv(ds_tmp, chi, Gamma(15)*timesI(psi), MINUS, 1);
				// match the total coeffient
				for(int id=0; id<Nd; ++id){
					ds_tmp[id] *= Real(-2.0);
					ds_u[id] -= mu * ds_tmp[id];
				}
			}

			// make the total coeffient right
			for(int id=0; id<Nd; ++id){
				ds_u[id] *= Real(-0.5);
			}

			getFermBC().zero(ds_u);
		}

	void ShiftSymEvenOddPrecCloverLinOp::
		derivMultipole(P& ds_u, const multi1d<T>& chi, const multi1d<T>& psi,
				enum PlusMinus isign) const
		{
			if( chi.size() != psi.size() ) {
				QDPIO::cerr << "Incompatible array sizes in ShiftSymEvenOddPrecCloverOp::derivMultipole" << std::endl;
				QDP_abort(1);
			}
			int n_poles= psi.size();
			multi1d<T> M_eo_psi(n_poles);
			multi1d<T> M_oe_dag_chi(n_poles);
			T tmp;
			multi1d<T> Ainv_left(n_poles);
			multi1d<T> Ainv_right(n_poles);

			if(isign == PLUS){
				for(int i=0; i<n_poles; ++i){
					D.apply(tmp, psi[i], PLUS, 0);
					invclov.apply(M_eo_psi[i], tmp, PLUS, 0);

					invclov.apply(tmp, chi[i], MINUS, 1);
					D.apply(M_oe_dag_chi[i], tmp, MINUS, 0);

					(M_eo_psi[i])[rb[0]] *= Real(-0.5);
					(M_oe_dag_chi[i])[rb[0]] *= Real(-0.5);
				}
				ds_u.resize(Nd);
				ds_u = zero;
				P ds_tmp;
				ds_tmp.resize(Nd);

				// total deriv (chain rule)
				for(int i=0; i<n_poles; ++i){
					D.apply(tmp, M_eo_psi[i], PLUS, 1);
					invclov.apply(Ainv_right[i], tmp, PLUS, 1);
					invclov.apply(Ainv_left[i], chi[i], MINUS, 1);
				}
				clov.derivMultipole(ds_u, Ainv_left, Ainv_right, PLUS, 1);

				for(int i=0; i<n_poles; ++i){
					invclov.apply(tmp, chi[i], MINUS, 1);
					D.deriv(ds_tmp, tmp, M_eo_psi[i], PLUS, 1);
					ds_u -= ds_tmp;
				}

				for(int i=0; i<n_poles; ++i){
					invclov.apply(Ainv_left[i], M_oe_dag_chi[i], PLUS, 0);
				}
				clov.derivMultipole(ds_tmp, Ainv_left, M_eo_psi, PLUS, 0);
				// match the total coeffient
				for(int id=0; id<Nd; ++id){
					ds_tmp[id] *= Real(-2.0);
				}
				ds_u += ds_tmp;

				for(int i=0; i<n_poles; ++i){
					D.deriv(ds_tmp, Ainv_left[i], psi[i], PLUS, 0);
					ds_u -= ds_tmp;
				}
				// shifted mass term
				multi1d<T> psi_timesI;
				for(int i=0; i<n_poles; ++i){
					psi_timesI[i] = Gamma(15)*timesI(psi[i]);
				}
				clov.derivMultipole(ds_tmp, chi, psi_timesI, PLUS, 1);
				// match the total coeffient
				for(int id=0; id<Nd; ++id){
					ds_tmp[id] *= Real(-2.0);
					ds_u[id] += mu * ds_tmp[id];
				}
			}else{
				for(int i=0; i<n_poles; ++i){
					invclov.apply(tmp, psi[i], MINUS, 1);
					D.apply(M_eo_psi[i], tmp, MINUS, 0);

					D.apply(tmp, chi[i], PLUS, 0);
					invclov.apply(M_oe_dag_chi[i], tmp, PLUS, 0);

					(M_eo_psi[i])[rb[0]] *= Real(-0.5);
					(M_oe_dag_chi[i])[rb[0]] *= Real(-0.5);
				}
				ds_u.resize(Nd);
				ds_u = zero;
				P ds_tmp;
				ds_tmp.resize(Nd);

				// total deriv (chain rule)
				for(int i=0; i<n_poles; ++i){
					invclov.apply(Ainv_right[i], M_eo_psi[i], MINUS, 0);
				}
				clov.derivMultipole(ds_u, M_oe_dag_chi, Ainv_right, MINUS, 0);
				// match the total coeffient
				for(int id=0; id<Nd; ++id){
					ds_u[id] *= Real(-2.0);
				}

				for(int i=0; i<n_poles; ++i){
					invclov.apply(tmp, M_eo_psi[i], MINUS, 0);
					D.deriv(ds_tmp, chi[i], tmp, MINUS, 1);
					ds_u -= ds_tmp;

					invclov.apply(tmp, psi[i], MINUS, 1);
					D.deriv(ds_tmp, M_oe_dag_chi[i], tmp, MINUS, 0);
					ds_u -= ds_tmp;
				}

				for(int i=0; i<n_poles; ++i){
					D.apply(tmp, M_oe_dag_chi[i], PLUS, 1);
					invclov.apply(Ainv_left[i], tmp, PLUS, 1);
					invclov.apply(Ainv_right[i], psi[i], MINUS, 1);
				}
				clov.derivMultipole(ds_tmp, Ainv_left, Ainv_right, MINUS, 1);

				// shifted mass term
				multi1d<T> psi_timesI;
				for(int i=0; i<n_poles; ++i){
					psi_timesI[i] = Gamma(15)*timesI(psi[i]);
				}
				clov.derivMultipole(ds_tmp, chi, psi_timesI, MINUS, 1);
				// match the total coeffient
				for(int id=0; id<Nd; ++id){
					ds_tmp[id] *= Real(-2.0);
					ds_u[id] -= mu *  ds_tmp[id];
				}
			}
			for(int id=0; id<Nd; ++id){
				ds_u[id] *= Real(-0.5);
			}
			getFermBC().zero(ds_u);

		}

	//! Return flops performed by the operator()
	// TODO:: need modify
	unsigned long ShiftSymEvenOddPrecCloverLinOp::nFlops() const
	{
		unsigned long cbsite_flops = 2*D.nFlops()+2*clov.nFlops()+4*Nc*Ns;

		return cbsite_flops*(Layout::sitesOnNode()/2);
	}

	//! Get the log det of the even even part
	// BUt for now, return zero for testing.
	Double ShiftSymEvenOddPrecCloverLinOp::logDetEvenEvenLinOp(void) const  {
		return invclov.cholesDet(0);
	}

	//! Get the log det of the odd odd part
	// BUt for now, return zero for testing.
	Double ShiftSymEvenOddPrecCloverLinOp::logDetOddOddLinOp(void) const  {
		return invclov.cholesDet(1);
	}


}// End Namespace Chroma

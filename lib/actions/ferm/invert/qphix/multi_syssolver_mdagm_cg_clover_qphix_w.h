// -*- C++ -*-
/*! \file
 *  \brief Solve a M*psi=chi linear system by BiCGStab
 */

#ifndef __multi_syssolver_mdagm_clover_qphix_h__
#define __multi_syssolver_mdagm_clover_qphix_h__

#include "chroma_config.h"

#ifdef BUILD_QPHIX
#include "handle.h"
#include "state.h"
#include "syssolver.h"
#include "linearop.h"
#include "lmdagm.h"
#include "actions/ferm/invert/multi_syssolver_mdagm.h"
#include "actions/ferm/fermbcs/simple_fermbc.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/invert/qphix/multi_syssolver_qphix_clover_params.h"
#include "actions/ferm/linop/clover_term_w.h"
#include "actions/ferm/linop/eoprec_clover_linop_w.h"
#include "io/aniso_io.h"
#include <string>

#include "chroma_config.h"
#include "util/gauge/reunit.h"
#include "io/aniso_io.h"

// Header files from the dslash package
#include "qphix/qphix_config.h"
#include "qphix/geometry.h"
#include "qphix/dslash_utils.h"
#include "qphix/qdp_packer.h"
#include "qphix/clover.h"
#include "qphix/minvcg.h"
#include "actions/ferm/invert/qphix/qphix_vec_traits.h"
#include "qphix_singleton.h"
using namespace QDP;

namespace Chroma
{

//! Richardson system solver namespace
namespace MdagMMultiSysSolverQPhiXCloverEnv
{
//! Register the syssolver
bool registerAll();
}


//! Solve a Clover Fermion System using the QPhiX inverter
/*! \ingroup invert
 *** WARNING THIS SOLVER WORKS FOR Clover FERMIONS ONLY ***
 */
using namespace QPhiXVecTraits;
template<typename T, typename U>
class MdagMMultiSysSolverQPhiXClover : public MdagMMultiSystemSolver<T>
{
public:
	typedef multi1d<U> Q;
	typedef typename WordType<T>::Type_t REALT;
	typedef typename QPhiX::Geometry<REALT,VecTraits<REALT>::Vec,VecTraits<REALT>::Soa,VecTraits<REALT>::compress12>::FourSpinorBlock QPhiX_Spinor;
	typedef typename QPhiX::Geometry<REALT,VecTraits<REALT>::Vec,VecTraits<REALT>::Soa,VecTraits<REALT>::compress12>::SU3MatrixBlock QPhiX_Gauge;
	typedef typename QPhiX::Geometry<REALT,VecTraits<REALT>::Vec,VecTraits<REALT>::Soa,VecTraits<REALT>::compress12>::CloverBlock QPhiX_Clover;


	//! Constructor
	/*!
	 * \param M_        Linear operator ( Read )
	 * \param invParam  inverter parameters ( Read )
	 */
	MdagMMultiSysSolverQPhiXClover(Handle< LinearOperator<T> > A_,
			Handle< FermState<T,Q,Q> > state_,
			const MultiSysSolverQPhiXCloverParams& invParam_) :
				A(A_), invParam(invParam_),
				clov(new CloverTermT<T, U>()),
				invclov(new CloverTermT<T, U>())
	{
		QDPIO::cout << "MdagMMultiSysSolverQPhiXClover:" << std::endl;
		QDPIO::cout << "AntiPeriodicT is: " << invParam.AntiPeriodicT << std::endl;


		QDPIO::cout << "Veclen is " << VecTraits<REALT>::Vec << std::endl;
		QDPIO::cout << "Soalen is " << VecTraits<REALT>::Soa << std::endl;

		if ( VecTraits<REALT>::Soa > VecTraits<REALT>::Vec ) {
			QDPIO::cerr << "PROBLEM: Soalen > Veclen. Please set soalen appropriately (<=VECLEN) at compile time" << std::endl;
			QDP_abort(1);
		}

		Q u(Nd);
		for(int mu=0; mu < Nd; mu++) {
			u[mu] = state_->getLinks()[mu];
		}

		// Set up aniso coefficients
		multi1d<Real> aniso_coeffs(Nd);
		for(int mu=0; mu < Nd; mu++) aniso_coeffs[mu] = Real(1);

		bool anisotropy = invParam.CloverParams.anisoParam.anisoP;
		if( anisotropy ) {
			aniso_coeffs = makeFermCoeffs( invParam.CloverParams.anisoParam );
		}

		double t_boundary=(double)(1);
		// NB: In this case, the state will have boundaries applied.
		// So we only need to apply our own boundaries if Compression is enabled
		//
		if (invParam.AntiPeriodicT) {
			t_boundary=(double)(-1);
			// Flip off the boundaries -- Dslash expects them off..
			u[3] *= where(Layout::latticeCoordinate(3) == (Layout::lattSize()[3]-1),
					Real(t_boundary), Real(1));
		}

		const QPhiX::QPhiXCLIArgs& QPhiXParams = TheQPhiXParams::Instance();

		QDPIO::cout << "About to grap a Dslash" << std::endl;
		geom = new QPhiX::Geometry<REALT,
				VecTraits<REALT>::Vec,
				VecTraits<REALT>::Soa,
				VecTraits<REALT>::compress12>(Layout::subgridLattSize().slice(),
						QPhiXParams.getBy(),
						QPhiXParams.getBz(),
						QPhiXParams.getNCores(),
						QPhiXParams.getSy(),
						QPhiXParams.getSz(),
						QPhiXParams.getPxy(),
						QPhiXParams.getPxyz(),
						QPhiXParams.getMinCt());

		QDPIO::cout << " Allocating Clover" << std::endl << std::flush ;
		QPhiX_Clover* A_cb0=(QPhiX_Clover *)geom->allocCBClov();
		QPhiX_Clover* A_cb1=(QPhiX_Clover *)geom->allocCBClov();
		clov_packed[0] = A_cb0;
		clov_packed[1] = A_cb1;

		QDPIO::cout << " Allocating CloverInv " << std::endl << std::flush ;
		QPhiX_Clover* A_inv_cb0=(QPhiX_Clover *)geom->allocCBClov();
		QPhiX_Clover* A_inv_cb1=(QPhiX_Clover *)geom->allocCBClov();
		invclov_packed[0] = A_inv_cb0;
		invclov_packed[1] = A_inv_cb1;

		// Pack the gauge field
		QDPIO::cout << "Packing gauge field..."  << std::endl << std::flush ;
		QPhiX_Gauge* packed_gauge_cb0=(QPhiX_Gauge *)geom->allocCBGauge();
		QPhiX_Gauge* packed_gauge_cb1=(QPhiX_Gauge *)geom->allocCBGauge();

		QPhiX::qdp_pack_gauge<>(u, packed_gauge_cb0,packed_gauge_cb1, *geom);
		u_packed[0] = packed_gauge_cb0;
		u_packed[1] = packed_gauge_cb1;


		QDPIO::cout << "Creating Clover Term" << std::endl;
		CloverTerm clov_qdp;
		clov->create(state_, invParam.CloverParams);
		QDPIO::cout << "Inverting Clover Term" << std::endl;
		invclov->create(state_, invParam.CloverParams, (*clov));
		for(int cb=0; cb < 2; cb++) {
			invclov->choles(cb);
		}
		QDPIO::cout << "Done" << std::endl;
		QDPIO::cout << "Packing Clover term..." << std::endl;

		for(int cb=0; cb < 2; cb++) {
			QPhiX::qdp_pack_clover<>((*invclov), invclov_packed[cb], *geom, cb);
		}

		for(int cb=0; cb < 2; cb++) {
			QPhiX::qdp_pack_clover<>((*clov), clov_packed[cb], *geom, cb);
		}
		QDPIO::cout << "Done" << std::endl;

		QDPIO::cout << "Creating the Even Odd Operator" << std::endl;
		M=new QPhiX::EvenOddCloverOperator<REALT,
				VecTraits<REALT>::Vec,
				VecTraits<REALT>::Soa,
				VecTraits<REALT>::compress12>(u_packed,
						clov_packed[1],
						invclov_packed[0],
						geom,
						t_boundary,
						toDouble(aniso_coeffs[0]),
						toDouble(aniso_coeffs[3]));


		mcg_solver = new QPhiX::MInvCG<REALT,VecTraits<REALT>::Vec, VecTraits<REALT>::Soa, VecTraits<REALT>::compress12>((*M), invParam.MaxIter, invParam.MaxShifts);
		if( invParam.TuneP ) mcg_solver->tune();
	}


	//! Destructor is automatic
	~MdagMMultiSysSolverQPhiXClover()
	{

		// Need to unalloc all the memory...
		QDPIO::cout << "Destructing" << std::endl;

		geom->free(invclov_packed[0]);
		geom->free(invclov_packed[1]);
		invclov_packed[0] = 0x0;
		invclov_packed[1] = 0x0;

		geom->free(clov_packed[0]);
		geom->free(clov_packed[1]);
		clov_packed[0] = 0x0;
		clov_packed[1] = 0x0;

		geom->free(u_packed[0]);
		geom->free(u_packed[1]);
		u_packed[0] = 0x0;
		u_packed[1] = 0x0;

		delete geom;

		// Evertyhing else freed as handles freed.
	}

	//! Return the subset on which the operator acts
	const Subset& subset() const {return A->subset();}

	virtual SystemSolverResults_t operator()(multi1d<T>& psi,
			const multi1d<Real>& shifts,
			const T& chi) const
	{
		SystemSolverResults_t res;

		START_CODE();
		StopWatch swatch;
		swatch.reset(); swatch.start();

		int n_shift = shifts.size();
		QDPIO::cout << "operator(): n_shift = " << n_shift << std::endl;

		// Sanity check 1:
		if (n_shift > invParam.MaxShifts) {
			QDPIO::cerr << "n_shift="<<n_shift<<" is greater than MaxShifts="
					<< invParam.MaxShifts << std::endl;
			QDP_abort(1);
		}


		// Allocate and pack the RHS
		QPhiX_Spinor* chi_qphix;
		chi_qphix =(QPhiX_Spinor *)geom->allocCBFourSpinor();
		if (chi_qphix == 0x0) {
			QDPIO::cout << "Unable to allocate chi_qphix" << std::endl;
			QDP_abort(1);
		}
		QDPIO::cout << "Packing chi" << std::endl;
		QPhiX::qdp_pack_cb_spinor<>(chi, chi_qphix, *geom, 1);


		// Allocate the solutions
		QPhiX_Spinor** psi_qphix;
		psi_qphix=(QPhiX_Spinor**)ALIGNED_MALLOC(n_shift*sizeof(QPhiX_Spinor*),
				QPHIX_LLC_CACHE_ALIGN);
		if( psi_qphix == 0x0 ) {
			QDPIO::cout << "Unable toallocate psi_qphix" << std::endl;
			QDP_abort(1);
		}
		for(int s =0; s < n_shift; s++)  {
			psi_qphix[s] = (QPhiX_Spinor *)geom->allocCBFourSpinor();
			if( psi_qphix[s] == 0x0 ) {
				QDPIO::cout << "Unable to allocate psi_qphix["<<s<<"]" << std::endl;
				QDP_abort(1);
			}
		}

		// Initialize the solutions to zero.
		{
			QDPIO::cout << "Zeroing solutions" << std::endl;
			LatticeFermion zerotmp;
			zerotmp[ A->subset() ]  = zero;
			for(int s=0; s < n_shift; s++) {
				QPhiX::qdp_pack_cb_spinor<>(zerotmp,psi_qphix[s], *geom, 1);
			}
			// zerotmp goes away
		}

		multi1d<double> rsd_final(n_shift);
		multi1d<double> qphix_shifts(n_shift);
		for(int s=0; s < n_shift; s++) { qphix_shifts[s] = toDouble(shifts[s]); }

		multi1d<double> qphix_rsd_target(n_shift);

		if( invParam.RsdTarget.size() == n_shift ) {
			// There is a target for each pole. Just copy them over.
			for(int s=0; s < n_shift; s++) {
				qphix_rsd_target[s] = toDouble(invParam.RsdTarget[s]);
			}
		}
		else {

			if( invParam.RsdTarget.size() == 1 ) {
				// There is only 1 target, for all the poles. Replicate
				for(int s=0; s < n_shift; s++) {
					qphix_rsd_target[s] = toDouble(invParam.RsdTarget[0]);
				}
			}
			else{
				// Blow up!
				QDPIO::cout << "Size Mismatch between invParam.rsdTarget (size="
						<< invParam.RsdTarget.size()
						<< ") and n_shift=" << n_shift << std::endl;
				QDP_abort(1);
			}
		}

		// Zero counters
		unsigned long site_flops=0;
		unsigned long mv_apps=0;
		int my_isign=1;

		// Call solver
		double start = omp_get_wtime();
		(*mcg_solver)((QPhiX_Spinor **)psi_qphix,
				(const QPhiX_Spinor *)chi_qphix,
				(const int)n_shift,
				(const double*)qphix_shifts.slice(),
				(const double*)qphix_rsd_target.slice(),
				(int &)res.n_count,
				(double *)rsd_final.slice(),
				(unsigned long &)site_flops,
				(unsigned long &)mv_apps,
				(int)my_isign,
				(bool)invParam.VerboseP);
		double end = omp_get_wtime();

		// delete chi_qphix -- we are done with it.
		geom->free(chi_qphix);

		// Unpack solutions
		psi.resize(n_shift);
		for(int s=0; s < n_shift; s++) {
			QPhiX::qdp_unpack_cb_spinor<>(psi_qphix[s],psi[s], *geom, 1);
			geom->free(psi_qphix[s]); // done with it.
		}
		ALIGNED_FREE(psi_qphix); // done with it.

		// Report the residuum for the smallest shift
		// tho since other shifts can have looser criteria, this is
		// somewhat arbitrary
		int smallest = 0;
		for (int s=1; s < n_shift; s++) {
			if ( qphix_shifts[s] < qphix_shifts[smallest] ) smallest = s;
		}
		res.resid = rsd_final[smallest];

		QDPIO::cout << "QPHIX_RESIDUUM_REPORT:" << std::endl;
		for(int s=0; s < n_shift; s++) {
			QDPIO::cout << "\t Red_Final["<<s<<"]=" << rsd_final[s] <<std::endl;
		}

		// check results -- if native code is slow, this may be slow
		// so I give the option to turn it off...
		if ( invParam.SolutionCheckP ) {
			Double b2 = norm2(chi, A->subset());
			QDPIO::cout << "QPHIX_CLOVER_MULTI_SHIFT_CG_MDAGM_SOLVER: Residua Check: " << std::endl;
			for( int s=0; s < n_shift; s++) {
				T r = chi;
				T tmp,tmp2;

				(*A)(tmp, psi[s], PLUS);
				(*A)(tmp2, tmp, MINUS);
				tmp[A->subset()] = shifts[s]*psi[s] + tmp2; // Shift it.

				r[ A->subset()] -= tmp;
				Double r2 = norm2(r,A->subset());
				Double rel_resid = sqrt(r2/b2);
				QDPIO::cout << "\t shift["<<s<<"]  Actual || r || / || b || = "
						<< rel_resid << std::endl;
				if ( toDouble(rel_resid) > qphix_rsd_target[s]*toDouble(invParam.RsdToleranceFactor) )  {
					QDPIO::cout << "SOLVE FAILED: rel_resid=" << rel_resid
							<< " target=" << qphix_rsd_target[s]
															  << " tolerance_factor=" << invParam.RsdToleranceFactor
															  << " max tolerated=" << qphix_rsd_target[s]*toDouble(invParam.RsdToleranceFactor) << std::endl;
					QDP_abort(1);
				}
			}
		}

		int num_cb_sites = Layout::vol()/2;
		unsigned long total_flops = (site_flops + (1320+504+1320+504+48)*mv_apps)*num_cb_sites;
		double gflops = (double)(total_flops)/(1.0e9);

		double total_time = end - start;

		QDPIO::cout << "QPHIX_CLOVER_MULTI_SHIFT_CG_MDAGM_SOLVER: Iters=" << res.n_count << " Solver Time="<< total_time <<" (sec)  Performance=" << gflops / total_time << " GFLOPS" << std::endl;
		swatch.stop();
		QDPIO::cout << "QPHIX_CLOVER_MULTI_SHIFT_CG_MDAGM_SOLVER: total time: " << swatch.getTimeInSeconds() << " (sec)" << std::endl;
		END_CODE();
		return res;
	}

private:
	// Hide default constructor
	MdagMMultiSysSolverQPhiXClover() {}



	Handle< LinearOperator<T> > A;
	const MultiSysSolverQPhiXCloverParams invParam;
	Handle< CloverTermT<T, U> > clov;
	Handle< CloverTermT<T, U> > invclov;

	QPhiX::Geometry<REALT, VecTraits<REALT>::Vec, VecTraits<REALT>::Soa, VecTraits<REALT>::compress12>* geom;

	Handle< QPhiX::EvenOddCloverOperator<REALT, VecTraits<REALT>::Vec, VecTraits<REALT>::Soa, VecTraits<REALT>::compress12> > M;

	Handle< QPhiX::MInvCG<REALT,VecTraits<REALT>::Vec, VecTraits<REALT>::Soa, VecTraits<REALT>::compress12> > mcg_solver;


	QPhiX_Clover* invclov_packed[2];
	QPhiX_Clover* clov_packed[2];
	QPhiX_Gauge* u_packed[2];


};


} // End namespace

#endif // BUILD_QPHIX
#endif 


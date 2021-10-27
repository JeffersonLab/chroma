// -*- C++ -*-
/*! \file
 *  \QUDA MULTIGRID MdagM Clover solver.
 */

#ifndef __syssolver_mdagm_quda_multigrid_clover_h__
#define __syssolver_mdagm_quda_multigrid_clover_h__

#include "chroma_config.h"
#include "chromabase.h"
#include <cfloat>
#include <cstdio>

using namespace QDP;



#include "handle.h"
#include "state.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/fermbcs/simple_fermbc.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/invert/quda_solvers/syssolver_quda_multigrid_clover_params.h"
#include "actions/ferm/linop/clover_term_w.h"
#include "meas/gfix/temporal_gauge.h"
#include "io/aniso_io.h"
#include <string>
#include <sstream>

#include "lmdagm.h"
#include "util/gauge/reunit.h"
#include "actions/ferm/invert/quda_solvers/quda_mg_utils.h"
#include "actions/ferm/invert/mg_solver_exception.h"

//#include <util_quda.h>
#ifdef BUILD_QUDA
#include <quda.h>
#ifdef QDP_IS_QDPJIT
#include "actions/ferm/invert/quda_solvers/qdpjit_memory_wrapper.h"
#endif

#include "update/molecdyn/predictor/zero_guess_predictor.h"
#include "update/molecdyn/predictor/quda_predictor.h"
#include "meas/inline/io/named_objmap.h"

namespace Chroma
{

namespace MdagMSysSolverQUDAMULTIGRIDCloverEnv
{
//! Register the syssolver
bool registerAll();

}

class MdagMSysSolverQUDAMULTIGRIDClover : public MdagMSystemSolver<LatticeFermion>
{
public:
	typedef LatticeFermion T;
	typedef LatticeColorMatrix U;
	typedef multi1d<LatticeColorMatrix> Q;

	typedef LatticeFermionF TF;
	typedef LatticeColorMatrixF UF;
	typedef multi1d<LatticeColorMatrixF> QF;

	typedef LatticeFermionF TD;
	typedef LatticeColorMatrixF UD;
	typedef multi1d<LatticeColorMatrixF> QD;

	typedef WordType<T>::Type_t REALT;



	MdagMSysSolverQUDAMULTIGRIDClover(Handle< LinearOperator<T> > A_,
			Handle< FermState<T,Q,Q> > state_,
			const SysSolverQUDAMULTIGRIDCloverParams& invParam_) :
	  A(A_), gstate(state_), invParam(invParam_), clov(new CloverTermT<T, U>() ), invclov(new CloverTermT<T, U>())
	{
		StopWatch init_swatch;
		init_swatch.reset(); init_swatch.start();

		// Set the solver string
		{
			std::ostringstream solver_string_stream;
			solver_string_stream << "QUDA_MULTIGRID_CLOVER_MDAGM_SOLVER( Mass = " << invParam.CloverParams.Mass <<" , Id = "
					<< invParam.SaveSubspaceID << " ): ";
			solver_string = solver_string_stream.str();

		}
		QDPIO::cout << solver_string << "Initializing" << std::endl;

		// Check free mem
#if 0
		size_t free_mem = QUDAMGUtils::getCUDAFreeMem();
		QDPIO::cout << solver_string << "MEMCHECK: free mem = " << free_mem << std::endl;
#endif
		// FOLLOWING INITIALIZATION in test QUDA program

		// 1) work out cpu_prec, cuda_prec, cuda_prec_sloppy
		int s = sizeof( WordType<T>::Type_t );
		if (s == 4) {
			cpu_prec = QUDA_SINGLE_PRECISION;
		}
		else {
			cpu_prec = QUDA_DOUBLE_PRECISION;
		}

		// Work out GPU precision
		switch( invParam.cudaPrecision ) {
		case HALF:
			gpu_prec = QUDA_HALF_PRECISION;
			break;
		case SINGLE:
			gpu_prec = QUDA_SINGLE_PRECISION;
			break;
		case DOUBLE:
			gpu_prec = QUDA_DOUBLE_PRECISION;
			break;
		default:
			gpu_prec = cpu_prec;
			break;
		}

		// Work out GPU Sloppy precision
		// Default: No Sloppy
		switch( invParam.cudaSloppyPrecision ) {
		case HALF:
			gpu_half_prec = QUDA_HALF_PRECISION;
			break;
		case SINGLE:
			gpu_half_prec = QUDA_SINGLE_PRECISION;
			break;
		case DOUBLE:
			gpu_half_prec = QUDA_DOUBLE_PRECISION;
			break;
		default:
			gpu_half_prec = gpu_prec;
			break;
		}

		// 2) pull 'new; GAUGE and Invert params
		q_gauge_param = newQudaGaugeParam();
		quda_inv_param = newQudaInvertParam();

		// 3) set lattice size
		const multi1d<int>& latdims = Layout::subgridLattSize();

		q_gauge_param.X[0] = latdims[0];
		q_gauge_param.X[1] = latdims[1];
		q_gauge_param.X[2] = latdims[2];
		q_gauge_param.X[3] = latdims[3];

		// 4) - deferred (anisotropy)

		// 5) - set QUDA_WILSON_LINKS, QUDA_GAUGE_ORDER
		q_gauge_param.type = QUDA_WILSON_LINKS;
#ifndef BUILD_QUDA_DEVIFACE_GAUGE
		q_gauge_param.gauge_order = QUDA_QDP_GAUGE_ORDER; // gauge[mu], p
#else
		q_gauge_param.location = QUDA_CUDA_FIELD_LOCATION;
		q_gauge_param.gauge_order = QUDA_QDPJIT_GAUGE_ORDER;
#endif

		// 6) - set t_boundary
		// Convention: BC has to be applied already
		// This flag just tells QUDA that this is so,
		// so that QUDA can take care in the reconstruct
		if( invParam.AntiPeriodicT ) {
			q_gauge_param.t_boundary = QUDA_ANTI_PERIODIC_T;
		}
		else {
			q_gauge_param.t_boundary = QUDA_PERIODIC_T;
		}

		// Set cpu_prec, cuda_prec, reconstruct and sloppy versions
		q_gauge_param.cpu_prec = cpu_prec;
		q_gauge_param.cuda_prec = gpu_prec;

		switch( invParam.cudaReconstruct ) {
		case RECONS_NONE:
			q_gauge_param.reconstruct = QUDA_RECONSTRUCT_NO;
			break;
		case RECONS_8:
			q_gauge_param.reconstruct = QUDA_RECONSTRUCT_8;
			break;
		case RECONS_12:
			q_gauge_param.reconstruct = QUDA_RECONSTRUCT_12;
			break;
		default:
			q_gauge_param.reconstruct = QUDA_RECONSTRUCT_12;
			break;
		};

		q_gauge_param.cuda_prec_sloppy = gpu_half_prec;
		q_gauge_param.cuda_prec_precondition = gpu_half_prec;

		switch( invParam.cudaSloppyReconstruct ) {
		case RECONS_NONE:
			q_gauge_param.reconstruct_sloppy = QUDA_RECONSTRUCT_NO;
			break;
		case RECONS_8:
			q_gauge_param.reconstruct_sloppy = QUDA_RECONSTRUCT_8;
			break;
		case RECONS_12:
			q_gauge_param.reconstruct_sloppy = QUDA_RECONSTRUCT_12;
			break;
		default:
			q_gauge_param.reconstruct_sloppy = QUDA_RECONSTRUCT_12;
			break;
		};

		// This may be overridden later
		q_gauge_param.reconstruct_precondition=q_gauge_param.reconstruct_sloppy;

		// Gauge fixing:

		// These are the links
		// They may be smeared and the BC's may be applied
		Q links_single(Nd);

		// Now downcast to single prec fields.
		for(int mu=0; mu < Nd; mu++) {
			links_single[mu] = (state_->getLinks())[mu];
		}

		// GaugeFix
		if( invParam.axialGaugeP ) {
			temporalGauge(links_single, GFixMat, Nd-1);
			for(int mu=0; mu < Nd; mu++) {
				links_single[mu] = GFixMat*(state_->getLinks())[mu]*adj(shift(GFixMat, FORWARD, mu));
			}
			q_gauge_param.gauge_fix = QUDA_GAUGE_FIXED_YES;
		}
		else {
			// No GaugeFix
			q_gauge_param.gauge_fix = QUDA_GAUGE_FIXED_NO;// No Gfix yet
		}

		// deferred 4) Gauge Anisotropy
		const AnisoParam_t& aniso = invParam.CloverParams.anisoParam;
		if( aniso.anisoP ) {                     // Anisotropic case
			Real gamma_f = aniso.xi_0 / aniso.nu;
			q_gauge_param.anisotropy = toDouble(gamma_f);
		}
		else {
			q_gauge_param.anisotropy = 1.0;
		}

		// MAKE FSTATE BEFORE RESCALING links_single
		// Because the clover term expects the unrescaled links...
		Handle<FermState<T,Q,Q> > fstate( new PeriodicFermState<T,Q,Q>(links_single));

		if( aniso.anisoP ) {                     // Anisotropic case
			multi1d<Real> cf=makeFermCoeffs(aniso);
			for(int mu=0; mu < Nd; mu++) {
				links_single[mu] *= cf[mu];
			}
		}

		// Now onto the inv param:
		// Dslash type
		quda_inv_param.dslash_type = QUDA_CLOVER_WILSON_DSLASH;

		// Hardwire to GCR
		quda_inv_param.inv_type = QUDA_GCR_INVERTER;
		quda_inv_param.compute_true_res = 0;

		quda_inv_param.kappa = 0.5;
		quda_inv_param.clover_coeff = 1.0; // Dummy, not used
                quda_inv_param.clover_csw = 0.0;   // Use zero and QUDA will ignore it 
		quda_inv_param.Ls=1;
		quda_inv_param.tol = toDouble(invParam.RsdTarget);
		quda_inv_param.maxiter = invParam.MaxIter;
		quda_inv_param.reliable_delta = toDouble(invParam.Delta);
		quda_inv_param.pipeline = invParam.Pipeline;

		quda_inv_param.solution_type = QUDA_MATPC_SOLUTION;
		quda_inv_param.solve_type = QUDA_DIRECT_PC_SOLVE;

		
		quda_inv_param.matpc_type = QUDA_MATPC_ODD_ODD;

		quda_inv_param.dagger = QUDA_DAG_NO;
		quda_inv_param.mass_normalization = QUDA_KAPPA_NORMALIZATION;

		quda_inv_param.cpu_prec = cpu_prec;
		quda_inv_param.cuda_prec = gpu_prec;
		quda_inv_param.cuda_prec_sloppy = gpu_half_prec;
		quda_inv_param.cuda_prec_precondition = gpu_half_prec;
		quda_inv_param.preserve_source = QUDA_PRESERVE_SOURCE_NO;
		quda_inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;

#ifndef BUILD_QUDA_DEVIFACE_SPINOR
		quda_inv_param.dirac_order = QUDA_DIRAC_ORDER;
		quda_inv_param.input_location = QUDA_CPU_FIELD_LOCATION;
		quda_inv_param.output_location = QUDA_CPU_FIELD_LOCATION;

#else
		quda_inv_param.dirac_order = QUDA_QDPJIT_DIRAC_ORDER;
		quda_inv_param.input_location = QUDA_CUDA_FIELD_LOCATION;
		quda_inv_param.output_location = QUDA_CUDA_FIELD_LOCATION;
#endif

		// Autotuning
		if( invParam.tuneDslashP ) {
			quda_inv_param.tune = QUDA_TUNE_YES;
		}
		else {
			quda_inv_param.tune = QUDA_TUNE_NO;
		}

		// Setup padding
		multi1d<int> face_size(4);
		face_size[0] = latdims[1]*latdims[2]*latdims[3]/2;
		face_size[1] = latdims[0]*latdims[2]*latdims[3]/2;
		face_size[2] = latdims[0]*latdims[1]*latdims[3]/2;
		face_size[3] = latdims[0]*latdims[1]*latdims[2]/2;

		int max_face = face_size[0];
		for(int i=1; i <=3; i++) {
			if ( face_size[i] > max_face ) {
				max_face = face_size[i];
			}
		}

		q_gauge_param.ga_pad = max_face;
		// PADDING
		quda_inv_param.sp_pad = 0;
		quda_inv_param.cl_pad = 0;

		// Clover precision and order
		quda_inv_param.clover_cpu_prec = cpu_prec;
		quda_inv_param.clover_cuda_prec = gpu_prec;
		quda_inv_param.clover_cuda_prec_sloppy = gpu_half_prec;
		quda_inv_param.clover_cuda_prec_precondition = gpu_half_prec;

		if( !invParam.MULTIGRIDParamsP )  {
			QDPIO::cout << solver_string << "ERROR: MG Solver had MULTIGRIDParamsP set to false" << std::endl;
			QDP_abort(1);
		}

		// Dereference handle
		const MULTIGRIDSolverParams& ip = *(invParam.MULTIGRIDParams);


		// Set preconditioner precision
		switch( ip.prec ) {
		case HALF:
			quda_inv_param.cuda_prec_precondition = QUDA_HALF_PRECISION;
			quda_inv_param.clover_cuda_prec_precondition = QUDA_HALF_PRECISION;
			q_gauge_param.cuda_prec_precondition = QUDA_HALF_PRECISION;
			break;

		case SINGLE:
			quda_inv_param.cuda_prec_precondition = QUDA_SINGLE_PRECISION;
			quda_inv_param.clover_cuda_prec_precondition = QUDA_SINGLE_PRECISION;
			q_gauge_param.cuda_prec_precondition = QUDA_SINGLE_PRECISION;
			break;

		case DOUBLE:
			quda_inv_param.cuda_prec_precondition = QUDA_DOUBLE_PRECISION;
			quda_inv_param.clover_cuda_prec_precondition = QUDA_DOUBLE_PRECISION;
			q_gauge_param.cuda_prec_precondition = QUDA_DOUBLE_PRECISION;
			break;
		default:
			quda_inv_param.cuda_prec_precondition = QUDA_HALF_PRECISION;
			quda_inv_param.clover_cuda_prec_precondition = QUDA_HALF_PRECISION;
			q_gauge_param.cuda_prec_precondition = QUDA_HALF_PRECISION;
			break;
		}

		switch( ip.reconstruct ) {
		case RECONS_NONE:
			q_gauge_param.reconstruct_precondition = QUDA_RECONSTRUCT_NO;
			break;
		case RECONS_8:
			q_gauge_param.reconstruct_precondition = QUDA_RECONSTRUCT_8;
			break;
		case RECONS_12:
			q_gauge_param.reconstruct_precondition = QUDA_RECONSTRUCT_12;
			break;
		default:
			q_gauge_param.reconstruct_precondition = QUDA_RECONSTRUCT_12;
			break;
		};

		// Set up the links
		void* gauge[4];

#ifndef BUILD_QUDA_DEVIFACE_GAUGE
		for(int mu=0; mu < Nd; mu++) {
			gauge[mu] = (void *)&(links_single[mu].elem(all.start()).elem().elem(0,0).real());
		}
#else
		GetMemoryPtrGauge(gauge,links_single);
		// std::vector<int> ids;
		// for(int mu=0; mu < Nd; mu++) 
		//   ids.push_back( links_single[mu].getId() );
		// std::vector<void*> dev_ptr = GetMemoryPtr( ids );
		// for(int mu=0; mu < Nd; mu++) 
		//   gauge[mu] = dev_ptr[mu];
#endif

		loadGaugeQuda((void *)gauge, &q_gauge_param);


		quda_inv_param.tol_precondition = toDouble(ip.tol[0]);
		quda_inv_param.maxiter_precondition = ip.maxIterations[0];
		quda_inv_param.gcrNkrylov = ip.outer_gcr_nkrylov;
		quda_inv_param.residual_type = static_cast<QudaResidualType>(QUDA_L2_RELATIVE_RESIDUAL);

		//Replacing above with what's in the invert test.
		switch( ip.schwarzType ) {
		case ADDITIVE_SCHWARZ :
			quda_inv_param.schwarz_type = QUDA_ADDITIVE_SCHWARZ;
			break;
		case MULTIPLICATIVE_SCHWARZ :
			quda_inv_param.schwarz_type = QUDA_MULTIPLICATIVE_SCHWARZ;
			break;
		default:
			quda_inv_param.schwarz_type = QUDA_ADDITIVE_SCHWARZ;
			break;
		}
		quda_inv_param.precondition_cycle = 1;
		//Invert test always sets this to 1.


		if( invParam.verboseP ) {
			quda_inv_param.verbosity = QUDA_VERBOSE;
		}
		else {
			quda_inv_param.verbosity = QUDA_SUMMARIZE;
		}

		quda_inv_param.verbosity_precondition = QUDA_SILENT;

		quda_inv_param.inv_type_precondition = QUDA_MG_INVERTER;

		//      Setup the clover term...
		QDPIO::cout <<solver_string << "Creating CloverTerm" << std::endl;
		clov->create(fstate, invParam_.CloverParams);

		// Don't recompute, just copy
		invclov->create(fstate, invParam_.CloverParams);

		QDPIO::cout <<solver_string<< "Inverting CloverTerm" << std::endl;
		invclov->choles(0);
		invclov->choles(1);

#ifndef BUILD_QUDA_DEVIFACE_CLOVER
#warning "NOT USING QUDA DEVICE IFACE"
		quda_inv_param.clover_order = QUDA_PACKED_CLOVER_ORDER;

		multi1d<QUDAPackedClovSite<REALT> > packed_clov;

		packed_clov.resize(all.siteTable().size());

		clov->packForQUDA(packed_clov, 0);
		clov->packForQUDA(packed_clov, 1);

		// Always need inverse
		multi1d<QUDAPackedClovSite<REALT> > packed_invclov(all.siteTable().size());
		invclov->packForQUDA(packed_invclov, 0);
		invclov->packForQUDA(packed_invclov, 1);

		loadCloverQuda(&(packed_clov[0]), &(packed_invclov[0]), &quda_inv_param);

#else

#warning "USING QUDA DEVICE IFACE"

		quda_inv_param.clover_location = QUDA_CUDA_FIELD_LOCATION;
		quda_inv_param.clover_order = QUDA_QDPJIT_CLOVER_ORDER;

		void *clover[2];
		void *cloverInv[2];

		GetMemoryPtrClover(clov->getOffId(),clov->getDiaId(),invclov->getOffId(),invclov->getDiaId());

		loadCloverQuda( (void*)(clover), (void *)(cloverInv), &quda_inv_param);
#endif

		quda_inv_param.omega = toDouble(ip.relaxationOmegaOuter);

// Copy ThresholdCount from invParams into threshold_counts.
threshold_counts = invParam.ThresholdCount;

if(TheNamedObjMap::Instance().check(invParam.SaveSubspaceID))
{
	StopWatch update_swatch;
	update_swatch.reset(); update_swatch.start();
	// Subspace ID exists add it to mg_state
	QDPIO::cout<< solver_string <<"Recovering subspace..."<<std::endl;
	subspace_pointers = TheNamedObjMap::Instance().getData< QUDAMGUtils::MGSubspacePointers* >(invParam.SaveSubspaceID);
	for(int j=0; j < ip.mg_levels-1;++j) {
		(subspace_pointers->mg_param).setup_maxiter_refresh[j] = 0;
	}
	updateMultigridQuda(subspace_pointers->preconditioner, &(subspace_pointers->mg_param));
	update_swatch.stop();

	QDPIO::cout << solver_string << " subspace_update_time = "
			<< update_swatch.getTimeInSeconds() << " sec. " << std::endl;
}
else
{
	// Create the subspace.
	StopWatch create_swatch;
	create_swatch.reset(); create_swatch.start();
	QDPIO::cout << solver_string << "Creating Subspace" << std::endl;
	subspace_pointers = QUDAMGUtils::create_subspace<T>(invParam);
	XMLBufferWriter file_xml;
	push(file_xml, "FileXML");
	pop(file_xml);

	int foo = 5;

	XMLBufferWriter record_xml;
	push(record_xml, "RecordXML");
	write(record_xml, "foo", foo);
	pop(record_xml);


	TheNamedObjMap::Instance().create< QUDAMGUtils::MGSubspacePointers* >(invParam.SaveSubspaceID);
	TheNamedObjMap::Instance().get(invParam.SaveSubspaceID).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(invParam.SaveSubspaceID).setRecordXML(record_xml);

	TheNamedObjMap::Instance().getData< QUDAMGUtils::MGSubspacePointers* >(invParam.SaveSubspaceID) = subspace_pointers;
	create_swatch.stop();
	QDPIO::cout << solver_string << " subspace_create_time = "
			<< create_swatch.getTimeInSeconds() << " sec. " << std::endl;

}
quda_inv_param.preconditioner = subspace_pointers->preconditioner;

init_swatch.stop();
QDPIO::cout << solver_string << " init_time = "
		<< init_swatch.getTimeInSeconds() << " sec. "
		<< std::endl;

	}

	//! Destructor is not automatic
	~MdagMSysSolverQUDAMULTIGRIDClover()
	{

		quda_inv_param.preconditioner = nullptr;
		subspace_pointers = nullptr;
		freeGaugeQuda();
		freeCloverQuda();
	}

	//! Return the subset on which the operator acts
	const Subset& subset() const {return A->subset();}

	//! Solver the linear system
	/*!
	 * \param psi      solution ( Modify )
	 * \param chi      source ( Read )
	 * \return syssolver results
	 */
	SystemSolverResults_t operator() (T& psi, const T& chi) const
	{
		SystemSolverResults_t res1;
		SystemSolverResults_t res2;
		SystemSolverResults_t res;

		START_CODE();
		StopWatch swatch;
		swatch.start();

		// I want to use the predictor versions of the code as they have been made robust.
		// So I should use either a null predictor or a zero guess predictor here.
		// The MG two step solve logic is quite complicated and may need to reinit the fields.
		// I don't want to triplicate that logic so I'll just use a dummy predictor and call through.
		ZeroGuess4DChronoPredictor dummy_predictor;
		res = (*this)(psi, chi, dummy_predictor);


		END_CODE();
		return res;
	}



	SystemSolverResults_t operator() (T& psi, const T& chi, Chroma::AbsTwoStepChronologicalPredictor4D<T>& predictor ) const
	{

		START_CODE();

		StopWatch swatch;
		swatch.start();

		MULTIGRIDSolverParams& ip = *(invParam.MULTIGRIDParams);
		// Use this in residuum checks.
		Double norm2chi=sqrt(norm2(chi, A->subset()));

		// Allow QUDA to use initial guess
		QudaUseInitGuess old_guess_policy = quda_inv_param.use_init_guess;
		quda_inv_param.use_init_guess = QUDA_USE_INIT_GUESS_YES;


		SystemSolverResults_t res;
		SystemSolverResults_t res1;
		SystemSolverResults_t res2;

		// Create MdagM op
		Handle< LinearOperator<T> > MdagM( new MdagMLinOp<T>(A) );


		QDPIO::cout << solver_string <<"Two Step Solve" << std::endl;


		// Try to cast the predictor to a two step predictor
		StopWatch X_prediction_timer; X_prediction_timer.reset();
		StopWatch Y_prediction_timer; Y_prediction_timer.reset();
		StopWatch Y_solve_timer;      Y_solve_timer.reset();
		StopWatch X_solve_timer;      X_solve_timer.reset();
		StopWatch Y_predictor_add_timer; Y_predictor_add_timer.reset();
		StopWatch X_predictor_add_timer; X_predictor_add_timer.reset();
		StopWatch X_refresh_timer; X_refresh_timer.reset();
		StopWatch Y_refresh_timer; Y_refresh_timer.reset();

		QDPIO::cout << solver_string << "Predicting Y" << std::endl;
		Y_prediction_timer.start();
		T Y_prime = zero;
		{
			T tmp_vec = psi;
			predictor.predictY(tmp_vec, *A, chi); // Predicts for M^\dagger Y = chi

			// We are going to solve M \gamma
			Y_prime = Gamma(Nd*Nd-1)*tmp_vec;
		}
		Y_prediction_timer.stop();

		// Y solve: M^\dagger Y = chi
		//        g_5 M g_5 Y = chi
		//     =>    M Y' = chi'  with chi' = gamma_5*chi

		Y_solve_timer.start();


		T g5chi = zero;
		T Y = zero;
		g5chi[rb[1]]= Gamma(Nd*Nd-1)*chi;

		// Y solve at 0.5 * Target Residuum -- Evan's bound
		quda_inv_param.tol = toDouble(Real(0.5)*invParam.RsdTarget);
		if( invParam.asymmetricP == true ) {
			res1 = qudaInvert(*clov,
					*invclov,
					g5chi,
					Y_prime);
			Y[rb[1]] = Gamma(Nd*Nd -1)*Y_prime;
		}
		else {
			T tmp = zero;
			invclov->apply(tmp,g5chi,MINUS,1);

			res1 = qudaInvert(*clov,
						*invclov,
						tmp,
						Y_prime);
#ifdef QUDA_DEBUG
			{
			   char Y_prime_norm[256];
 			   char Y_prime_norm_full[256];
			   std::sprintf(Y_prime_norm, "%.*e", DECIMAL_DIG, toDouble(norm2(Y_prime, A->subset())));
			   std::sprintf(Y_prime_norm_full, "%.*e", DECIMAL_DIG, toDouble(norm2(Y_prime)));
			   QDPIO::cout << "Y solution: norm2(subset) = " << Y_prime_norm << " norm(full) = " << Y_prime_norm_full << std::endl;
		        }
#endif
			tmp[rb[1]] = Gamma(Nd*Nd-1)*Y_prime;
			clov->apply(Y,tmp,MINUS,1);

		}


		bool solution_good = true;

		// Check solution
		{
			T r=zero;
			r[A->subset()]=chi;
			T tmp;
			(*A)(tmp, Y, MINUS);
			r[A->subset()] -= tmp;

			res1.resid = sqrt(norm2(r, A->subset()));
			QDPIO::cout << "Y-solve: ||r||=" << res1.resid << "   ||r||/||b||="
					<< res1.resid/sqrt(norm2(chi,rb[1])) << std::endl;
			if ( toBool( res1.resid/norm2chi > invParam.RsdToleranceFactor * invParam.RsdTarget ) ) {
				solution_good = false;
			}
		}

		if ( solution_good ) {
			if( res1.n_count >= threshold_counts ) {
				QDPIO::cout << solver_string << "Iteration Threshold Exceeded! Y Solver iters = " << res1.n_count << " Threshold=" << threshold_counts << std::endl;
				QDPIO::cout << solver_string << "Refreshing Subspace" << std::endl;

				Y_refresh_timer.start();
				// refresh the subspace
				// Setup the number of subspace Iterations
				for(int j=0; j < ip.mg_levels-1; j++) {
					(subspace_pointers->mg_param).setup_maxiter_refresh[j] = ip.maxIterSubspaceRefresh[j];
				}
				updateMultigridQuda(subspace_pointers->preconditioner, &(subspace_pointers->mg_param));
				for(int j=0; j < ip.mg_levels-1; j++) {
					(subspace_pointers->mg_param).setup_maxiter_refresh[j] = 0;
				}
				Y_refresh_timer.stop();
				QDPIO::cout << solver_string << "Subspace Refresh Time = " << Y_refresh_timer.getTimeInSeconds() << " secs\n";
			}
		}
		else {
			QDPIO::cout << solver_string << "Y-Solve failed (seq: "<< seqno <<"). Blowing away and reiniting subspace" << std::endl;
			StopWatch reinit_timer; reinit_timer.reset();
			reinit_timer.start();

			// Delete the saved subspace completely
			QUDAMGUtils::delete_subspace(invParam.SaveSubspaceID);

			// Recreate the subspace
			bool saved_value = ip.check_multigrid_setup;
			ip.check_multigrid_setup = true;
			subspace_pointers = QUDAMGUtils::create_subspace<T>(invParam);
			ip.check_multigrid_setup = saved_value;

			// Make subspace XML snippets
			XMLBufferWriter file_xml;
			push(file_xml, "FileXML");
			pop(file_xml);

			int foo = 5;
			XMLBufferWriter record_xml;
			push(record_xml, "RecordXML");
			write(record_xml, "foo", foo);
			pop(record_xml);


			// Create named object entry.
			TheNamedObjMap::Instance().create< QUDAMGUtils::MGSubspacePointers* >(invParam.SaveSubspaceID);
			TheNamedObjMap::Instance().get(invParam.SaveSubspaceID).setFileXML(file_xml);
			TheNamedObjMap::Instance().get(invParam.SaveSubspaceID).setRecordXML(record_xml);

			// Assign the pointer into the named object
			TheNamedObjMap::Instance().getData< QUDAMGUtils::MGSubspacePointers* >(invParam.SaveSubspaceID) = subspace_pointers;
			quda_inv_param.preconditioner = subspace_pointers->preconditioner;

			reinit_timer.stop();
			QDPIO::cout << solver_string << "Subspace Reinit Time: " << reinit_timer.getTimeInSeconds() << " sec."  << std::endl;

			// Re-solve
			QDPIO::cout << solver_string << "Re-Solving for Y with zero guess" << std::endl;
			SystemSolverResults_t res_tmp;

			Y_prime = zero;
			if( invParam.asymmetricP == true ) {
				res_tmp = qudaInvert(*clov,
						*invclov,
						g5chi,
						Y_prime);
				Y[rb[1]] = Gamma(Nd*Nd -1)*Y_prime;
			}
			else {
				T tmp = zero;
				invclov->apply(tmp,g5chi,MINUS,1);
				res_tmp = qudaInvert(*clov,
							*invclov,
							tmp,
							Y_prime);

#ifdef QUDA_DEBUG
			{
                           char Y_prime_norm[256];
                           char Y_prime_norm_full[256];
                           std::sprintf(Y_prime_norm, "%.*e", DECIMAL_DIG, toDouble(norm2(Y_prime, A->subset())));
                           std::sprintf(Y_prime_norm_full, "%.*e", DECIMAL_DIG, toDouble(norm2(Y_prime)));
                           QDPIO::cout << "Y solution: norm2(subset) = " << Y_prime_norm << " norm(full) = " << Y_prime_norm_full << std::endl;
                        }
#endif

				tmp[rb[1]] = Gamma(Nd*Nd-1)*Y_prime;
				clov->apply(Y,tmp,MINUS,1);

			}

			// Check solution
			{
				T r=zero;
				r[A->subset()]=chi;
				T tmp;
				(*A)(tmp, Y,MINUS);
				r[A->subset()] -= tmp;

				res_tmp.resid = sqrt(norm2(r, A->subset()));
				if ( toBool( res_tmp.resid/sqrt(norm2(chi)) > invParam.RsdToleranceFactor * invParam.RsdTarget ) ) {
					QDPIO::cout << solver_string << "Re Solve for Y Failed (seq: " << seqno << " ) Rsd = " << res_tmp.resid/norm2chi << " RsdTarget = " << invParam.RsdTarget << std::endl;
					QDPIO::cout << solver_string << "Throwing Exception! This will ABORT" << std::endl;

					dumpYSolver(g5chi,Y_prime);

					MGSolverException convergence_fail(invParam.CloverParams.Mass,
							invParam.SaveSubspaceID,
							res_tmp.n_count,
							Real(res_tmp.resid/norm2chi),
							invParam.RsdTarget*invParam.RsdToleranceFactor);
					throw convergence_fail;
				}
			} // Check solution

			// threhold count is good, and solution is good
			res1.n_count += res_tmp.n_count; // Add resolve iterations
			res1.resid  = res_tmp.resid; // Copy new residuum.

		}
		Y_solve_timer.stop();

		// At this point we should have a good solution.
		Y_predictor_add_timer.start();
		predictor.newYVector(Y);
		Y_predictor_add_timer.stop();

		X_prediction_timer.start();
		// Can predict psi in the usual way without reference to Y
		predictor.predictX(psi, (*MdagM), chi);
		X_prediction_timer.stop();

		// Restore resid target for X solve
		quda_inv_param.tol = toDouble(invParam.RsdTarget);
		X_solve_timer.start();
		// Solve for psi
		res2 = qudaInvert(*clov,
				*invclov,
				Y,
				psi);
#ifdef QUDA_DEBUG
			{
                           char X_prime_norm[256];
                           char X_prime_norm_full[256];
                           std::sprintf(X_prime_norm, "%.*e", DECIMAL_DIG, toDouble(norm2(psi, A->subset())));
                           std::sprintf(X_prime_norm_full, "%.*e", DECIMAL_DIG, toDouble(norm2(psi)));
                           QDPIO::cout << "X solution: norm2(subset) = " << X_prime_norm << " norm(full) = " << X_prime_norm_full << std::endl;
                        }
#endif
		solution_good = true;

		// Check solution
		{
			T r;
			r[A->subset()]=chi;
			T tmp;
			(*MdagM)(tmp, psi, PLUS);
			r[A->subset()] -= tmp;

			res2.resid = sqrt(norm2(r, A->subset()));
			if ( toBool( res2.resid/norm2chi > invParam.RsdToleranceFactor * invParam.RsdTarget ) ) {
				solution_good = false;
			}
		}

		if( solution_good )  {
			if( res2.n_count >= threshold_counts ) {
				QDPIO::cout << solver_string <<"Threshold Reached! X Solver iters = " << res2.n_count << " Threshold=" << threshold_counts << std::endl;
				QDPIO::cout << solver_string << "Refreshing Subspace" << std::endl;

				X_refresh_timer.start();
				// refresh the subspace
				// Regenerate space. Destroy and recreate
				// Setup the number of subspace Iterations
				for(int j=0; j < ip.mg_levels-1; j++) {
					(subspace_pointers->mg_param).setup_maxiter_refresh[j] = ip.maxIterSubspaceRefresh[j];
				}
				updateMultigridQuda(subspace_pointers->preconditioner, &(subspace_pointers->mg_param));
				for(int j=0; j < ip.mg_levels-1; j++) {
					(subspace_pointers->mg_param).setup_maxiter_refresh[j] = 0;
				}
				X_refresh_timer.stop();

				QDPIO::cout << solver_string << "X Subspace Refresh Time = " << X_refresh_timer.getTimeInSeconds() << " secs\n";
			}
		}
		else {

			QDPIO::cout << solver_string << "X-Solve failed (seq: "<<seqno<<") . Blowing away and reiniting subspace" << std::endl;
			StopWatch reinit_timer; reinit_timer.reset();
			reinit_timer.start();

			// Delete the saved subspace completely
			QUDAMGUtils::delete_subspace(invParam.SaveSubspaceID);

			// Recreate the subspace
			bool saved_value = ip.check_multigrid_setup;
                        ip.check_multigrid_setup = true;
                        subspace_pointers = QUDAMGUtils::create_subspace<T>(invParam);
                        ip.check_multigrid_setup = saved_value;


			// Make subspace XML snippets
			XMLBufferWriter file_xml;
			push(file_xml, "FileXML");
			pop(file_xml);

			int foo = 5;
			XMLBufferWriter record_xml;
			push(record_xml, "RecordXML");
			write(record_xml, "foo", foo);
			pop(record_xml);


			// Create named object entry.
			TheNamedObjMap::Instance().create< QUDAMGUtils::MGSubspacePointers* >(invParam.SaveSubspaceID);
			TheNamedObjMap::Instance().get(invParam.SaveSubspaceID).setFileXML(file_xml);
			TheNamedObjMap::Instance().get(invParam.SaveSubspaceID).setRecordXML(record_xml);

			// Assign the pointer into the named object
			TheNamedObjMap::Instance().getData< QUDAMGUtils::MGSubspacePointers* >(invParam.SaveSubspaceID) = subspace_pointers;
			quda_inv_param.preconditioner = subspace_pointers->preconditioner;
			reinit_timer.stop();
			QDPIO::cout << solver_string << "Subspace Reinit Time: " << reinit_timer.getTimeInSeconds() << " sec."  << std::endl;

			// Re-solve
			QDPIO::cout << solver_string << "Re-Solving for X with zero guess" << std::endl;
			SystemSolverResults_t res_tmp;
			psi = zero;
			res_tmp = qudaInvert(*clov,
					*invclov,
					Y,
					psi);
#ifdef QUDA_DEBUG
                        {
                           char X_prime_norm[256];
                           char X_prime_norm_full[256];
                           std::sprintf(X_prime_norm, "%.*e", DECIMAL_DIG, toDouble(norm2(psi, A->subset())));
                           std::sprintf(X_prime_norm_full, "%.*e", DECIMAL_DIG, toDouble(norm2(psi)));
                           QDPIO::cout << "X solution: norm2(subset) = " << X_prime_norm << " norm(full) = " << X_prime_norm_full << std::endl;
                        }
#endif


			// Check solution
			{
				T r;
				r[A->subset()]=chi;
				T tmp;
				(*MdagM)(tmp, psi, PLUS);
				r[A->subset()] -= tmp;

				res_tmp.resid = sqrt(norm2(r, A->subset()));
				if ( toBool( res_tmp.resid/norm2chi > invParam.RsdToleranceFactor * invParam.RsdTarget ) ) {
					QDPIO::cout << solver_string << "Re Solve for X Failed (seq: " << seqno << " ) Rsd = " << res_tmp.resid/norm2chi << " RsdTarget = " << invParam.RsdTarget << std::endl;
					QDPIO::cout << solver_string << "Throwing Exception! This will ABORT" << std::endl;

					dumpXSolver(chi,Y,psi);

					MGSolverException convergence_fail(invParam.CloverParams.Mass,
							invParam.SaveSubspaceID,
							res_tmp.n_count,
							Real(res_tmp.resid/norm2chi),
							invParam.RsdTarget*invParam.RsdToleranceFactor);
					throw convergence_fail;

					QDP_abort(1);
				}
			}
			// At this point the solution is good
			res2.n_count += res_tmp.n_count;
			res2.resid = res_tmp.resid;

		}
		X_solve_timer.stop();

		X_predictor_add_timer.start();
		predictor.newXVector(psi);
		X_predictor_add_timer.stop();
		swatch.stop();
		double time = swatch.getTimeInSeconds();

		res.n_count = res1.n_count + res2.n_count;
		res.resid = res2.resid;

		Double rel_resid = res.resid/norm2chi;

		QDPIO::cout <<  solver_string   << " seq: " << (seqno++) << " iterations: " << res1.n_count << " + "
				<<  res2.n_count << " = " << res.n_count
				<<  " Rsd = " << res.resid << " Relative Rsd = " << rel_resid << std::endl;

		QDPIO::cout <<solver_string  << "Time: Y_predict: " << Y_prediction_timer.getTimeInSeconds() << " (s) "
				<< "Y_solve: " << Y_solve_timer.getTimeInSeconds() << " (s) "
				<< "Y_register: " << Y_predictor_add_timer.getTimeInSeconds() << " (s) "
				<< "X_predict: " << X_prediction_timer.getTimeInSeconds() << " (s) "
				<< "X_solve: " << X_solve_timer.getTimeInSeconds() << " (s) "
				<< "X_register: " << X_predictor_add_timer.getTimeInSeconds() << " (s) "
				<< "Total time: " << time <<  "(s)" << std::endl;

		quda_inv_param.use_init_guess = old_guess_policy;

		return res;
	}

	SystemSolverResults_t operator() (T& psi, const T& chi, Chroma::QUDA4DChronoPredictor& predictor ) const
	{


		START_CODE();

		StopWatch swatch;
		swatch.start();

		MULTIGRIDSolverParams& ip = *(invParam.MULTIGRIDParams);
		// Use this in residuum checks.
		Double norm2chi=sqrt(norm2(chi, A->subset()));

		// Allow QUDA to use initial guess
		QudaUseInitGuess old_guess_policy = quda_inv_param.use_init_guess;
		quda_inv_param.use_init_guess = QUDA_USE_INIT_GUESS_YES;


		SystemSolverResults_t res;
		SystemSolverResults_t res1;
		SystemSolverResults_t res2;

		// Create MdagM op
		Handle< LinearOperator<T> > MdagM( new MdagMLinOp<T>(A) );


		QDPIO::cout << solver_string <<"Two Step Solve" << std::endl;


		// Try to cast the predictor to a two step predictor
		StopWatch Y_solve_timer;      Y_solve_timer.reset();
		StopWatch X_solve_timer;      X_solve_timer.reset();
		StopWatch Y_refresh_timer;    Y_refresh_timer.reset();
		StopWatch X_refresh_timer;    X_refresh_timer.reset();
		int X_index=predictor.getXIndex();
		int Y_index=predictor.getYIndex();

		QDPIO::cout << "Two Step Solve using QUDA predictor: (Y_index,X_index) = ( " << Y_index << " , " << X_index << " ) \n";


		// Select the channel for QUDA's predictor here.
		//
		//

		quda_inv_param.chrono_max_dim = predictor.getMaxChrono();
		quda_inv_param.chrono_index = Y_index;
		quda_inv_param.chrono_make_resident = true;
		quda_inv_param.chrono_use_resident = true;
		quda_inv_param.chrono_replace_last = false;

		// Y solve is at 0.5*RsdTarget  -- Evan's analysis
		quda_inv_param.tol = toDouble(Real(0.5)*invParam.RsdTarget);
		if ( predictor.getChronoPrecision() == DEFAULT ) {
			QDPIO::cout << "Setting Default Chrono precision of " << cpu_prec << std::endl;
			quda_inv_param.chrono_precision = cpu_prec;
		}
		else {
			quda_inv_param.chrono_precision = theChromaToQudaPrecisionTypeMap::Instance()[ predictor.getChronoPrecision() ];
			QDPIO::cout << "Setting Chrono precision of " << quda_inv_param.chrono_precision << std::endl;
		}

		/// channel set done
		T Y_prime = zero;
		T Y = zero;
		// Y solve: M^\dagger Y = chi
		//        g_5 M g_5 Y = chi
		//     =>    M Y' = chi'  with chi' = gamma_5*chi
		Y_solve_timer.start();
		T g5chi = zero;
		g5chi[rb[1]]= Gamma(Nd*Nd-1)*chi;
		if( invParam.asymmetricP == true ) {
			res1 = qudaInvert(*clov,
					*invclov,
					g5chi,
					Y_prime);
			Y[rb[1]] = Gamma(Nd*Nd -1)*Y_prime;
		}
		else {
			T tmp = zero;
			invclov->apply(tmp,g5chi,MINUS,1);

			res1 = qudaInvert(*clov,
						*invclov,
						tmp,
						Y_prime);

#ifdef QUDA_DEBUG
			{
                           char Y_prime_norm[256];
                           char Y_prime_norm_full[256];
                           std::sprintf(Y_prime_norm, "%.*e", DECIMAL_DIG, toDouble(norm2(Y_prime, A->subset())));
                           std::sprintf(Y_prime_norm_full, "%.*e", DECIMAL_DIG, toDouble(norm2(Y_prime)));
                           QDPIO::cout << "Y solution: norm2(subset) = " << Y_prime_norm << " norm(full) = " << Y_prime_norm_full << std::endl;
                        }
#endif

			tmp[rb[1]] = Gamma(Nd*Nd-1)*Y_prime;
			clov->apply(Y,tmp,MINUS,1);

		}

		bool solution_good = true;

		// Check solution
		{
			T r;
			r[A->subset()]=chi;
			T tmp;
			(*A)(tmp, Y, MINUS);
			r[A->subset()] -= tmp;

			res1.resid = sqrt(norm2(r, A->subset()));
			if ( toBool( res1.resid/norm2chi > invParam.RsdToleranceFactor * invParam.RsdTarget ) ) {
				solution_good = false;
			}
		}

		if ( solution_good ) {
			if( res1.n_count >= threshold_counts ) {
				QDPIO::cout << solver_string << "Iteration Threshold Exceeded:Y Solver iters = " << res1.n_count << " Threshold=" << threshold_counts << std::endl;
				QDPIO::cout << solver_string << "Refreshing Subspace" << std::endl;

				Y_refresh_timer.start();
				// refresh the subspace
				// Setup the number of subspace Iterations
				for(int j=0; j < ip.mg_levels-1; j++) {
					(subspace_pointers->mg_param).setup_maxiter_refresh[j] = ip.maxIterSubspaceRefresh[j];
				}
				updateMultigridQuda(subspace_pointers->preconditioner, &(subspace_pointers->mg_param));
				for(int j=0; j < ip.mg_levels-1; j++) {
					(subspace_pointers->mg_param).setup_maxiter_refresh[j] = 0;
				}
				Y_refresh_timer.stop();
				QDPIO::cout << solver_string << "Y Subspace Refresh Time = " << Y_refresh_timer.getTimeInSeconds() << " secs\n";
			}
		}
		else {
			QDPIO::cout << solver_string << "Y-Solve failed (seq: "<<seqno<<"). Blowing away and reiniting subspace" << std::endl;
			StopWatch reinit_timer; reinit_timer.reset();
			reinit_timer.start();
			// BLow away subspace, re-set it up and then re-solve
			QUDAMGUtils::delete_subspace(invParam.SaveSubspaceID);

			// Recreate the subspace
			bool saved_value = ip.check_multigrid_setup;
                        ip.check_multigrid_setup = true;
                        subspace_pointers = QUDAMGUtils::create_subspace<T>(invParam);
                        ip.check_multigrid_setup = saved_value;


			// Make subspace XML snippets
			XMLBufferWriter file_xml;
			push(file_xml, "FileXML");
			pop(file_xml);

			int foo = 5;
			XMLBufferWriter record_xml;
			push(record_xml, "RecordXML");
			write(record_xml, "foo", foo);
			pop(record_xml);


			// Create named object entry.
			TheNamedObjMap::Instance().create< QUDAMGUtils::MGSubspacePointers* >(invParam.SaveSubspaceID);
			TheNamedObjMap::Instance().get(invParam.SaveSubspaceID).setFileXML(file_xml);
			TheNamedObjMap::Instance().get(invParam.SaveSubspaceID).setRecordXML(record_xml);

			// Assign the pointer into the named object
			TheNamedObjMap::Instance().getData< QUDAMGUtils::MGSubspacePointers* >(invParam.SaveSubspaceID) = subspace_pointers;
			quda_inv_param.preconditioner = subspace_pointers->preconditioner;
			reinit_timer.stop();
			QDPIO::cout << solver_string << "Subspace Reinit Time: " << reinit_timer.getTimeInSeconds() << " sec."  << std::endl;

			// Re-solve
			// This is a re-solve. So use_resident=false means used my initial guess
			// (do not repredict)
			quda_inv_param.chrono_use_resident = false;

			// The last solve, stored a chrono vector. We will overwrite this
			// thanks to the setting below
			quda_inv_param.chrono_replace_last = true;

			QDPIO::cout << solver_string << "Re-Solving for Y (zero guess)" << std::endl;
			SystemSolverResults_t res_tmp;
			Y_prime = zero;

			if( invParam.asymmetricP == true ) {
				res_tmp = qudaInvert(*clov,
						*invclov,
						g5chi,
						Y_prime);
				Y[rb[1]] = Gamma(Nd*Nd -1)*Y_prime;
			}
			else {
				T tmp = zero;
				invclov->apply(tmp,g5chi,MINUS,1);

				res_tmp = qudaInvert(*clov,
							*invclov,
							tmp,
							Y_prime);

#ifdef QUDA_DEBUG
				{
                           char Y_prime_norm[256];
                           char Y_prime_norm_full[256];
                           std::sprintf(Y_prime_norm, "%.*e", DECIMAL_DIG, toDouble(norm2(Y_prime, A->subset())));
                           std::sprintf(Y_prime_norm_full, "%.*e", DECIMAL_DIG, toDouble(norm2(Y_prime)));
                           QDPIO::cout << "Y solution: norm2(subset) = " << Y_prime_norm << " norm(full) = " << Y_prime_norm_full << std::endl;
                       	        }
#endif
				tmp[rb[1]] = Gamma(Nd*Nd-1)*Y_prime;
				clov->apply(Y,tmp,MINUS,1);
			}

			// Check solution
			{
				T r;
				r[A->subset()]=chi;
				T tmp;
				(*A)(tmp, Y, MINUS);
				r[A->subset()] -= tmp;

				res_tmp.resid = sqrt(norm2(r, A->subset()));
				if ( toBool( res_tmp.resid/norm2chi > invParam.RsdToleranceFactor * invParam.RsdTarget ) ) {
					// If we fail on the resolve then barf
					QDPIO::cout << solver_string << "Re Solve for Y Failed (seq: " << seqno << " )  Rsd = " << res_tmp.resid/norm2chi << " RsdTarget = " << invParam.RsdTarget << std::endl;

					dumpYSolver(g5chi,Y_prime);

					QDPIO::cout << solver_string << "Throwing Exception! This will ABORT" << std::endl;

					MGSolverException convergence_fail(invParam.CloverParams.Mass,
							invParam.SaveSubspaceID,
							res_tmp.n_count,
							Real(res_tmp.resid/norm2chi),
							invParam.RsdTarget*invParam.RsdToleranceFactor);
					throw convergence_fail;
				}
			}

			// At this point solution should be good again and subspace should be reinited
			res1.n_count += res_tmp.n_count; // Add resolve iterations
			res1.resid  = res_tmp.resid; // Copy new residuum.

		}
		Y_solve_timer.stop();

		// At this point we should have a good solution.
		// After the good solve, solution will be added to the right channel
		// by QUDA
		// Some diagnostics would be nice


		// Now select QUDA Chrono Index here
		quda_inv_param.chrono_max_dim = predictor.getMaxChrono();
		quda_inv_param.chrono_index = X_index;
		quda_inv_param.chrono_make_resident = true;
		quda_inv_param.chrono_use_resident = true;
		quda_inv_param.chrono_replace_last = false;

		// Reset Target Residuum for X solve
		quda_inv_param.tol = toDouble(invParam.RsdTarget);
		X_solve_timer.start();
		//psi[A->subset()]=zero;
		psi = zero;
		// Solve for psi
		res2 = qudaInvert(*clov,
				*invclov,
				Y,
				psi);
#ifdef QUDA_DEBUG
                   {
                           char X_prime_norm[256];
                           char X_prime_norm_full[256];
                           std::sprintf(X_prime_norm, "%.*e", DECIMAL_DIG, toDouble(norm2(psi, A->subset())));
                           std::sprintf(X_prime_norm_full, "%.*e", DECIMAL_DIG, toDouble(norm2(psi)));
                           QDPIO::cout << "X solution: norm2(subset) = " << X_prime_norm << " norm(full) = " << X_prime_norm_full << std::endl;
                        }
#endif


		solution_good = true;
		// Check solution
		{   
			T r=zero;
			r[A->subset()]=Y;
			T tmp=zero;
			// Checkin MX = Y solve
		        (*A)(tmp, psi, PLUS);
			r[ A->subset() ] -= tmp;
			Double resid_MXY = sqrt(norm2(r,A->subset()));
			Double normY = sqrt(norm2(Y,A->subset()));
			QDPIO::cout << "X solve: || Y - MX || / || Y || = " << resid_MXY/normY << std::endl;
			r[A->subset()]=chi;
			(*MdagM)(tmp, psi, PLUS);
			r[A->subset()] -= tmp;

			res2.resid = sqrt(norm2(r, A->subset()));
			if ( toBool( res2.resid/norm2chi > invParam.RsdToleranceFactor * invParam.RsdTarget ) ) {
				solution_good = false;
			}
		}


		if( solution_good )  {
			if( res2.n_count >= threshold_counts ) {
				QDPIO::cout << solver_string <<"Threshold Reached: X Solver iters = " << res2.n_count << " Threshold=" << threshold_counts << std::endl;
				QDPIO::cout << solver_string << "Refreshing Subspace" << std::endl;

				X_refresh_timer.start();
				// refresh the subspace
				// Regenerate space. Destroy and recreate
				// Setup the number of subspace Iterations
				for(int j=0; j < ip.mg_levels-1; j++) {
					(subspace_pointers->mg_param).setup_maxiter_refresh[j] = ip.maxIterSubspaceRefresh[j];
				}
				updateMultigridQuda(subspace_pointers->preconditioner, &(subspace_pointers->mg_param));
				for(int j=0; j < ip.mg_levels-1; j++) {
					(subspace_pointers->mg_param).setup_maxiter_refresh[j] = 0;
				}
				X_refresh_timer.stop();

				QDPIO::cout << solver_string << "Subspace Refresh Time = " << X_refresh_timer.getTimeInSeconds() << " secs\n";
			}
		}
		else {

			QDPIO::cout << solver_string << "X-Solve failed (seq: "<<seqno<<")  Blowing away and reiniting subspace" << std::endl;
			StopWatch reinit_timer; reinit_timer.reset();
			reinit_timer.start();

			// Delete the saved subspace completely
			QUDAMGUtils::delete_subspace(invParam.SaveSubspaceID);

			// Recreate the subspace
		        bool saved_value = ip.check_multigrid_setup;
                        ip.check_multigrid_setup = true;
                        subspace_pointers = QUDAMGUtils::create_subspace<T>(invParam);
                        ip.check_multigrid_setup = saved_value;
	
			// Make subspace XML snippets
			XMLBufferWriter file_xml;
			push(file_xml, "FileXML");
			pop(file_xml);

			int foo = 5;
			XMLBufferWriter record_xml;
			push(record_xml, "RecordXML");
			write(record_xml, "foo", foo);
			pop(record_xml);


			// Create named object entry.
			TheNamedObjMap::Instance().create< QUDAMGUtils::MGSubspacePointers* >(invParam.SaveSubspaceID);
			TheNamedObjMap::Instance().get(invParam.SaveSubspaceID).setFileXML(file_xml);
			TheNamedObjMap::Instance().get(invParam.SaveSubspaceID).setRecordXML(record_xml);

			// Assign the pointer into the named object
			TheNamedObjMap::Instance().getData< QUDAMGUtils::MGSubspacePointers* >(invParam.SaveSubspaceID) = subspace_pointers;
			quda_inv_param.preconditioner = subspace_pointers->preconditioner;	
			reinit_timer.stop();
			QDPIO::cout << solver_string << "Subspace Reinit Time: " << reinit_timer.getTimeInSeconds() << " sec."  << std::endl;

			// Re-solve
			// This is a re-solve. So use_resident=false means used my initial guess
			// (do not repredict)
			quda_inv_param.chrono_use_resident = false;

			// The last solve, stored a chrono vector. We will overwrite this
			// thanks to the setting below
			quda_inv_param.chrono_replace_last = true;

			QDPIO::cout << solver_string << "Re-Solving for X (zero guess)" << std::endl;

			SystemSolverResults_t res_tmp;
			
		 //	psi[rb[1]] = zero;
			psi = zero;
			res_tmp = qudaInvert(*clov,
					*invclov,
					Y,
					psi);

#ifdef QUDA_DEBUG
                   {
                           char X_prime_norm[256];
                           char X_prime_norm_full[256];
                           std::sprintf(X_prime_norm, "%.*e", DECIMAL_DIG, toDouble(norm2(psi, A->subset())));
                           std::sprintf(X_prime_norm_full, "%.*e", DECIMAL_DIG, toDouble(norm2(psi)));
                           QDPIO::cout << "X solution: norm2(subset) = " << X_prime_norm << " norm(full) = " << X_prime_norm_full << std::endl;
                        }
#endif
			// Check solution
			{
				T r=zero;
             			r[A->subset()]=Y;
                        	T tmp=zero;
                        	// Checkin MX = Y solve
                        	(*A)(tmp, psi, PLUS);
                        	r[ A->subset() ] -= tmp;
                       		Double resid_MXY = sqrt(norm2(r,A->subset()));
                        	Double normY = sqrt(norm2(Y,A->subset()));
                        	QDPIO::cout << "X re-solve: || Y - MX || / || Y || = " << resid_MXY/normY << std::endl;
				r[A->subset()]=chi;
				(*MdagM)(tmp, psi, PLUS);
				r[A->subset()] -= tmp;

				res_tmp.resid = sqrt(norm2(r, A->subset()));
				if ( toBool( res_tmp.resid/norm2chi > invParam.RsdToleranceFactor * invParam.RsdTarget ) ) {
					QDPIO::cout << solver_string << "Re Solve for X Failed (seq: " << seqno << " ) Rsd = " << res_tmp.resid/norm2chi << " RsdTarget = " << invParam.RsdTarget << std::endl;

					QDPIO::cout << "Dumping state (solve seqno : " << seqno << " ) " << std::endl;
					dumpXSolver(chi,Y,psi);


					QDPIO::cout << solver_string << "Throwing Exception! This will ABORT" << std::endl;
					MGSolverException convergence_fail(invParam.CloverParams.Mass,
							invParam.SaveSubspaceID,
							res_tmp.n_count,
							Real(res_tmp.resid/norm2chi),
							invParam.RsdTarget*invParam.RsdToleranceFactor);
					throw convergence_fail;
				}
			}
			// At this point the solution is good
			res2.n_count += res_tmp.n_count;
			res2.resid = res_tmp.resid;

		}
		X_solve_timer.stop();
		swatch.stop();
		double time = swatch.getTimeInSeconds();



		// Stats and done
		res.n_count = res1.n_count + res2.n_count;
		res.resid = res2.resid;

		Double rel_resid = res.resid/norm2chi;

		QDPIO::cout <<  solver_string  << " seq: " << (seqno++) << " iterations: " << res1.n_count << " + "
				<<  res2.n_count << " = " << res.n_count
				<<  " Rsd = " << res.resid << " Relative Rsd = " << rel_resid << std::endl;

		QDPIO::cout << "Y_solve: " << Y_solve_timer.getTimeInSeconds() << " (s) "
				<< "X_solve: " << X_solve_timer.getTimeInSeconds() << " (s) "
				<< "Total time: " << time <<  "(s)" << std::endl;

		quda_inv_param.use_init_guess = old_guess_policy;

		// Turn off chrono. Next solve can turn it on again
		quda_inv_param.chrono_make_resident = false;
		quda_inv_param.chrono_use_resident = false;
		quda_inv_param.chrono_replace_last = false;


		return res;
	}


	SystemSolverResults_t operator() (T& psi, const T& chi, Chroma::AbsChronologicalPredictor4D<T>& predictor ) const
	{
		SystemSolverResults_t res;

		// Try using QUDA predictor
		try {
			Chroma::QUDA4DChronoPredictor& quda_pred =
					dynamic_cast<Chroma::QUDA4DChronoPredictor&>(predictor);

			res = (*this)(psi,chi,quda_pred);
			return res;
		}
		catch(MGSolverException &e) {
			throw;
		}
		catch(...) {
			QDPIO::cout << "Failed to cast predictor to QUDA predictor"
					<< std::endl;
		}

		// QUDA Predictor failed -- Try abs 2 step
		try {
			Chroma::AbsTwoStepChronologicalPredictor4D<T>& two_step_pred =
					dynamic_cast< Chroma::AbsTwoStepChronologicalPredictor4D<T>&>(predictor);

			res = (*this)(psi,chi,two_step_pred);
			return res;
		}
		catch(MGSolverException &e) {
			throw;
		}
		catch(...) {
			QDPIO::cout << "Failed to cast predictor to QUDA or Two Step  predictor"
					<< std::endl;
			QDP_abort(1);
		}
	}


private:
	// Hide default constructor
	MdagMSysSolverQUDAMULTIGRIDClover() {}

#if 1
	Q links_orig;
#endif

	U GFixMat;
	QudaPrecision_s cpu_prec;
	QudaPrecision_s gpu_prec;
	QudaPrecision_s gpu_half_prec;

        Handle< LinearOperator<T> > A;
        Handle< FermState<T,Q,Q> > gstate;
	mutable SysSolverQUDAMULTIGRIDCloverParams invParam;
	QudaGaugeParam q_gauge_param;
	mutable QudaInvertParam quda_inv_param;
	mutable QUDAMGUtils::MGSubspacePointers* subspace_pointers;


	Handle< CloverTermT<T, U> > clov;
	Handle< CloverTermT<T, U> > invclov;

	SystemSolverResults_t qudaInvert(const CloverTermT<T, U>& clover,
			const CloverTermT<T, U>& inv_clov,
			const T& chi_s,
			T& psi_s
	)const;

	std::string solver_string;
	int threshold_counts;

        void dumpYSolver(const LatticeFermion& chi,
			 const LatticeFermion& Y) const;

        void dumpXSolver(const LatticeFermion& chi,
			 const LatticeFermion& Y,
			 const LatticeFermion& X) const;

	static unsigned long seqno;

};

} // End namespace

#endif // BUILD_QUDA
#endif 


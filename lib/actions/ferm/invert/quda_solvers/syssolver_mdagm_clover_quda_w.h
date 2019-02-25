// -*- C++ -*-
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by BiCGStab
 */

#ifndef __syssolver_mdagm_quda_clover_w_h__
#define __syssolver_mdagm_quda_clover_w_h__

#include "chroma_config.h"


#ifdef BUILD_QUDA

#include "handle.h"
#include "state.h"
#include "syssolver.h"
#include "linearop.h"
#include "lmdagm.h"
#include "actions/ferm/fermbcs/simple_fermbc.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/invert/quda_solvers/syssolver_quda_clover_params.h"
#include "actions/ferm/linop/clover_term_w.h"
#include "meas/gfix/temporal_gauge.h"
#include "io/aniso_io.h"
#include <string>
#include "update/molecdyn/predictor/null_predictor.h"
#include "update/molecdyn/predictor/quda_predictor.h"
#include "util/gauge/reunit.h"

#include <quda.h>
#ifdef QDP_IS_QDPJIT
#include "actions/ferm/invert/quda_solvers/qdpjit_memory_wrapper.h"
#endif

//#undef BUILD_QUDA_DEVIFACE_GAUGE
//#undef BUILD_QUDA_DEVIFACE_SPINOR
//#undef BUILD_QUDA_DEVIFACE_CLOVER

//#include <util_quda.h>

namespace Chroma
{

//! Richardson system solver namespace
namespace MdagMSysSolverQUDACloverEnv
{
//! Register the syssolver
bool registerAll();
}



//! Solve a Clover Fermion System using the QUDA inverter
/*! \ingroup invert
 *** WARNING THIS SOLVER WORKS FOR Clover FERMIONS ONLY ***
 */

class MdagMSysSolverQUDAClover : public MdagMSystemSolver<LatticeFermion>
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
	//! Constructor
	/*!
	 * \param M_        Linear operator ( Read )
	 * \param invParam  inverter parameters ( Read )
	 */
	MdagMSysSolverQUDAClover(Handle< LinearOperator<T> > A_,
			Handle< FermState<T,Q,Q> > state_,
			const SysSolverQUDACloverParams& invParam_) :
				A(A_), state(state_), invParam(invParam_), clov(new CloverTermT<T, U>()), invclov(new CloverTermT<T, U>())
	{
		QDPIO::cout << "MdagMSysSolverQUDAClover:" << std::endl;

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
		//QDPIO::cout << "MDAGM Using QDP-JIT gauge order" << std::endl;
		q_gauge_param.location    = QUDA_CUDA_FIELD_LOCATION;
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

		// Default for no preconditioning.
		// Data read from the innner preconditioner struct
		// can overwrite this later
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

		// Default: may be overwritten later
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
			QDPIO::cout << "Fixing Temporal Gauge" << std::endl;
			temporalGauge(links_single, GFixMat, Nd-1);
			for(int mu=0; mu < Nd; mu++){
				links_single[mu] = GFixMat*(state_->getLinks())[mu]*adj(shift(GFixMat, FORWARD, mu));
			}
			q_gauge_param.gauge_fix = QUDA_GAUGE_FIXED_YES;
		}
		else {
			// No GaugeFix
			q_gauge_param.gauge_fix = QUDA_GAUGE_FIXED_NO;  // No Gfix yet
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

		// Invert type:
		switch( invParam.solverType ) {
		case CG:
			quda_inv_param.inv_type = QUDA_CG_INVERTER;
			solver_string = "CG";
			break;
		case BICGSTAB:
			quda_inv_param.inv_type = QUDA_BICGSTAB_INVERTER;
			solver_string = "BICGSTAB";
			break;
		case GCR:
			quda_inv_param.inv_type = QUDA_GCR_INVERTER;
			solver_string = "GCR";
			break;
		default:
			quda_inv_param.inv_type = QUDA_CG_INVERTER;
			solver_string = "CG";
			break;
		}

		// Turn off true residuum calculation for MdagM solves
		// Since Chroma will check on this.
		quda_inv_param.compute_true_res = 0;

		// Mass

		// Fiendish idea from Ron. Set the kappa=1/2 and use
		// unmodified clover term, and ask for Kappa normalization
		// This should give us A - (1/2)D as the unpreconditioned operator
		// and probabl 1 - {1/4} A^{-1} D A^{-1} D as the preconditioned
		// op. Apart from the A_oo stuff on the antisymmetric we have
		// nothing to do...
		quda_inv_param.kappa = 0.5;

		// Dummy parameter to pass check_params
		// FIXME: If we ever use QUDA to compute our clover term we have to sort out anisotropy
		// RIght now this won't do anything since we pass in the clover term
		quda_inv_param.clover_coeff = 1.0;
		quda_inv_param.tol = toDouble(invParam.RsdTarget);
		quda_inv_param.maxiter = invParam.MaxIter;
		quda_inv_param.reliable_delta = toDouble(invParam.Delta);
		quda_inv_param.pipeline = invParam.Pipeline;

		// Solution type
		quda_inv_param.solution_type = QUDA_MATPCDAG_MATPC_SOLUTION;

		// Solve type
		switch( invParam.solverType ) {
		case CG:
			quda_inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;
			break;
		case BICGSTAB:
			quda_inv_param.solve_type = QUDA_DIRECT_PC_SOLVE;
			break;

		case GCR:
			quda_inv_param.solve_type = QUDA_DIRECT_PC_SOLVE;
			break;
		case CA_GCR:
			quda_inv_param.solve_type = QUDA_DIRECT_PC_SOLVE;
			break;
		case MR:
			quda_inv_param.solve_type = QUDA_DIRECT_PC_SOLVE;
			break;

		default:
			quda_inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;

			break;
		}

		
		if( invParam.asymmetricP ) {
		  QDPIO::cout << "Working with Asymmetric LinOp" << std::endl;
		  quda_inv_param.matpc_type = QUDA_MATPC_ODD_ODD_ASYMMETRIC;
		}
		else { 
		  QDPIO::cout << "Working with Symmetric LinOp" << std::endl;
		  quda_inv_param.matpc_type = QUDA_MATPC_ODD_ODD;
                }

		quda_inv_param.dagger = QUDA_DAG_NO;
		quda_inv_param.mass_normalization = QUDA_KAPPA_NORMALIZATION;

		quda_inv_param.cpu_prec = cpu_prec;
		quda_inv_param.cuda_prec = gpu_prec;
		quda_inv_param.cuda_prec_sloppy = gpu_half_prec;

		// Default for no preconditioining
		// Data read from inner parameters can override this later.
		quda_inv_param.cuda_prec_precondition = gpu_half_prec;

		quda_inv_param.preserve_source = QUDA_PRESERVE_SOURCE_NO;
		quda_inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;

#ifndef BUILD_QUDA_DEVIFACE_SPINOR
		quda_inv_param.dirac_order = QUDA_DIRAC_ORDER;
		quda_inv_param.input_location = QUDA_CPU_FIELD_LOCATION;
		quda_inv_param.output_location = QUDA_CPU_FIELD_LOCATION;
#else
		//QDPIO::cout << "MDAGM Using QDP-JIT spinor order" << std::endl;
		quda_inv_param.dirac_order    = QUDA_QDPJIT_DIRAC_ORDER;
		quda_inv_param.input_location = QUDA_CUDA_FIELD_LOCATION;
		quda_inv_param.output_location = QUDA_CUDA_FIELD_LOCATION;
#endif

		// Clover precision and order
		quda_inv_param.clover_cpu_prec = cpu_prec;
		quda_inv_param.clover_cuda_prec = gpu_prec;
		quda_inv_param.clover_cuda_prec_sloppy = gpu_half_prec;

		// Default for no preconditioning
		// Data read from iner parameters can overwrite this later.
		quda_inv_param.clover_cuda_prec_precondition = gpu_half_prec;

#ifndef BUILD_QUDA_DEVIFACE_CLOVER
		quda_inv_param.clover_order = QUDA_PACKED_CLOVER_ORDER;
#else
		//QDPIO::cout << "MDAGM Clover CUDA location\n";
		quda_inv_param.clover_location = QUDA_CUDA_FIELD_LOCATION;
		quda_inv_param.clover_order = QUDA_QDPJIT_CLOVER_ORDER;
#endif


		// Autotuning
		if( invParam.tuneDslashP ) {
			QDPIO::cout << "Enabling Dslash Autotuning" << std::endl;

			quda_inv_param.tune = QUDA_TUNE_YES;
		}
		else {
			QDPIO::cout << "Disabling Dslash Autotuning" << std::endl;

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
		quda_inv_param.sp_pad = 0;
		quda_inv_param.cl_pad = 0;

		if( invParam.innerParamsP ) {
			QDPIO::cout << "Setting inner solver params" << std::endl;
			// Dereference handle
			GCRInnerSolverParams ip = *(invParam.innerParams);

			switch( ip.precPrecondition ) {
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

			switch( ip.reconstructPrecondition ) {
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

			quda_inv_param.tol_precondition = toDouble(ip.tolPrecondition);
			quda_inv_param.maxiter_precondition = ip.maxIterPrecondition;
			quda_inv_param.gcrNkrylov = ip.gcrNkrylov;
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
			quda_inv_param.precondition_cycle = ip.preconditionCycle;

			if( ip.verboseInner ) {
				quda_inv_param.verbosity_precondition = QUDA_VERBOSE;
			}
			else {
				quda_inv_param.verbosity_precondition = QUDA_SILENT;
			}

			switch( ip.invTypePrecondition ) {
			case CG:
				quda_inv_param.inv_type_precondition = QUDA_CG_INVERTER;
				break;
			case BICGSTAB:
				quda_inv_param.inv_type_precondition = QUDA_BICGSTAB_INVERTER;
				break;
			case CA_GCR:
				quda_inv_param.inv_type_precondition= QUDA_CA_GCR_INVERTER;
				break;
			case MR:
				quda_inv_param.inv_type_precondition= QUDA_MR_INVERTER;
				break;

			default:
				quda_inv_param.inv_type_precondition = QUDA_MR_INVERTER;
				break;
			}
		}
		else {
			QDPIO::cout << "Setting Precondition stuff to defaults for not using" << std::endl;
			quda_inv_param.inv_type_precondition= QUDA_INVALID_INVERTER;
			quda_inv_param.tol_precondition = 1.0e-1;
			quda_inv_param.maxiter_precondition = 1000;
			quda_inv_param.verbosity_precondition = QUDA_SILENT;
			quda_inv_param.gcrNkrylov = 1;

		}

		if( invParam.verboseP ) {
			quda_inv_param.verbosity = QUDA_VERBOSE;
		}
		else {
			quda_inv_param.verbosity = QUDA_SUMMARIZE;
		}

		// Set up the links
		void* gauge[4];

#ifndef BUILD_QUDA_DEVIFACE_GAUGE
		for(int mu=0; mu < Nd; mu++) {
			gauge[mu] = (void *)&(links_single[mu].elem(all.start()).elem().elem(0,0).real());
		}
#else
		GetMemoryPtrGauge(gauge,links_single);
		//gauge[mu] = GetMemoryPtr( links_single[mu].getId() );
			//std::cout << "MDAGM CUDA gauge[" << mu << "] in = " << gauge[mu] << "\n";
#endif


		loadGaugeQuda((void *)gauge, &q_gauge_param);

		// Setup Clover Term
		QDPIO::cout << "Creating CloverTerm" << std::endl;
		clov->create(fstate, invParam_.CloverParams);
		// Don't recompute, just copy
		invclov->create(fstate, invParam_.CloverParams);

		QDPIO::cout << "Inverting CloverTerm" << std::endl;
		invclov->choles(0);
		invclov->choles(1);

#ifndef BUILD_QUDA_DEVIFACE_CLOVER
		multi1d<QUDAPackedClovSite<REALT> > packed_clov;

		packed_clov.resize(all.siteTable().size());

		clov->packForQUDA(packed_clov, 0);
		clov->packForQUDA(packed_clov, 1);

		// Always need inverse
		multi1d<QUDAPackedClovSite<REALT> > packed_invclov(all.siteTable().size());
		invclov->packForQUDA(packed_invclov, 0);
		invclov->packForQUDA(packed_invclov, 1);

		loadCloverQuda(&(packed_clov[0]), &(packed_invclov[0]),&quda_inv_param);
#else
		void *clover[2];
		void *cloverInv[2];

		GetMemoryPtrClover(clov->getOffId(),clov->getDiaId(),invclov->getOffId(),invclov->getDiaId());

		
		loadCloverQuda( (void*)(clover) , (void*)(cloverInv) ,&quda_inv_param);
#endif
		if ( invParam.CloverParams.twisted_m_usedP == true ) {
			quda_inv_param.dslash_type = QUDA_CLOVER_HASENBUSCH_TWIST_DSLASH;
			if ( invParam.asymmetricP ) {
				// - comes from add/subtract convention
				quda_inv_param.m5 = -toDouble(invParam.CloverParams.twisted_m);
			}
			else {
				// - comes from add/subtract convention
				// 2 comes from clover normalization
				quda_inv_param.m5 = -2.0*toDouble(invParam.CloverParams.twisted_m);
			}
		}
	}


	//! Destructor is automatic
	~MdagMSysSolverQUDAClover()
	{
		QDPIO::cout << "Destructing" << std::endl;
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
	SystemSolverResults_t operator() (T& psi, const T& chi ) const
	{
		SystemSolverResults_t res;
		Null4DChronoPredictor not_predicting;
		(*this)(psi,chi, not_predicting);
		START_CODE();
		END_CODE();
		return res;
	}


	SystemSolverResults_t operator() (T& psi, const T& chi, Chroma::QUDA4DChronoPredictor& predictor ) const
	{
		SystemSolverResults_t res;
		SystemSolverResults_t res1;
		SystemSolverResults_t res2;

		START_CODE();
		// This is a solve with initial guess. So reset the policy  (quda_inv_param is mutable)
		QudaUseInitGuess old_guess_policy = quda_inv_param.use_init_guess;
		quda_inv_param.use_init_guess = QUDA_USE_INIT_GUESS_YES;


		StopWatch swatch;
		swatch.start();

		// Create MdagM op
		Handle< LinearOperator<T> > MdagM( new MdagMLinOp<T>(A) );

		StopWatch Y_solve_timer;      Y_solve_timer.reset();
		StopWatch X_solve_timer;      X_solve_timer.reset();


		// TWO STEP SOLVE
		QDPIO::cout << "Two Step Solve using QUDA predictor: (X_index,Y_index) = ( " << predictor.getXIndex() << " , " << predictor.getYIndex() << " ) \n";
		quda_inv_param.chrono_max_dim = predictor.getMaxChrono();
		quda_inv_param.chrono_index = predictor.getYIndex();
		quda_inv_param.chrono_make_resident = true;
		quda_inv_param.chrono_use_resident = true;
		quda_inv_param.chrono_replace_last = false;

		T Y;
		Y[ A->subset() ] = psi; // Y is initial guess

		// Set up prediction for Y in QUDA here.

		Y_solve_timer.start();
		// Now want to solve with M^\dagger on this.
		quda_inv_param.solution_type = QUDA_MATPC_SOLUTION;
		quda_inv_param.dagger = QUDA_DAG_YES;
		res1 = qudaInvert(*clov,
				*invclov,
				chi,
				Y);
		Y_solve_timer.stop();


		// Set up prediction for X in QUDA here
		quda_inv_param.chrono_index = predictor.getXIndex();
		quda_inv_param.chrono_make_resident = true;
		quda_inv_param.chrono_use_resident = true;
		quda_inv_param.chrono_replace_last = false;


		// Step 2: Solve M X = Y
		// Predict X

		X_solve_timer.start();
		quda_inv_param.dagger = QUDA_DAG_NO; // Solve without dagger
		res2 = qudaInvert(*clov,
				*invclov,
				Y,
				psi);
		X_solve_timer.stop();
		swatch.stop();
		double time = swatch.getTimeInSeconds();

		// reset init guess policy

		// Check solution
		{
			T r;
			r[A->subset()]=chi;
			T tmp;
			(*MdagM)(tmp, psi, PLUS);
			r[A->subset()] -= tmp;
			res.resid = sqrt(norm2(r, A->subset()));
		}

		Double rel_resid = res.resid/sqrt(norm2(chi,A->subset()));
		res.n_count = res1.n_count + res2.n_count;

		QDPIO::cout << "QUDA_"<< solver_string <<"_CLOVER_SOLVER: " << res.n_count << " iterations. Rsd = " << res.resid << " Relative Rsd = " << rel_resid << std::endl
				;

		QDPIO::cout << "QUDA_"<< solver_string <<"_CLOVER_SOLVER: Time: "
				<< "Y_solve: " << Y_solve_timer.getTimeInSeconds() << " (s) "
				<< "X_solve: " << X_solve_timer.getTimeInSeconds() << " (s) "
				<< "Total time: " << time <<  "(s)" << std::endl;

		quda_inv_param.use_init_guess = old_guess_policy;
		quda_inv_param.chrono_make_resident = false;
		quda_inv_param.chrono_use_resident = false;
		quda_inv_param.chrono_replace_last = false;

		// Check for failure
		if (  toBool( rel_resid >  invParam.RsdToleranceFactor*invParam.RsdTarget) ) {
			QDPIO::cout << "QUDA_" << solver_string << "_CLOVER_SOLVER: SOLVE Failed to converge" << std::endl;
			QDP_abort(1);
		}



		END_CODE();
		return res;

	}


	SystemSolverResults_t operator() (T& psi, const T& chi, Chroma::AbsTwoStepChronologicalPredictor4D<T>& predictor ) const
	{
		SystemSolverResults_t res;
		SystemSolverResults_t res1;
		SystemSolverResults_t res2;

		START_CODE();
		// This is a solve with initial guess. So reset the policy  (quda_inv_param is mutable)
		QudaUseInitGuess old_guess_policy = quda_inv_param.use_init_guess;
		quda_inv_param.use_init_guess = QUDA_USE_INIT_GUESS_YES;


		StopWatch swatch;
		swatch.start();

		// This is a solver whic will use an initial guess:

		// Create MdagM op
		Handle< LinearOperator<T> > MdagM( new MdagMLinOp<T>(A) );

		StopWatch X_prediction_timer; X_prediction_timer.reset();
		StopWatch Y_prediction_timer; Y_prediction_timer.reset();
		StopWatch Y_solve_timer;      Y_solve_timer.reset();
		StopWatch X_solve_timer;      X_solve_timer.reset();
		StopWatch Y_predictor_add_timer; Y_predictor_add_timer.reset();
		StopWatch X_predictor_add_timer; X_predictor_add_timer.reset();

		QDPIO::cout << "Two Step Solve" << std::endl;

		T Y;
		Y[ A->subset() ] = psi; // Y is initial guess

		// Predict Y
		Y_prediction_timer.start();
		predictor.predictY(Y,*A,chi);
		Y_prediction_timer.stop();

		Y_solve_timer.start();
		// Now want to solve with M^\dagger on this.
		quda_inv_param.solution_type = QUDA_MATPC_SOLUTION;
		quda_inv_param.dagger = QUDA_DAG_YES;
		res1 = qudaInvert(*clov,
				*invclov,
				chi,
				Y);
		Y_solve_timer.stop();

		Y_predictor_add_timer.start();
		predictor.newYVector(Y);
		Y_predictor_add_timer.stop();


		// Step 2: Solve M X = Y
		// Predict X
		X_prediction_timer.start();
		predictor.predictX(psi,*MdagM, chi);
		X_prediction_timer.stop();
		X_solve_timer.start();
		quda_inv_param.dagger = QUDA_DAG_NO; // Solve without dagger
		res2 = qudaInvert(*clov,
				*invclov,
				Y,
				psi);
		X_solve_timer.stop();
		X_predictor_add_timer.start();
		predictor.newXVector(psi);
		X_predictor_add_timer.stop();

		res.n_count = res1.n_count + res2.n_count;  // Two step solve so combine iteration count

		swatch.stop();
		double time = swatch.getTimeInSeconds();

		// reset init guess policy
		quda_inv_param.use_init_guess = old_guess_policy;

		// Check solution
		{
			T r;
			r[A->subset()]=chi;
			T tmp;
			(*MdagM)(tmp, psi, PLUS);
			r[A->subset()] -= tmp;
			res.resid = sqrt(norm2(r, A->subset()));
		}

		Double rel_resid = res.resid/sqrt(norm2(chi,A->subset()));

		QDPIO::cout << "QUDA_"<< solver_string <<"_CLOVER_SOLVER: " << res.n_count << " iterations. Rsd = " << res.resid << " Relative Rsd = " << rel_resid << std::endl
				;




		QDPIO::cout << "QUDA_"<< solver_string <<"_CLOVER_SOLVER: Time: Y_predict: " << Y_prediction_timer.getTimeInSeconds() << " (s) "
				<< "Y_solve: " << Y_solve_timer.getTimeInSeconds() << " (s) "
				<< "Y_register: " << Y_predictor_add_timer.getTimeInSeconds() << " (s) "
				<< "X_predict: " << X_prediction_timer.getTimeInSeconds() << " (s) "
				<< "X_solve: " << X_solve_timer.getTimeInSeconds() << " (s) "
				<< "X_register: " << X_predictor_add_timer.getTimeInSeconds() << " (s) "
				<< "Total time: " << time <<  "(s)" << std::endl;

		if (  toBool( rel_resid >  invParam.RsdToleranceFactor*invParam.RsdTarget) ) {
			QDPIO::cout << "QUDA_"<< solver_string <<"_CLOVER_SOLVER: Solver Failed to Converge! Aborting" << std::endl;
			QDP_abort(1);
		}


		END_CODE();
		return res;
	}

	SystemSolverResults_t operator() (T& psi, const T& chi, Chroma::AbsChronologicalPredictor4D<T>& predictor ) const
	{
		SystemSolverResults_t res;
		START_CODE();


		StopWatch swatch;

		// Determine if 2 step solve is needed
		bool two_step=false;
		switch(  invParam.solverType  ) {
		case BICGSTAB:
			two_step = true;
			break;
		case GCR:
			two_step = true;
			break;
		case CA_GCR:
			two_step = true;
			break;
		case MR:
			two_step = true;
			break;
		case CG:
			two_step = false;
			break;
		default:
			two_step = false;
			break;
		};

		// This is a solver whic will use an initial guess:


		if ( two_step ) {
			// User must supply a two step predictor. This can be
			// a general one or a specific one. Try to cast to the specific one first.
			try {
				QUDA4DChronoPredictor& quda_predictor =
						dynamic_cast<QUDA4DChronoPredictor&>(predictor);

				res = (*this)(psi,chi,quda_predictor);
				END_CODE();
				return res;

			}
			catch(std::bad_cast) {} // We have another to try

			// Try a general two step predictor cast
			try {
				AbsTwoStepChronologicalPredictor4D<T>& abs_two_step_predictor =
						dynamic_cast<AbsTwoStepChronologicalPredictor4D<T>&>(predictor);
				res = (*this)(psi,chi,abs_two_step_predictor);
				END_CODE();
				return res;
			}
			catch(std::bad_cast) {
				// Now we are stuck.
				QDPIO::cout << "QUDA_"<< solver_string <<"_CLOVER_SOLVER: Two Step Solver with two step incapable predictor" << std::endl;
				QDP_abort(1);
			}

		}


		StopWatch X_prediction_timer; X_prediction_timer.reset();
		StopWatch X_solve_timer;      X_solve_timer.reset();
		StopWatch X_predictor_add_timer; X_predictor_add_timer.reset();

		// This is a solve with initial guess. So reset the policy  (quda_inv_param is mutable)
		QudaUseInitGuess old_guess_policy = quda_inv_param.use_init_guess;
		quda_inv_param.use_init_guess = QUDA_USE_INIT_GUESS_YES;


		// Create MdagM op
		Handle< LinearOperator<T> > MdagM( new MdagMLinOp<T>(A) );

		// If we are here it is not a two step.
		swatch.start();

		// Single Step Solve
		QDPIO::cout << "Single Step Solve" << std::endl;
		X_prediction_timer.start();
		predictor(psi, (*MdagM), chi);
		X_prediction_timer.stop();

		X_solve_timer.start();
		res = qudaInvert(*clov,
				*invclov,
				chi,
				psi);
		X_solve_timer.stop();

		X_predictor_add_timer.start();
		predictor.newVector(psi);
		X_predictor_add_timer.stop();

		// reset init guess policy
		quda_inv_param.use_init_guess = old_guess_policy;

		// Check solution
		{
			T r;
			r[A->subset()]=chi;
			T tmp;
			(*MdagM)(tmp, psi, PLUS);
			r[A->subset()] -= tmp;
			res.resid = sqrt(norm2(r, A->subset()));
		}
		swatch.stop();
		double time=swatch.getTimeInSeconds();
		Double rel_resid = res.resid/sqrt(norm2(chi,A->subset()));

		QDPIO::cout << "QUDA_"<< solver_string <<"_CLOVER_SOLVER: " << res.n_count << " iterations. Rsd = " << res.resid << " Relative Rsd = " << rel_resid << std::endl
				;
		QDPIO::cout << "QUDA_"<< solver_string <<"_CLOVER_SOLVER: Time: X_predict: " << X_prediction_timer.getTimeInSeconds() << " (s) "
				<< "X_solve: " << X_solve_timer.getTimeInSeconds() << " (s) "
				<< "X_register: " << X_predictor_add_timer.getTimeInSeconds() << " (s) "
				<< "Total time: " << time <<  "(s)" << std::endl;

		if (  toBool( rel_resid >  invParam.RsdToleranceFactor*invParam.RsdTarget) ) {
			QDPIO::cout << "QUDA_"<< solver_string <<"_CLOVER_SOLVER: Solver Failed to Converge! Aborting" << std::endl;
			QDP_abort(1);
		}

		END_CODE();
		return res;

	}


private:
	// Hide default constructor
	MdagMSysSolverQUDAClover() {}

	Q links_orig;

	U GFixMat;
	QudaPrecision_s cpu_prec;
	QudaPrecision_s gpu_prec;
	QudaPrecision_s gpu_half_prec;

	Handle< LinearOperator<T> > A;
	Handle< FermState<T,Q,Q> > state;
	const SysSolverQUDACloverParams invParam;
	QudaGaugeParam q_gauge_param;
	mutable QudaInvertParam quda_inv_param;

	Handle< CloverTermT<T, U> > clov;
	Handle< CloverTermT<T, U> > invclov;

	SystemSolverResults_t qudaInvert(const CloverTermT<T, U>& clover,
			const CloverTermT<T, U>& inv_clov,
			const T& chi_s,
			T& psi_s
	)const ;

	std::string solver_string;
};


} // End namespace

#endif // BUILD_QUDA
#endif 


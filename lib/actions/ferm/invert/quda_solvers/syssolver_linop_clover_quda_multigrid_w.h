// -*- C++ -*-
/*! \file
 *  \QUDA MULTIGRID Clover solver.
 */

#ifndef __syssolver_linop_quda_multigrid_clover_h__
#define __syssolver_linop_quda_multigrid_clover_h__

#include "chroma_config.h"

#ifdef BUILD_QUDA
#include <quda.h>

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
#include "quda_mg_utils.h"
#include <string>
#include <sstream>
#include "util/gauge/reunit.h"
#ifdef QDP_IS_QDPJIT
#include "actions/ferm/invert/quda_solvers/qdpjit_memory_wrapper.h"
#endif

//#include <util_quda.h>

namespace Chroma
{

	//! Richardson system solver namespace
	namespace LinOpSysSolverQUDAMULTIGRIDCloverEnv
	{
		//! Register the syssolver
		bool registerAll();
	}

	//! Solve a Clover Fermion System using the QUDA inverter
	/*! \ingroup invert
	 *** WARNING THIS SOLVER WORKS FOR Clover FERMIONS ONLY ***
	 */

	class LinOpSysSolverQUDAMULTIGRIDClover : public LinOpSystemSolver<LatticeFermion>
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
		LinOpSysSolverQUDAMULTIGRIDClover(Handle< LinearOperator<T> > A_,
				Handle< FermState<T,Q,Q> > state_,
				const SysSolverQUDAMULTIGRIDCloverParams& invParam_) :
		A(A_), invParam(invParam_), clov(new CloverTermT<T, U>() ), invclov(new CloverTermT<T, U>())
		{
		  StopWatch init_swatch;
		  init_swatch.reset(); init_swatch.start();
		   // Set the solver string
		    {
		      std::ostringstream solver_string_stream;
		      solver_string_stream << "QUDA_MULTIGRID_CLOVER_LINOP_SOLVER( "
		          << invParam.SaveSubspaceID << " ): ";
		      solver_string = solver_string_stream.str();

		    }
		    QDPIO::cout << solver_string << "Initializing" << std::endl;

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

			// Default. This may be overwritten by inner params
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
		    	q_gauge_param.reconstruct_precondition = q_gauge_param.reconstruct_sloppy;
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


			quda_inv_param.kappa = 0.5;
			quda_inv_param.clover_coeff = 1.0; // Dummy, not used
			quda_inv_param.Ls = 1;

			quda_inv_param.tol = toDouble(invParam.RsdTarget);
			quda_inv_param.maxiter = invParam.MaxIter;
			quda_inv_param.reliable_delta = toDouble(invParam.Delta);
			quda_inv_param.pipeline = invParam.Pipeline;

			// Solution type
			//quda_inv_param.solution_type = QUDA_MATPC_SOLUTION;
			//Taken from invert test.
			quda_inv_param.solution_type = QUDA_MATPC_SOLUTION;
			quda_inv_param.solve_type = QUDA_DIRECT_PC_SOLVE;

			// Solve type
			/*switch( invParam.solverType ) {
			 case CG:
			 quda_inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;
			 break;
			 case BICGSTAB:
			 quda_inv_param.solve_type = QUDA_DIRECT_PC_SOLVE;
			 break;
			 case GCR:
			 quda_inv_param.solve_type = QUDA_DIRECT_PC_SOLVE;
			 break;

			 case MR:
			 quda_inv_param.solve_type = QUDA_DIRECT_PC_SOLVE;
			 break;

			 default:
			 quda_inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;

			 break;
			 }*/

			quda_inv_param.matpc_type = QUDA_MATPC_ODD_ODD;

			quda_inv_param.dagger = QUDA_DAG_NO;
			quda_inv_param.mass_normalization = QUDA_KAPPA_NORMALIZATION;

			quda_inv_param.cpu_prec = cpu_prec;
			quda_inv_param.cuda_prec = gpu_prec;
			quda_inv_param.cuda_prec_sloppy = gpu_half_prec;
			quda_inv_param.cuda_prec_precondition = gpu_half_prec;

			//
			//Done...
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
		        quda_inv_param.clover_cuda_prec_precondition  = gpu_half_prec;
			if( invParam.MULTIGRIDParamsP ) {
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
			}
			// Set up the links
			void* gauge[4];

#ifndef BUILD_QUDA_DEVIFACE_GAUGE
			for(int mu=0; mu < Nd; mu++) {
				gauge[mu] = (void *)&(links_single[mu].elem(all.start()).elem().elem(0,0).real());
			}
#else
				//gauge[mu] = GetMemoryPtr( links_single[mu].getId() );
				GetMemoryPtrGauge(gauge,links_single);
#endif

			loadGaugeQuda((void *)gauge, &q_gauge_param);

			MULTIGRIDSolverParams ip = *(invParam.MULTIGRIDParams);
			//
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

			QDPIO::cout<< solver_string << "Basic MULTIGRID params copied."<<std::endl;

			//      Setup the clover term...
			QDPIO::cout << solver_string << "Creating CloverTerm" << std::endl;
			clov->create(fstate, invParam_.CloverParams);

			// Don't recompute, just copy
			invclov->create(fstate, invParam_.CloverParams);

			QDPIO::cout << solver_string << "Inverting CloverTerm" << std::endl;
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


// merged from mdgam_clover_quda_multigrid, begin
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
// merged from mdgam_clover_quda_multigrid, end  

			init_swatch.stop();
			QDPIO::cout << solver_string << "init_time = "
			            << init_swatch.getTimeInSeconds() << " sec. " << std::endl;

		}

		//! Destructor is automatic
		~LinOpSysSolverQUDAMULTIGRIDClover()
		{
			QDPIO::cout << solver_string << "Destructing" << std::endl;
			quda_inv_param.preconditioner = nullptr;
			subspace_pointers = nullptr;
			freeGaugeQuda();
			freeCloverQuda();
//			destroyMultigridQuda(quda_inv_param.preconditioner);
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
			SystemSolverResults_t res;

			START_CODE();
			StopWatch swatch;
			swatch.start();

			psi = zero; // Zero initial guess
			//    T MdagChi;

			// This is a CGNE. So create new RHS
			//      (*A)(MdagChi, chi, MINUS);
			// Handle< LinearOperator<T> > MM(new MdagMLinOp<T>(A));
			if ( invParam.axialGaugeP ) {
				T g_chi,g_psi;

				// Gauge Fix source and initial guess
				g_chi[ rb[1] ] = GFixMat * chi;
				g_psi[ rb[1] ] = GFixMat * psi;
				res = qudaInvert(*clov,
						*invclov,
						g_chi,
						g_psi);
				psi[ rb[1]] = adj(GFixMat)*g_psi;

			}
			else {
				res = qudaInvert(*clov,
						*invclov,
						chi,
						psi);
			}

			swatch.stop();
			Double rel_resid;

			if( invParam.SolutionCheckP )  {
			  
			  {
			    T r;
			    r[A->subset()]=chi;
			    T tmp;
			    (*A)(tmp, psi, PLUS);
			    r[A->subset()] -= tmp;
			    res.resid = sqrt(norm2(r, A->subset()));
			  }

			  rel_resid = res.resid/sqrt(norm2(chi,A->subset()));

			  QDPIO::cout << solver_string  << res.n_count << " iterations. Rsd = " << res.resid << " Relative Rsd = " << rel_resid << std::endl;
			}
			else { 
			  // just believe the QUDA residuum.
			  // which is always a true reiduum
			  rel_resid = res.resid;
			}

			// Convergence Check/Blow Up
			if ( ! invParam.SilentFailP ) {
				if ( toBool( rel_resid > invParam.RsdToleranceFactor*invParam.RsdTarget) ) {
					QDPIO::cerr << solver_string << "ERROR:   Solver residuum is outside tolerance: QUDA resid="<< rel_resid << " Desired =" << invParam.RsdTarget << " Max Tolerated = " << invParam.RsdToleranceFactor*invParam.RsdTarget << std::endl;
					QDP_abort(1);
				}
			}

			END_CODE();
			return res;
		}

	private:
		// Hide default constructor
		LinOpSysSolverQUDAMULTIGRIDClover() {}

#if 1
		Q links_orig;
#endif

		U GFixMat;
		QudaPrecision_s cpu_prec;
		QudaPrecision_s gpu_prec;
		QudaPrecision_s gpu_half_prec;

		Handle< LinearOperator<T> > A;
		const SysSolverQUDAMULTIGRIDCloverParams invParam;
		QudaGaugeParam q_gauge_param;
		QudaInvertParam quda_inv_param;
		mutable QUDAMGUtils::MGSubspacePointers* subspace_pointers;

		Handle< CloverTermT<T, U> > clov;
		Handle< CloverTermT<T, U> > invclov;

		SystemSolverResults_t qudaInvert(const CloverTermT<T, U>& clover,
				const CloverTermT<T, U>& inv_clov,
				const T& chi_s,
				T& psi_s
		)const;

		std::string solver_string;
	};

} // End namespace

#endif // BUILD_QUDA
#endif 


// -*- C++ -*-
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by BiCGStab
 */

#ifndef __syssolver_linop_quda_clover_h__
#define __syssolver_linop_quda_clover_h__

#include "chroma_config.h"

#ifdef BUILD_QUDA
#include <quda.h>

#include "handle.h"
#include "state.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/fermbcs/simple_fermbc.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/invert/quda_solvers/syssolver_quda_clover_params.h"
#include "actions/ferm/linop/clover_term_w.h"
#include "meas/gfix/temporal_gauge.h"
#include "io/aniso_io.h"
#include <string>
#include <vector>

#include "util/gauge/reunit.h"
#ifdef QDP_IS_QDPJIT
#include "actions/ferm/invert/quda_solvers/qdpjit_memory_wrapper.h"
#endif

//#include <util_quda.h>

namespace Chroma
{

    //! Richardson system solver namespace
    namespace LinOpSysSolverQUDACloverEnv
    {
        //! Register the syssolver
        bool registerAll();
    }



    //! Solve a Clover Fermion System using the QUDA inverter
    /*! \ingroup invert
     *** WARNING THIS SOLVER WORKS FOR Clover FERMIONS ONLY ***
     */

    class LinOpSysSolverQUDAClover : public LinOpSystemSolver<LatticeFermion>
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
            LinOpSysSolverQUDAClover(Handle< LinearOperator<T> > A_,
                    Handle< FermState<T,Q,Q> > state_,
                    const SysSolverQUDACloverParams& invParam_) :
                A(A_), invParam(invParam_), clov(new CloverTermT<T, U>() ), invclov(new CloverTermT<T, U>())
        {
            QDPIO::cout << "LinOpSysSolverQUDAClover:" << std::endl;

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
            QDPIO::cout << "MDAGM Using QDP-JIT gauge order" << std::endl;
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

            // Default for no preconditioner -- may be overwritten based
            // on innerParams
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

            // Default. This may be overrridden later.
            q_gauge_param.reconstruct_precondition = q_gauge_param.reconstruct_sloppy;
            // Gauge fixing:

            // These are the links
            // They may be smeared and the BC's may be applied
            links_single.resize(Nd);

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
                    QDPIO::cerr << "Unknown SOlver type" << std::endl;
                    QDP_abort(1);
                    break;
            }

            // Mass

            // Fiendish idea from Ron. Set the kappa=1/2 and use
            // unmodified clover term, and ask for Kappa normalization
            // This should give us A - (1/2)D as the unpreconditioned operator
            // and probabl 1 - {1/4} A^{-1} D A^{-1} D as the preconditioned
            // op. Apart from the A_oo stuff on the antisymmetric we have
            // nothing to do...
            quda_inv_param.kappa = 0.5;

            // FIXME: If we want QUDA to compute the clover coeff, we need to be able to deal
            // with awfuless of anisotropy
            // The value below is a dummy one.
            quda_inv_param.clover_coeff = 1.0; // Dummy tree level value. Not used
            quda_inv_param.tol = toDouble(invParam.RsdTarget);
            quda_inv_param.maxiter = invParam.MaxIter;
            quda_inv_param.reliable_delta = toDouble(invParam.Delta);
            quda_inv_param.pipeline = invParam.Pipeline;

            // Solution type
            quda_inv_param.solution_type = QUDA_MATPC_SOLUTION;

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
                QDPIO::cout << "Using Asymmetric Linop: A_oo - D A^{-1}_ee D" << std::endl;
                quda_inv_param.matpc_type = QUDA_MATPC_ODD_ODD_ASYMMETRIC;
            }
            else {
                QDPIO::cout << "Using Symmetric Linop: 1 - A^{-1}_oo D A^{-1}_ee D" << std::endl;
                quda_inv_param.matpc_type = QUDA_MATPC_ODD_ODD;
            }

            quda_inv_param.dagger = QUDA_DAG_NO;
            quda_inv_param.mass_normalization = QUDA_KAPPA_NORMALIZATION;

            quda_inv_param.cpu_prec = cpu_prec;
            quda_inv_param.cuda_prec = gpu_prec;
            quda_inv_param.cuda_prec_sloppy = gpu_half_prec;

            // Default. May be overridden by inner params
            quda_inv_param.cuda_prec_precondition = gpu_half_prec;


            quda_inv_param.preserve_source = QUDA_PRESERVE_SOURCE_NO;
            quda_inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;

#ifndef BUILD_QUDA_DEVIFACE_SPINOR
            quda_inv_param.dirac_order = QUDA_DIRAC_ORDER;
#else
            QDPIO::cout << "MDAGM Using QDP-JIT spinor order" << std::endl;
            quda_inv_param.dirac_order    = QUDA_QDPJIT_DIRAC_ORDER;
            quda_inv_param.input_location = QUDA_CUDA_FIELD_LOCATION;
            quda_inv_param.output_location = QUDA_CUDA_FIELD_LOCATION;
#endif


            // Clover precision and order
            quda_inv_param.clover_cpu_prec = cpu_prec;
            quda_inv_param.clover_cuda_prec = gpu_prec;
            quda_inv_param.clover_cuda_prec_sloppy = gpu_half_prec;

            // Default. may be overrridden by inner params
            quda_inv_param.clover_cuda_prec_precondition = gpu_half_prec;

#ifndef BUILD_QUDA_DEVIFACE_CLOVER
#warning "NOT USING QUDA DEVICE IFACE"
            quda_inv_param.clover_order = QUDA_PACKED_CLOVER_ORDER;
#else
#warning "USING QUDA DEVICE IFACE"
            QDPIO::cout << "MDAGM clover CUDA location\n";
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

            if( invParam.innerParamsP ) {
                QDPIO::cout << "Setting inner solver params" << std::endl;
                // Dereference handle
                const GCRInnerSolverParams& ip = *(invParam.innerParams);

                // Set preconditioner precision
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
                q_gauge_param.reconstruct_precondition = QUDA_RECONSTRUCT_NO;
                quda_inv_param.gcrNkrylov = 1;
            }


            if( invParam.verboseP ) {
                quda_inv_param.verbosity = QUDA_VERBOSE;
            }
            else {
                quda_inv_param.verbosity = QUDA_SUMMARIZE;
            }

#ifndef BUILD_QUDA_DEVIFACE_GAUGE
            for(int mu=0; mu < Nd; mu++) {
                gauge[mu] = (void *)&(links_single[mu].elem(all.start()).elem().elem(0,0).real());
            }
#else
            GetMemoryPtrGauge(gauge,links_single);
            // gauge[mu] = GetMemoryPtr( links_single[mu].getId() );
            // QDPIO::cout << "MDAGM CUDA gauge[" << mu << "] in = " << gauge[mu] << "\n";
#endif

            loadGaugeQuda((void *)gauge, &q_gauge_param);

            //      Setup the clover term...
            QDPIO::cout << "Creating CloverTerm" << std::endl;
            clov->create(fstate, invParam_.CloverParams);
            // Don't recompute, just copy
            invclov->create(fstate, invParam_.CloverParams);

            QDPIO::cout << "Inverting CloverTerm" << std::endl;
            invclov->choles(0);
            invclov->choles(1);


#ifndef BUILD_QUDA_DEVIFACE_CLOVER
            packed_clov.resize(all.siteTable().size());
            packed_invclov.resize(all.siteTable().size());
            clov->packForQUDA(packed_clov, 0);
            clov->packForQUDA(packed_clov, 1);

            invclov->packForQUDA(packed_invclov, 0);
            invclov->packForQUDA(packed_invclov, 1);


            loadCloverQuda(&(packed_clov[0]), &(packed_invclov[0]),&quda_inv_param);

#else
            void *clover[2];
            void *cloverInv[2];

            // This is a yucky macro and needs the existence of 'clover' and 'cloverInv' to work
            GetMemoryPtrClover(clov->getOffId(),clov->getDiaId(),invclov->getOffId(),invclov->getDiaId());

            loadCloverQuda( (void*)(clover) , (void*)(cloverInv) ,&quda_inv_param);
#endif



        }


            //! Destructor is automatic
            ~LinOpSysSolverQUDAClover()
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
            SystemSolverResults_t operator() (T& psi, const T& chi) const
            {
                SystemSolverResults_t res;

                START_CODE();
                StopWatch swatch;
                swatch.start();

                //    T MdagChi;

                // This is a CGNE. So create new RHS
                //      (*A)(MdagChi, chi, MINUS);
                // Handle< LinearOperator<T> > MM(new MdagMLinOp<T>(A));
                if ( invParam.axialGaugeP ) {
                    T g_chi,g_psi;

                    // Gauge Fix source and initial guess
                    QDPIO::cout << "Gauge Fixing source and initial guess" << std::endl;
                    g_chi[ rb[1] ]  = GFixMat * chi;
                    g_psi[ rb[1] ]  = GFixMat * psi;
                    QDPIO::cout << "Solving" << std::endl;
                    res = qudaInvert( g_chi, g_psi);
                    QDPIO::cout << "Untransforming solution." << std::endl;
                    psi[ rb[1]]  = adj(GFixMat)*g_psi;

                }
                else {
                    res = qudaInvert( chi, psi);
                }

                swatch.stop();


                {
                    T r;
                    r[A->subset()]=chi;
                    T tmp;
                    (*A)(tmp, psi, PLUS);
                    r[A->subset()] -= tmp;
                    res.resid = sqrt(norm2(r, A->subset()));
                }

                Double rel_resid = res.resid/sqrt(norm2(chi,A->subset()));

                QDPIO::cout << "QUDA_"<< solver_string <<"_CLOVER_SOLVER: " << res.n_count << " iterations. Rsd = " << res.resid << " Relative Rsd = " << rel_resid << std::endl;

                // Convergence Check/Blow Up
                if ( ! invParam.SilentFailP ) {
                    if (  toBool( rel_resid >  invParam.RsdToleranceFactor*invParam.RsdTarget) ) {
                        QDPIO::cerr << "ERROR: QUDA Solver residuum is outside tolerance: QUDA resid="<< rel_resid << " Desired =" << invParam.RsdTarget << " Max Tolerated = " << invParam.RsdToleranceFactor*invParam.RsdTarget << std::endl;
                        QDP_abort(1);
                    }
                }

                END_CODE();
                return res;
            }



            std::vector<SystemSolverResults_t> operator() (const std::vector<std::shared_ptr<T>>& psi, const std::vector<std::shared_ptr<const T>>& chi) const override
            {

                START_CODE();
                QDPIO::cout << "Entering MRHS solution: N_src = " << chi.size() << "\n";

                std::vector<SystemSolverResults_t> res(chi.size());
                if( psi.size() != chi.size() ) {
                    QDPIO::cout << "Number of sources does not match number of solutions\n";
                    QDPIO::cout << "psi.size() = " << psi.size() << " but chi.size() = " << chi.size() << "\n";
                    QDP_abort(1);
                }

                StopWatch swatch;
                swatch.start();

                if ( invParam.axialGaugeP ) {
                    QDPIO::cerr << "Multi RHS solve in axial gauge not yet implemented\n";
                    QDP_abort(1);
                }

                qudaInvertMultiSrc( psi, chi, res);

                swatch.stop();

                // Check solutions
                for(int soln =0; soln < psi.size(); soln++)	{
                    T r;
                    r[A->subset()]=*(chi[ soln ]);
                    T tmp;
                    (*A)(tmp, *(psi[soln]), PLUS);
                    r[A->subset()] -= tmp;
                    res[soln].resid = sqrt(norm2(r, A->subset()));

                    Double rel_resid = res[soln].resid/sqrt(norm2(*(chi[soln]),A->subset()));

                    QDPIO::cout << "QUDA_"<< solver_string <<"_CLOVER_SOLVER: solution " << soln << 
                        " : "  << res[soln].n_count << " iterations. Rsd = " << res[soln].resid 
                        << " Relative Rsd = " << rel_resid << std::endl;

                    // Convergence Check/Blow Up
                    if ( ! invParam.SilentFailP ) {
                        if (  toBool( rel_resid >  invParam.RsdToleranceFactor*invParam.RsdTarget) ) {
                            QDPIO::cerr << "ERROR: QUDA Solver residuum for solution " << soln 
                                << " is outside tolerance: QUDA resid="<< rel_resid << " Desired =" 
                                << invParam.RsdTarget << " Max Tolerated = " 
                                << invParam.RsdToleranceFactor*invParam.RsdTarget << std::endl;
                            QDP_abort(1);
                        }
                    }
                }

                END_CODE();
                return res;
            }

        private:
            // Hide default constructor
            LinOpSysSolverQUDAClover() {}

            Q links_orig;
            Q links_single;
            void *gauge[4];

            U GFixMat;
            QudaPrecision_s cpu_prec;
            QudaPrecision_s gpu_prec;
            QudaPrecision_s gpu_half_prec;

            Handle< LinearOperator<T> > A;
            const SysSolverQUDACloverParams invParam;
            QudaGaugeParam q_gauge_param;
            QudaInvertParam quda_inv_param;

            Handle< CloverTermT<T, U> > clov;
            Handle< CloverTermT<T, U> > invclov;



#ifndef BUILD_QUDA_DEVIFACE_CLOVER
            multi1d<QUDAPackedClovSite<REALT> > packed_clov;
            multi1d<QUDAPackedClovSite<REALT> > packed_invclov;
#endif

            SystemSolverResults_t qudaInvert( const T& chi_s, T& psi_s) const ;

            void qudaInvertMultiSrc( const std::vector<std::shared_ptr<T>>& psi, 
                                     const std::vector<std::shared_ptr<const T>>& chi,
                                     std::vector<SystemSolverResults_t>& res) const;

            std::string solver_string;
    };


} // End namespace

#endif // BUILD_QUDA
#endif 


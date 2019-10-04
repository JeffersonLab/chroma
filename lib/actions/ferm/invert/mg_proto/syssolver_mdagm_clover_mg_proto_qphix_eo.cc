/*
 * solve a MdagM*psi=chi system by mg_proto
 */

#include "chromabase.h"
#include "handle.h"
#include "state.h"
#include "actions/ferm/invert/syssolver_mdagm_factory.h"
#include "actions/ferm/invert/syssolver_mdagm_aggregate.h"
#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/mg_proto/syssolver_mdagm_clover_mg_proto_qphix_eo.h"
#include "lattice/qphix/qphix_blas_wrappers.h"

using namespace QDP;

namespace Chroma
{
    namespace MdagMSysSolverMGProtoQPhiXEOCloverEnv
    {

        //! Anonymous namespace
        namespace
        {
            //! Name to be used
            const std::string name("MG_PROTO_QPHIX_EO_CLOVER_INVERTER");

            //! Local registration flag
            bool registered = false;
        }



        // Double precision
        MdagMSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
                const std::string& path,
                Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state,
                Handle< LinearOperator<LatticeFermion> > A)
        {
            return new MdagMSysSolverMGProtoQPhiXEOClover(A,state,MGProtoSolverParams(xml_in, path));
        }

        //! Register all the factories
        bool registerAll()
        {
            bool success = true;
            if (! registered)
            {

                success &= Chroma::TheMdagMFermSystemSolverFactory::Instance().registerObject(name, createFerm);
                registered = true;
            }
            return success;
        }
    };

    // Save me typing, by exposing this file level from here
    using T = MdagMSysSolverMGProtoQPhiXEOClover::T;
    using Q = MdagMSysSolverMGProtoQPhiXEOClover::Q;


    // Constructor
    MdagMSysSolverMGProtoQPhiXEOClover::MdagMSysSolverMGProtoQPhiXEOClover(Handle< LinearOperator<T> > A_,
            Handle< FermState<T,Q,Q> > state_,
            const MGProtoSolverParams& invParam_) :
        A(A_), state(state_), invParam(invParam_), subspaceId(invParam_.SubspaceId){
            using namespace MGProtoHelpersQPhiX;

            // Create the new Dirac operator with latest gauge field, which is necessary for
            // every MD step in HMC, while keep the same multigrid near null subspace
            // until the solver failed, then recalculate the subspace with current gauge field.
            // This may save time for the multigrid setup stage.
            auto mg_levels = std::make_shared<typename MGPreconditionerEO::LevelT>();

            IndexArray latdims = {{QDP::Layout::subgridLattSize()[0],
                QDP::Layout::subgridLattSize()[1],
                QDP::Layout::subgridLattSize()[2],
                QDP::Layout::subgridLattSize()[3]}};

            auto info = std::make_shared<LatticeInfo>(latdims, 4, 3, NodeInfo());
            M_ptr = createFineLinOpT<typename MGPreconditionerEO::LinOpT>(invParam, state->getLinks(), *info);

            mg_pointer = MGProtoHelpersQPhiX::getMGPreconditionerEO(subspaceId);
            if ( ! mg_pointer ) {
                QDPIO::cout << "EO MG Preconditioner not found in Named Obj. Creating" << std::endl;

                // Check on the links -- they are ferm state and may already have BC's applied? need to figure that out.
                MGProtoHelpersQPhiX::createMGPreconditionerEO(invParam, state->getLinks());

                // Now get the setup
                mg_pointer = MGProtoHelpersQPhiX::getMGPreconditionerEO(subspaceId);
            }

            // Next step is to  create a solver instance:
            fine_solve_params.MaxIter=invParam.OuterSolverMaxIters;
            fine_solve_params.RsdTarget=toDouble(invParam.OuterSolverRsdTarget);
            fine_solve_params.VerboseP =invParam.OuterSolverVerboseP;
            fine_solve_params.NKrylov = invParam.OuterSolverNKrylov;

            // Internal one with EO preconditioning
            using EoFGMRES = const MG::FGMRESSolverQPhiX;

            eo_solver = std::make_shared<const EoFGMRES>(*M_ptr, fine_solve_params, (mg_pointer->v_cycle).get());
        }

    // Destructor
    MdagMSysSolverMGProtoQPhiXEOClover::~MdagMSysSolverMGProtoQPhiXEOClover(){}

    //! Return the subset on which the operator acts
    const Subset&
        MdagMSysSolverMGProtoQPhiXEOClover::subset(void) const
        {
            return A->subset();
        }

    SystemSolverResults_t
        MdagMSysSolverMGProtoQPhiXEOClover::operator()(T& psi, const T& chi) const
        {
            START_CODE();

            QDPIO::cout << "Jolly Greetings from Even-Odd Multigridland" << std::endl;
            const Subset& s = A->subset();
            StopWatch swatch;

            swatch.reset();
            swatch.start();

            QDPIO::cout << "DEBUG: Norm2 Chi Before=" << norm2(chi,s) << std::endl;
            const LatticeInfo& info = M_ptr->GetInfo();
            QPhiXSpinor qphix_in(info);
            QPhiXSpinor qphix_out(info);


            // Solve the MdagM psi = chi by two-step
            // For asymmetric conditioned operator M_a = A_oo - D_oe A^-1_ee D_eo
            // Use M_a^\dagger = \gamma_5 M_a \gamma_5
            //
            // First solve M_a Y = chi', with Y = \gamma_5 M_a psi and chi'=\gamma_5 chi
            // Then Solve M_a psi = Y', with Y' = \gamma_5 Y
            // psi is what we need
            QDPIO::cout<<"MG_PROTO Two Step Solver: "<<std::endl;

            SystemSolverResults_t res;

            QDPIO::cout<<"MG_PROTO_QPHIX_EO_CLOVER_INVERTER: Solve-Y"<<std::endl;
            T chi_prime = Gamma(Nd*Nd - 1) * chi;
            QDPSpinorToQPhiXSpinor(chi_prime, qphix_in);
            ZeroVec(qphix_out,SUBSET_ALL);

            MG::LinearSolverResults res1 = (*eo_solver)(qphix_out, qphix_in, RELATIVE);
            psi = zero;
            QPhiXSpinorToQDPSpinor(qphix_out, psi);

            double qphix_out_norm_cb0 = MG::Norm2Vec(qphix_out, SUBSET_EVEN);
            double qphix_out_norm_cb1 = MG::Norm2Vec(qphix_out, SUBSET_ODD);

            Double psi_norm_cb0 = norm2(psi,rb[0]);
            Double psi_norm_cb1 = norm2(psi,rb[1]);

            Double chi_norm_after = norm2(chi,s);
            QDPIO::cout << "DEBUG: After Solve-Y Norm2 chi = " << chi_norm_after << std::endl;
            QDPIO::cout << "DEBUG: Norm2 qphix_out_cb_0 = " << qphix_out_norm_cb0 << "   Norm psi[0]="<< psi_norm_cb0 << std::endl;
            QDPIO::cout << "DEBUG: Norm2 qphix_out_cb_1 = " << qphix_out_norm_cb1 << "   Norm psi[1]="<< psi_norm_cb1 << std::endl;

            res.n_count = res1.n_count;
            bool solution_good = true;
            {
                // Chroma level check (may be slow)
                T tmp;
                tmp = zero;
                (*A)(tmp, psi, PLUS);

                tmp[s] -= chi_prime;
                Double n2 = norm2(tmp, s);
                Double n2rel = n2 / norm2(chi_prime, s);
                QDPIO::cout << "MG_PROTO_QPHIX_EO_CLOVER_INVERTER: Solve-Y iters = "
                    << res1.n_count << " rel resid = " << sqrt(n2rel) << std::endl;
                if( toBool( sqrt(n2rel) > invParam.OuterSolverRsdTarget * invParam.RsdToleranceFactor ) ) {
                    QDPIO::cout<<"Error in MG_PROTO Solve-Y convergence, retrying..."<<std::endl;
                    solution_good = false;
                }
            }

            // Using subspace from previous MD step will save the setup time,
            // but may increase the total number of iterations, so recreate the subspace
            // if the #iteration exceeds some threshold.
            if(solution_good){
                if(res1.n_count >= invParam.ThresholdCount){
                    QDPIO::cout<<"Solve-Y Iteration Threshold Exceeded! iters = "<<res.n_count<<" Threshold = "<<invParam.ThresholdCount<<std::endl;
                    QDPIO::cout<<"Refreshing Subspace"<<std::endl;

                    StopWatch refresh;
                    refresh.reset();
                    refresh.start();

                    MGProtoHelpersQPhiX::createMGPreconditionerEO(invParam, state->getLinks());
                    mg_pointer = MGProtoHelpersQPhiX::getMGPreconditionerEO(subspaceId);
                    using EoFGMRES = const MG::FGMRESSolverQPhiX;
                    M_ptr = mg_pointer->M;
                    eo_solver = std::make_shared<const EoFGMRES>(*M_ptr, fine_solve_params, (mg_pointer->v_cycle).get());

                    refresh.stop();

                    QDPIO::cout<<"Subspace Refreshing Time = "<<refresh.getTimeInSeconds()<<" secs"<<std::endl;
                }
            } else {
                QDPIO::cout<<"Bazinga! MG_PROTO Solve-Y failed, retry with new multigrid subspace"<<std::endl;

                MGProtoHelpersQPhiX::createMGPreconditionerEO(invParam, state->getLinks());
                mg_pointer = MGProtoHelpersQPhiX::getMGPreconditionerEO(subspaceId);

                using EoFGMRES = const MG::FGMRESSolverQPhiX;
                M_ptr = mg_pointer->M;
                eo_solver = std::make_shared<const EoFGMRES>(*M_ptr, fine_solve_params, (mg_pointer->v_cycle).get());

                ZeroVec(qphix_out,SUBSET_ALL);

                res1 = (*eo_solver)(qphix_out, qphix_in, RELATIVE);
                psi = zero;
                QPhiXSpinorToQDPSpinor(qphix_out, psi);

                double qphix_out_norm_cb0 = MG::Norm2Vec(qphix_out, SUBSET_EVEN);
                double qphix_out_norm_cb1 = MG::Norm2Vec(qphix_out, SUBSET_ODD);

                Double psi_norm_cb0 = norm2(psi,rb[0]);
                Double psi_norm_cb1 = norm2(psi,rb[1]);

                Double chi_norm_after = norm2(chi,s);
                QDPIO::cout << "DEBUG: After Resolve-Y Norm2 chi = " << chi_norm_after << std::endl;
                QDPIO::cout << "DEBUG: Norm2 qphix_out_cb_0 = " << qphix_out_norm_cb0 << "   Norm psi[0]="<< psi_norm_cb0 << std::endl;
                QDPIO::cout << "DEBUG: Norm2 qphix_out_cb_1 = " << qphix_out_norm_cb1 << "   Norm psi[1]="<< psi_norm_cb1 << std::endl;

                {
                    // Chroma level check (may be slow)
                    T tmp;
                    tmp = zero;
                    (*A)(tmp, psi, PLUS);

                    tmp[s] -= chi_prime;
                    Double n2 = norm2(tmp, s);
                    Double n2rel = n2 / norm2(chi_prime, s);
                    QDPIO::cout << "MG_PROTO_QPHIX_EO_CLOVER_INVERTER: Resolve-Y iters = "<< res1.n_count<< " rel resid = " << sqrt(n2rel) << std::endl;
                    if( toBool( sqrt(n2rel) > invParam.OuterSolverRsdTarget * invParam.RsdToleranceFactor ) ) {
                        QDPIO::cout<<"Error in MG_PROTO Resolve-Y convergence, exiting"<<std::endl;
                        MGSolverException convergence_fail(invParam.CloverParams.Mass,
                                subspaceId,
                                res.n_count + res1.n_count,
                                Real(sqrt(n2rel)),
                                invParam.OuterSolverRsdTarget);
                        throw convergence_fail;
                    }
                }
                res.n_count += res1.n_count;
            }

            QDPIO::cout<<"MG_PROTO_QPHIX_EO_CLOVER_INVERTER: Solve-X"<<std::endl;

            T Y_prime = Gamma(Nd*Nd - 1) * psi;
            ZeroVec(qphix_in, SUBSET_ALL);
            QDPSpinorToQPhiXSpinor(Y_prime, qphix_in);

            ZeroVec(qphix_out, SUBSET_ALL);
            MG::LinearSolverResults res2 = (*eo_solver)(qphix_out, qphix_in, RELATIVE);
            psi = zero;
            QPhiXSpinorToQDPSpinor(qphix_out, psi);

            res.n_count += res2.n_count;
            res.resid = res2.resid;
            {
                // Chroma level check (may be slow)
                T tmp, tmp1;
                tmp = zero;
                (*A)(tmp, psi, PLUS);
                tmp1 = zero;
                (*A)(tmp1, tmp, MINUS);

                tmp1[s] -= chi;
                Double n2 = norm2(tmp1, s);
                Double n2rel = n2 / norm2(chi, s);
                QDPIO::cout << "MG_PROTO_QPHIX_EO_CLOVER_INVERTER: Solve-X iters = "
                    << res2.n_count << " rel resid = " << sqrt(n2rel) << std::endl;
                if( toBool( sqrt(n2rel) > invParam.OuterSolverRsdTarget * invParam.RsdToleranceFactor ) ) {
                    QDPIO::cout<<"Error in MG_PROTO Solve-X convergence, retrying..."<<std::endl;
                    solution_good = false;
                }
            }

            // Using subspace from previous MD step will save the setup time,
            // but may increase the total number of iterations,
            // at this step, just delete the existing subspace if #iterations exceeds threshold
            // and next MD step will generate with new gauge field
            if(solution_good){
                if(res2.n_count >= invParam.ThresholdCount){
                    QDPIO::cout<<"Solver-X Iteration Threshold Exceeded! iters = "<<res2.n_count<<" Threshold = "<<invParam.ThresholdCount<<std::endl;
                    QDPIO::cout<<"Deleting Subspace"<<std::endl;

                    StopWatch refresh;
                    refresh.reset();
                    refresh.start();
                    MGProtoHelpersQPhiX::deleteMGPreconditionerEO(subspaceId);
                    refresh.stop();

                    QDPIO::cout<<"Subspace Deleting Time = "<<refresh.getTimeInSeconds()<<" secs"<<std::endl;
                }
            } else {
                QDPIO::cout<<"Bazinga! MG_PROTO Solve-X failed, retry with new multigrid subspace"<<std::endl;

                MGProtoHelpersQPhiX::createMGPreconditionerEO(invParam, state->getLinks());
                mg_pointer = MGProtoHelpersQPhiX::getMGPreconditionerEO(subspaceId);

                using EoFGMRES = const MG::FGMRESSolverQPhiX;
                M_ptr = mg_pointer->M;
                eo_solver = std::make_shared<const EoFGMRES>(*M_ptr, fine_solve_params, (mg_pointer->v_cycle).get());

                ZeroVec(qphix_out, SUBSET_ALL);
                res2 = (*eo_solver)(qphix_out, qphix_in, RELATIVE);
                psi = zero;
                QPhiXSpinorToQDPSpinor(qphix_out, psi);

                double qphix_out_norm_cb0 = MG::Norm2Vec(qphix_out, SUBSET_EVEN);
                double qphix_out_norm_cb1 = MG::Norm2Vec(qphix_out, SUBSET_ODD);

                Double psi_norm_cb0 = norm2(psi,rb[0]);
                Double psi_norm_cb1 = norm2(psi,rb[1]);

                Double chi_norm_after = norm2(chi,s);
                QDPIO::cout << "DEBUG: After Resolve-X Norm2 chi = " << chi_norm_after << std::endl;
                QDPIO::cout << "DEBUG: Norm2 qphix_out_cb_0 = " << qphix_out_norm_cb0 << "   Norm psi[0]="<< psi_norm_cb0 << std::endl;
                QDPIO::cout << "DEBUG: Norm2 qphix_out_cb_1 = " << qphix_out_norm_cb1 << "   Norm psi[1]="<< psi_norm_cb1 << std::endl;

                {
                    // Chroma level check (may be slow)
                    T tmp, tmp1;
                    tmp = zero;
                    (*A)(tmp, psi, PLUS);
                    tmp1 = zero;
                    (*A)(tmp1, tmp, MINUS);

                    tmp1[s] -= chi;
                    Double n2 = norm2(tmp1, s);
                    Double n2rel = n2 / norm2(chi, s);
                    QDPIO::cout << "MG_PROTO_QPHIX_EO_CLOVER_INVERTER: Resolve-X iters = "<<res2.n_count<< " rel resid = " << sqrt(n2rel) << std::endl;
                    if( toBool( sqrt(n2rel) > invParam.OuterSolverRsdTarget * invParam.RsdToleranceFactor ) ) {
                        QDPIO::cout<<"Error in MG_PROTO Resolve-X convergence, exiting"<<std::endl;
                        MGSolverException convergence_fail(invParam.CloverParams.Mass,
                                subspaceId,
                                res2.n_count,
                                Real(sqrt(n2rel)),
                                invParam.OuterSolverRsdTarget);
                        throw convergence_fail;
                    }
                }
                res.n_count += res2.n_count;
                res.resid = res2.resid;
            }

            swatch.stop();
            QDPIO::cout << "MG_PROTO_QPHIX_EO_CLOVER_INVERTER: total_iters=" << res.n_count << " rel resid = " << res.resid << std::endl;
            QDPIO::cout << "MG_PROTO_QPHIX_EO_CLOVER_INVERTER: total_time=" << swatch.getTimeInSeconds() << " sec." << std::endl;

            END_CODE();
            return res;
        }

    SystemSolverResults_t
        MdagMSysSolverMGProtoQPhiXEOClover::operator()(T& psi, const T& chi, AbsChronologicalPredictor4D<T>& predictor) const {
            START_CODE();

            SystemSolverResults_t res;
            res = (*this)(psi, chi);

            END_CODE();
            return res;
        }

};


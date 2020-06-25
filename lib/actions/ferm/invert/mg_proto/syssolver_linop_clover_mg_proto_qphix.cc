/*
 * syssolver_linop_clover_mg_proto.cc
 *
 *  Created on: Mar 23, 2017
 *      Author: bjoo
 */

#include "chromabase.h"
#include "actions/ferm/invert/mg_proto/syssolver_linop_clover_mg_proto_qphix.h"
#include "handle.h"
#include "state.h"
#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/mg_proto/mg_proto_qphix_helpers.h"

#include "lattice/solver.h"
#include "lattice/fgmres_common.h"
#include "lattice/qphix/invfgmres_qphix.h"
#include "lattice/qphix/qphix_qdp_utils.h"
#include "lattice/qphix/qphix_clover_linear_operator.h"
#include "actions/ferm/invert/mg_solver_exception.h"

#include <memory>
using namespace QDP;

namespace Chroma
{
  namespace LinOpSysSolverMGProtoQPhiXCloverEnv
  {

    //! Anonymous namespace
    namespace
    {
      //! Name to be used
      const std::string name("MG_PROTO_QPHIX_CLOVER_INVERTER");

      //! Local registration flag
      bool registered = false;
    }



    // Double precision
    LinOpSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path,
						  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state,

						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new LinOpSysSolverMGProtoQPhiXClover(A,state,MGProtoSolverParams(xml_in, path));
    }

    //! Register all the factories
     bool registerAll()
     {
       bool success = true;
       if (! registered)
       {

    	   success &= Chroma::TheLinOpFermSystemSolverFactory::Instance().registerObject(name, createFerm);
    	   registered = true;
       }
       return success;
     }
  };

  // Save me typing, by exposing this file level from here
  using T = LinOpSysSolverMGProtoQPhiXClover::T;
  using Q = LinOpSysSolverMGProtoQPhiXClover::Q;

  // Constructor
  LinOpSysSolverMGProtoQPhiXClover::LinOpSysSolverMGProtoQPhiXClover(Handle< LinearOperator<T> > A_,
		  Handle< FermState<T,Q,Q> > state_,
		  const MGProtoSolverParams& invParam_) :  A(A_), state(state_), invParam(invParam_) {}

  // Destructor
  LinOpSysSolverMGProtoQPhiXClover::~LinOpSysSolverMGProtoQPhiXClover(){}

  //! Return the subset on which the operator acts
  const Subset&
  LinOpSysSolverMGProtoQPhiXClover::subset(void) const
  {
	  return A->subset();
  }

  SystemSolverResults_t
  LinOpSysSolverMGProtoQPhiXClover::operator()(T& psi, const T& chi) const
  {
	  QDPIO::cout << "Jolly Greetings from Multigridland" << std::endl;
	  StopWatch swatch;
	  StopWatch swatch2;

	  swatch.reset();
	  swatch.start();

	  // Let us see if we have a Multigrid setup lying around.
	  const std::string& subspaceId = invParam.SubspaceId;

	  std::shared_ptr<MGProtoHelpersQPhiX::MGPreconditioner> mg_pointer = MGProtoHelpersQPhiX::getMGPreconditioner(subspaceId);
	  if ( ! mg_pointer ) {
		  QDPIO::cout << "MG Preconditioner not found in Named Obj. Creating" << std::endl;

		  // Check on the links -- they are ferm state and may already have BC's applied? need to figure that out.
		  MGProtoHelpersQPhiX::createMGPreconditioner(invParam, state->getLinks());

		  // Now get the setup
		  mg_pointer = MGProtoHelpersQPhiX::getMGPreconditioner(subspaceId);
	  }

	  // Next step is to  create a solver instance:
	  MG::LinearSolverParamsBase fine_solve_params;
	  fine_solve_params.MaxIter=invParam.OuterSolverMaxIters;
	  fine_solve_params.RsdTarget=toDouble(invParam.OuterSolverRsdTarget);
	  fine_solve_params.VerboseP =invParam.OuterSolverVerboseP;
	  fine_solve_params.NKrylov = invParam.OuterSolverNKrylov;

	  std::shared_ptr<MG::QPhiXWilsonCloverLinearOperator > M_ptr = (mg_pointer->M);

	  const LatticeInfo& info = M_ptr->GetInfo();
	  QPhiXSpinor qphix_in(info);
	  QPhiXSpinor qphix_out(info);

	  MG::FGMRESSolverQPhiX FGMRESOuter(*M_ptr, fine_solve_params, (mg_pointer->v_cycle).get());

	  // Solve the system
	  QDPSpinorToQPhiXSpinor(chi,qphix_in,0);
	  ZeroVec(qphix_out);

	  swatch2.reset();
	  swatch2.start();
	  MG::LinearSolverResults res=FGMRESOuter(qphix_out,qphix_in, RELATIVE)[0];
	  swatch2.stop();

	  QPhiXSpinorToQDPSpinor(qphix_out,0,psi);

	  {
		  T tmp;
		  tmp = zero;
		  (*A)(tmp, psi, PLUS);
		  tmp -= chi;
		  Double n2 = norm2(tmp);
		  Double n2rel = n2 / norm2(chi);
		  QDPIO::cout << "MG_PROTO_CLOVER_INVERTER: iters = "<< res.n_count << " rel resid = " << sqrt(n2rel) << std::endl;
		  if( toBool( sqrt(n2rel) > invParam.OuterSolverRsdTarget ) ) {
		    MGSolverException convergence_fail(invParam.CloverParams.Mass, 
						       subspaceId,
						       res.n_count,
						       Real(sqrt(n2rel)),
						       invParam.OuterSolverRsdTarget);
		    throw convergence_fail;

		  }
	  }
	  swatch.stop();
	  QDPIO::cout << "MG_PROTO_CLOVER_INVERTER_TIME: call_time = "<< swatch2.getTimeInSeconds() << " sec.  total_time=" << swatch.getTimeInSeconds() << " sec." << std::endl;
  }
};


/*
 * syssolver_linop_clover_mg_proto.cc
 *
 *  Created on: Mar 23, 2017
 *      Author: bjoo
 */

#include "chromabase.h"
#include "actions/ferm/invert/mg_proto/syssolver_linop_clover_mg_proto_qphix_eo.h"
#include "handle.h"
#include "state.h"
#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/mg_proto/mg_proto_qphix_helpers.h"

#include "lattice/solver.h"
#include "lattice/fgmres_common.h"
#include "lattice/qphix/invfgmres_qphix.h"
#include "lattice/qphix/qphix_qdp_utils.h"
#include "lattice/qphix/qphix_eo_clover_linear_operator.h"
#include "actions/ferm/invert/mg_solver_exception.h"

#include <memory>
using namespace QDP;

namespace Chroma
{
  namespace LinOpSysSolverMGProtoQPhiXEOCloverEnv
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
    LinOpSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path,
						  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state,

						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new LinOpSysSolverMGProtoQPhiXEOClover(A,state,MGProtoSolverParams(xml_in, path));
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
  using T = LinOpSysSolverMGProtoQPhiXEOClover::T;
  using Q = LinOpSysSolverMGProtoQPhiXEOClover::Q;

  // Constructor
  LinOpSysSolverMGProtoQPhiXEOClover::LinOpSysSolverMGProtoQPhiXEOClover(Handle< LinearOperator<T> > A_,
		  Handle< FermState<T,Q,Q> > state_,
		  const MGProtoSolverParams& invParam_) :  A(A_), state(state_), invParam(invParam_) {}

  // Destructor
  LinOpSysSolverMGProtoQPhiXEOClover::~LinOpSysSolverMGProtoQPhiXEOClover(){}

  //! Return the subset on which the operator acts
  const Subset&
  LinOpSysSolverMGProtoQPhiXEOClover::subset(void) const
  {
	  return A->subset();
  }

  SystemSolverResults_t
  LinOpSysSolverMGProtoQPhiXEOClover::operator()(T& psi, const T& chi) const
  {
	  QDPIO::cout << "Jolly Greetings from Even-Odd Multigridland" << std::endl;
	  StopWatch swatch;
	  StopWatch swatch2;

	  swatch.reset();
	  swatch.start();

	  // Let us see if we have a Multigrid setup lying around.
	  const std::string& subspaceId = invParam.SubspaceId;

	  std::shared_ptr<MGProtoHelpersQPhiX::MGPreconditionerEO> mg_pointer = MGProtoHelpersQPhiX::getMGPreconditionerEO(subspaceId);
	  if ( ! mg_pointer ) {
		  QDPIO::cout << "EO MG Preconditioner not found in Named Obj. Creating" << std::endl;

		  // Check on the links -- they are ferm state and may already have BC's applied? need to figure that out.
		  MGProtoHelpersQPhiX::createMGPreconditionerEO(invParam, state->getLinks());

		  // Now get the setup
		  mg_pointer = MGProtoHelpersQPhiX::getMGPreconditionerEO(subspaceId);
	  }

	  // Next step is to  create a solver instance:
	  MG::FGMRESParams fine_solve_params;
	  fine_solve_params.MaxIter=invParam.OuterSolverMaxIters;
	  fine_solve_params.RsdTarget=toDouble(invParam.OuterSolverRsdTarget);
	  fine_solve_params.VerboseP =invParam.OuterSolverVerboseP;
	  fine_solve_params.NKrylov = invParam.OuterSolverNKrylov;

	  std::shared_ptr<MG::QPhiXWilsonCloverEOLinearOperator > M_ptr = (mg_pointer->M);

	  const LatticeInfo& info = M_ptr->GetInfo();
	  QPhiXSpinor qphix_in(info);
	  QPhiXSpinor qphix_out(info);

	  // Shorthand for the UnprecWrapper
	  using UnprecFGMRES =  MG::UnprecFGMRESSolverQPhiXWrapper;

	  // Internal one with EO preconditioning
	  using EoFGMRES = const MG::FGMRESSolverQPhiX;

	//  std::shared_ptr<const MG::FGMRESSolverQPhiX> fgmres_eo=std::make_shared<const MG::FGMRESSolverQPhiX>(*M_ptr, fine_solve_params, (mg_pointer->v_cycle).get());
	  UnprecFGMRES wrapped(std::make_shared<const EoFGMRES>(*M_ptr, fine_solve_params, (mg_pointer->v_cycle).get()), M_ptr);

	  // Solve the system
	  QDPSpinorToQPhiXSpinor(chi,qphix_in);
	  ZeroVec(qphix_out);

	  swatch2.reset();
	  swatch2.start();
	  MG::LinearSolverResults res=wrapped(qphix_out,qphix_in, RELATIVE);
	  swatch2.stop();

	  QPhiXSpinorToQDPSpinor(qphix_out,psi);

	  {
		  // Chroma level check (may be slow)
		  T tmp;
		  tmp = zero;
		  (*A)(tmp, psi, PLUS);
		  tmp -= chi;
		  Double n2 = norm2(tmp);
		  Double n2rel = n2 / norm2(chi);
		  QDPIO::cout << "MG_PROTO_QPHIX_EO_CLOVER_INVERTER: iters = "<< res.n_count << " rel resid = " << sqrt(n2rel) << std::endl;
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
	  QDPIO::cout << "MG_PROTO_QPHIX_EO_CLOVER_INVERTER_TIME: call_time = "<< swatch2.getTimeInSeconds() << " sec.  total_time=" << swatch.getTimeInSeconds() << " sec." << std::endl;
  }
};


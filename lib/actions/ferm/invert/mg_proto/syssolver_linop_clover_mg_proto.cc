/*
 * syssolver_linop_clover_mg_proto.cc
 *
 *  Created on: Mar 23, 2017
 *      Author: bjoo
 */

#include "chromabase.h"
#include "actions/ferm/invert/mg_proto/syssolver_linop_clover_mg_proto.h"
#include "handle.h"
#include "state.h"
#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/mg_proto/mg_proto_helpers.h"


#include "lattice/fgmres_common.h"
#include "lattice/fine_qdpxx/invfgmres.h"

using namespace QDP;

namespace Chroma
{
  namespace LinOpSysSolverMGProtoCloverEnv
  {

    //! Anonymous namespace
    namespace
    {
      //! Name to be used
      const std::string name("MG_PROTO_CLOVER_INVERTER");

      //! Local registration flag
      bool registered = false;
    }



    // Double precision
    LinOpSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path,
						  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state,

						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new LinOpSysSolverMGProtoClover(A,state,MGProtoSolverParams(xml_in, path));
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
  using T = LinOpSysSolverMGProtoClover::T;
  using Q = LinOpSysSolverMGProtoClover::Q;

  // Constructor
  LinOpSysSolverMGProtoClover::LinOpSysSolverMGProtoClover(Handle< LinearOperator<T> > A_,
		  Handle< FermState<T,Q,Q> > state_,
		  const MGProtoSolverParams& invParam_) :  A(A_), state(state_), invParam(invParam_) {}

  // Destructor
  LinOpSysSolverMGProtoClover::~LinOpSysSolverMGProtoClover(){}

  //! Return the subset on which the operator acts
  const Subset&
  LinOpSysSolverMGProtoClover::subset(void) const
  {
	  return A->subset();
  }

  SystemSolverResults_t
  LinOpSysSolverMGProtoClover::operator()(T& psi, const T& chi) const
  {
	  QDPIO::cout << "Jolly Greetings from Multigridland" << std::endl;

	  // Let us see if we have a Multigrid setup lying around.
	  const std::string& subspaceId = invParam.SubspaceId;

	  std::shared_ptr<MGProtoHelpers::MGPreconditioner> mg_pointer = MGProtoHelpers::getMGPreconditioner(subspaceId);
	  if ( ! mg_pointer ) {
		  QDPIO::cout << "MG Preconditioner not found in Named Obj. Creating" << std::endl;

		  // Check on the links -- they are ferm state and may already have BC's applied? need to figure that out.
		  MGProtoHelpers::createMGPreconditioner(invParam, state->getLinks());

		  // Now get the setup
		  mg_pointer = MGProtoHelpers::getMGPreconditioner(subspaceId);
	  }

	  // Next step is to  create a solver instance:
	  MG::FGMRESParams fine_solve_params;
	  fine_solve_params.MaxIter=invParam.OuterSolverMaxIters;
	  fine_solve_params.RsdTarget=toDouble(invParam.OuterSolverRsdTarget);
	  fine_solve_params.VerboseP =invParam.OuterSolverVerboseP;
	  fine_solve_params.NKrylov = invParam.OuterSolverNKrylov;

	  shared_ptr<const MG::QDPWilsonCloverLinearOperator > M_ptr = (mg_pointer->mg_levels->fine_level).M;

	  MG::FGMRESSolver FGMRESOuter(*M_ptr, fine_solve_params, (mg_pointer->v_cycle).get());

	  // Solve the system
	  FGMRESOuter(psi, chi, RELATIVE);



  }
};


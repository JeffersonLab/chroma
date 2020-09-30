/*
 * syssolver_linop_clover_mg_proto.cc
 *
 *  Created on: Mar 23, 2017
 *      Author: bjoo
 */

#include "chromabase.h"
#include "handle.h"
#include "state.h"
#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/mg_proto/syssolver_linop_ali_mg_proto_qphix_eo.h"

#include "lattice/qphix/qphix_blas_wrappers.h"
#include "MG_config.h"

using namespace QDP;

namespace Chroma
{
  namespace LinOpSysSolverMGProtoQPhiXALIEnv
  {

    //! Anonymous namespace
    namespace
    {
      //! Name to be used
      const std::string name("MG_PROTO_QPHIX_ALI_INVERTER");

      //! Local registration flag
      bool registered = false;
    }



    // Double precision
    LinOpSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path,
						  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state,

						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new LinOpSysSolverMGProtoQPhiXALI(A,state,MGProtoALIPrecParams(xml_in, path));
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
  using T = LinOpSysSolverMGProtoQPhiXALI::T;
  using Q = LinOpSysSolverMGProtoQPhiXALI::Q;


  // Constructor
  LinOpSysSolverMGProtoQPhiXALI::LinOpSysSolverMGProtoQPhiXALI(Handle< LinearOperator<T> > A_,
		  Handle< FermState<T,Q,Q> > state_,
		  const MGProtoALIPrecParams& param_) :
				    A(A_), state(state_), aliprec(MGProtoHelpersQPhiX::createALIPrec(param_, state_->getLinks())), invParam(param_), subspaceId(param_.SubspaceId)
  {
	  // Create double precision operator
	  M = MGProtoHelpersQPhiX::createFineEOLinOp(param_, state_->getLinks(), aliprec->GetInfo());

	  // Next step is to  create a solver instance:
	  fine_solve_params.MaxIter=invParam.OuterSolverMaxIters;
	  fine_solve_params.RsdTarget=toDouble(invParam.OuterSolverRsdTarget);
	  fine_solve_params.VerboseP =invParam.OuterSolverVerboseP;
	  fine_solve_params.NKrylov = invParam.OuterSolverNKrylov;
  }

  // Destructor
  LinOpSysSolverMGProtoQPhiXALI::~LinOpSysSolverMGProtoQPhiXALI(){}

  //! Return the subset on which the operator acts
  const Subset&
  LinOpSysSolverMGProtoQPhiXALI::subset(void) const
  {
	  return A->subset();
  }

  SystemSolverResults_t
  LinOpSysSolverMGProtoQPhiXALI::operator()(T& psi, const T& chi) const
  {
	  return (*this)(std::vector<std::shared_ptr<T>>(1, std::shared_ptr<T>(&psi, [](T* p){})),
	                 std::vector<std::shared_ptr<const T>>(1, std::shared_ptr<const T>(&chi, [](const T* p){})))[0];
  }

  std::vector<SystemSolverResults_t>
  LinOpSysSolverMGProtoQPhiXALI::operator()(const std::vector<std::shared_ptr<T>>& psi, const std::vector<std::shared_ptr<const T>>& chi) const
  {
	  assert(psi.size() == chi.size());
	  int ncols = psi.size();

	  QDPIO::cout << "Jolly Greetings from ALI Multigridland multi-right-hand-side" << std::endl;
	  const Subset& s = A->subset(); 
	  StopWatch swatch;
	  StopWatch swatch2;

	  swatch.reset();
	  swatch.start();


	  //QDPIO::cout << "DEBUG: Norm2 Chi Before=" << norm2(chi,s) << std::endl;
	  const LatticeInfo& info = aliprec->GetInfo();
	  QPhiXSpinor qphix_in(info, ncols);
	  QPhiXSpinor qphix_out(info, ncols);

	  // Solve the system
	  for (int col=0; col<ncols; ++col) QDPSpinorToQPhiXSpinor(*chi[col],qphix_in,col);
	  ZeroVec(qphix_out,SUBSET_ALL);

	  swatch2.reset();
	  swatch2.start();
          MG::Timer::TimerAPI::reset();

	  // Single precision solution
          std::vector<MG::LinearSolverResults> res_f;
	  {
	  QPhiXSpinorF qphix_in_f(info, ncols);
	  QPhiXSpinorF qphix_out_f(info, ncols);
	  ZeroVec(qphix_out_f,SUBSET_ALL);
	  ConvertSpinor(qphix_in,qphix_in_f,SUBSET_ODD);
	  MG::LinearSolverParamsBase solver_params_f(fine_solve_params);
	  solver_params_f.RsdTarget = std::max(3e-6, solver_params_f.RsdTarget);
	  MG::FGMRESSolverQPhiXF solver(*aliprec->GetM(), solver_params_f, aliprec.get());
          res_f = solver(qphix_out_f,qphix_in_f, RELATIVE, MG::InitialGuessNotGiven);
	  assert(res_f.size() == ncols);
	  ConvertSpinor(qphix_out_f,qphix_out,SUBSET_ODD);
	  }

	  // Double precision solution
	  std::vector<MG::LinearSolverResults> res;
	  {
	  MG::FGMRESSolverQPhiX solver(*M, fine_solve_params, aliprec.get());
          res = solver(qphix_out,qphix_in, RELATIVE, MG::InitialGuessGiven);
	  }

	  swatch2.stop();
          MG::Timer::TimerAPI::reportAllTimer();
	
	  for (int col=0; col<ncols; ++col) {
	    *psi[col] = zero;
	    QPhiXSpinorToQDPSpinor(qphix_out,col,*psi[col]);
	  }

	  std::vector<SystemSolverResults_t> r(ncols);
	  for (int col=0; col<ncols; ++col) {
		  // Chroma level check (may be slow)
		  T tmp;
		  tmp = zero;
		  (*A)(tmp, *psi[col], PLUS);
		  tmp[s]  -= *chi[col];
		  Double n2 = norm2(tmp,s);
		  Double n2rel = n2 / norm2(*chi[col],s);
		  r[col].n_count = res_f[col].n_count + res[col].n_count;
		  QDPIO::cout << "MG_PROTO_QPHIX_ALI_INVERTER: iters = "<< r[col].n_count << " rel resid = " << sqrt(n2rel) << std::endl;
		  if( toBool( sqrt(n2rel) > invParam.OuterSolverRsdTarget ) ) {
		    MGSolverException convergence_fail(invParam.CloverParams.Mass, 
		        			       subspaceId,
		        			       res[col].n_count,
		        			       Real(sqrt(n2rel)),
		        			       invParam.OuterSolverRsdTarget);
		    throw convergence_fail;

		  }
	  }
	  swatch.stop();
	  QDPIO::cout << "MG_PROTO_QPHIX_ALI_INVERTER_TIME: call_time = "<< swatch2.getTimeInSeconds() << " sec.  total_time=" << swatch.getTimeInSeconds() << " sec." << std::endl;
	  return r;
  }
};


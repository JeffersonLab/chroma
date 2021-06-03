/*
 * quda_mg_utils.h
 *
 *  Created on: Oct 18, 2016
 *      Author: bjoo
 */

#ifndef LIB_ACTIONS_FERM_INVERT_QUDA_SOLVERS_QUDA_MG_UTILS_H_
#define LIB_ACTIONS_FERM_INVERT_QUDA_SOLVERS_QUDA_MG_UTILS_H_

#include "chromabase.h"

#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/fermstates/simple_fermstate.h"
#include "actions/ferm/invert/quda_solvers/qdpjit_memory_wrapper.h"
#include "actions/ferm/invert/quda_solvers/syssolver_quda_multigrid_clover_params.h"
#include "actions/ferm/linop/clover_term_w.h"
#include "handle.h"
#include "meas/gfix/temporal_gauge.h"
#include "meas/inline/io/named_objmap.h"
#include "state.h"
#include <quda.h>

#include <cuda_runtime_api.h>

namespace Chroma
{

  namespace QUDAMGUtils
  {

    // A Strut for holding pointers essentially
    struct MGSubspacePointers {
      QudaInvertParam mg_inv_param;
      QudaMultigridParam mg_param;
      void* preconditioner;
      MGSubspacePointers()
      {
	mg_inv_param = newQudaInvertParam();
	mg_param = newQudaMultigridParam();
	preconditioner = nullptr;
      }
    };

    template <typename T>
    MGSubspacePointers* create_subspace(const SysSolverQUDAMULTIGRIDCloverParams& invParam)
    {
      MGSubspacePointers* ret_val = new MGSubspacePointers();

      // References so I can keep code below
      QudaInvertParam& mg_inv_param = ret_val->mg_inv_param;
      QudaMultigridParam& mg_param = ret_val->mg_param;

      QudaPrecision_s cpu_prec;
      QudaPrecision_s gpu_prec;
      QudaPrecision_s gpu_half_prec;
      int s = sizeof(typename WordType<T>::Type_t);
      if (s == 4)
      {
	cpu_prec = QUDA_SINGLE_PRECISION;
      }
      else
      {
	cpu_prec = QUDA_DOUBLE_PRECISION;
      }

      // CUDA_PRECISION
      {
	auto found = theChromaToQudaPrecisionTypeMap::Instance().find(invParam.cudaPrecision);
	if (found != theChromaToQudaPrecisionTypeMap::Instance().end())
	{
	  gpu_prec = found->second;
	}
	else
	{
	  // Not found
	  gpu_prec = cpu_prec;
	}
      }

#if 0
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
#endif

      // CUDA_PRECISION
      {
	auto found = theChromaToQudaPrecisionTypeMap::Instance().find(invParam.cudaSloppyPrecision);
	if (found != theChromaToQudaPrecisionTypeMap::Instance().end())
	{
	  gpu_half_prec = found->second;
	}
	else
	{
	  // Not found
	  gpu_half_prec = gpu_prec;
	}
      }

#if 0
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
#endif

      QDPIO::cout << "Creating MG subspace." << std::endl;
      //Taken from various places in the old constructor.

      mg_inv_param.dslash_type = QUDA_CLOVER_WILSON_DSLASH;
      mg_inv_param.inv_type = QUDA_GCR_INVERTER;
      mg_inv_param.tol = 1e-10;
      mg_inv_param.maxiter = 10000;
      mg_inv_param.reliable_delta = 1e-10;
      mg_inv_param.cpu_prec = cpu_prec;
      mg_inv_param.cuda_prec = gpu_prec;
      mg_inv_param.cuda_prec_sloppy = gpu_half_prec;
      mg_inv_param.cuda_prec_precondition = gpu_half_prec;
      //Clover stuff
      mg_inv_param.clover_cpu_prec = cpu_prec;
      mg_inv_param.clover_cuda_prec = gpu_prec;
      mg_inv_param.clover_cuda_prec_sloppy = gpu_half_prec;
      mg_inv_param.clover_cuda_prec_precondition = gpu_half_prec;
      mg_inv_param.clover_order = QUDA_PACKED_CLOVER_ORDER;
      //
      //Done...
      // Autotuning
      if (invParam.tuneDslashP)
      {
	QDPIO::cout << "Enabling MG Dslash Autotuning" << std::endl;
	mg_inv_param.tune = QUDA_TUNE_YES;
      }
      else
      {
	QDPIO::cout << "Disabling MG Dslash Autotuning" << std::endl;
	mg_inv_param.tune = QUDA_TUNE_NO;
      }
      if (invParam.MULTIGRIDParamsP)
      {
	QDPIO::cout << "Setting MULTIGRID solver params" << std::endl;
	// Dereference handle
	const MULTIGRIDSolverParams& ip = *(invParam.MULTIGRIDParams);

	// CUDA_PRECISION
	{
	  auto found = theChromaToQudaPrecisionTypeMap::Instance().find(ip.prec);
	  if (found != theChromaToQudaPrecisionTypeMap::Instance().end())
	  {
	    mg_inv_param.cuda_prec_precondition = found->second;
	    mg_inv_param.clover_cuda_prec_precondition = found->second;
	  }
	  else
	  {
	    // Default:
	    mg_inv_param.cuda_prec_precondition = QUDA_HALF_PRECISION;
	    mg_inv_param.clover_cuda_prec_precondition = QUDA_HALF_PRECISION;
	  }
	}

#if 0
				// Set preconditioner precision
				switch( ip.prec ) {
				case HALF:
					mg_inv_param.cuda_prec_precondition = QUDA_HALF_PRECISION;
					mg_inv_param.clover_cuda_prec_precondition = QUDA_HALF_PRECISION;
					break;
				case SINGLE:
					mg_inv_param.cuda_prec_precondition = QUDA_SINGLE_PRECISION;
					mg_inv_param.clover_cuda_prec_precondition = QUDA_SINGLE_PRECISION;
					break;
				case DOUBLE:
					mg_inv_param.cuda_prec_precondition = QUDA_DOUBLE_PRECISION;
					mg_inv_param.clover_cuda_prec_precondition = QUDA_DOUBLE_PRECISION;
					break;
				default:
					mg_inv_param.cuda_prec_precondition = QUDA_HALF_PRECISION;
					mg_inv_param.clover_cuda_prec_precondition = QUDA_HALF_PRECISION;
					break;
				}
#endif
	mg_inv_param.gcrNkrylov = ip.precond_gcr_nkrylov;
	if (ip.verbosity == true)
	{
	  mg_inv_param.verbosity = QUDA_VERBOSE;
	}
	else
	{
	  mg_inv_param.verbosity = QUDA_SUMMARIZE;
	}

	mg_inv_param.verbosity_precondition = QUDA_SILENT;

	mg_inv_param.sp_pad = 0;
	mg_inv_param.cl_pad = 0;
	mg_inv_param.preserve_source = QUDA_PRESERVE_SOURCE_NO;
	mg_inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
	mg_inv_param.dirac_order = QUDA_DIRAC_ORDER;
	mg_inv_param.input_location = QUDA_CPU_FIELD_LOCATION;
	mg_inv_param.output_location = QUDA_CPU_FIELD_LOCATION;
	//		mg_inv_param.use_init_guess = QUDA_USE_INIT_GUESS_NO;
	mg_inv_param.dagger = QUDA_DAG_NO;
	mg_inv_param.kappa = 0.5;
	mg_inv_param.mass_normalization = QUDA_KAPPA_NORMALIZATION;
	mg_inv_param.clover_coeff = 1.0; // Dummy not used
	mg_inv_param.matpc_type = QUDA_MATPC_ODD_ODD;
	mg_inv_param.solution_type = QUDA_MAT_SOLUTION;
	mg_inv_param.solve_type = QUDA_DIRECT_SOLVE;
	mg_param.invert_param = &mg_inv_param;
	mg_inv_param.Ls = 1;
	mg_param.n_level = ip.mg_levels;
	mg_param.run_low_mode_check = QUDA_BOOLEAN_NO;
	mg_param.run_oblique_proj_check = QUDA_BOOLEAN_NO;

	if (ip.check_multigrid_setup == true)
	{
	  mg_param.run_verify = QUDA_BOOLEAN_YES;
	}
	else
	{
	  mg_param.run_verify = QUDA_BOOLEAN_NO;
	}

	for (int i = 0; i < mg_param.n_level - 1; ++i)
	{
	  if (ip.setup_on_gpu[i])
	  {
	    mg_param.setup_location[i] = QUDA_CUDA_FIELD_LOCATION;
	  }
	  else
	  {
	    mg_param.setup_location[i] = QUDA_CPU_FIELD_LOCATION;
	  }
	}

	for (int i = 0; i < mg_param.n_level; i++)
	{
	  for (int j = 0; j < Nd; j++)
	  {
	    if (i < mg_param.n_level - 1)
	    {
	      mg_param.geo_block_size[i][j] = ip.blocking[i][j];
	    }
	    else
	    {
	      mg_param.geo_block_size[i][j] = 4;
	    }
	  }

	  QDPIO::cout << "Set Level " << i << " blocking as: ";
	  for (int j = 0; j < 4; ++j)
	  {
	    QDPIO::cout << " " << mg_param.geo_block_size[i][j];
	  }
	  QDPIO::cout << std::endl;

	  mg_param.spin_block_size[i] = 1;
	  // FIXME: Elevate ip.nvec, ip.nu_pre, ip.nu_post, ip.tol to arrays in the XML
	  if (i < mg_param.n_level - 1)
	  {
	    mg_param.n_vec[i] = ip.nvec[i];
	    mg_param.nu_pre[i] = ip.nu_pre[i];
	    mg_param.nu_post[i] = ip.nu_post[i];
	  }

	  mg_param.mu_factor[i] = 1.0; // default is one in QUDA test program

	  // Hardwire setup solver now
	  if (i < mg_param.n_level - 1)
	  {
	    mg_param.setup_inv_type[i] =
	      theChromaToQudaSolverTypeMap::Instance()[ip.subspaceSolver[i]];
	    mg_param.setup_tol[i] = toDouble(ip.rsdTargetSubspaceCreate[i]);
	    mg_param.setup_maxiter[i] = ip.maxIterSubspaceCreate[i];
	    mg_param.setup_maxiter_refresh[i] =
	      ip.maxIterSubspaceRefresh[i]; // Will set this from outside...
	    mg_param.num_setup_iter[i] = 1; // 1 refine for now
	    mg_param.precision_null[i] = mg_inv_param.cuda_prec_precondition;
	  }

	  // FIXME: Allow coarse_solver as an XML Parameter
	  // FIXME: What about bottom solver. Always GCR or should I make it BiCGStab?

	  if (i > 0)
	  {

	    switch (ip.coarseSolverType[i - 1])
	    {
	    case GCR: mg_param.coarse_solver[i] = QUDA_GCR_INVERTER; break;
	    case CA_GCR:
	      mg_param.coarse_solver[i] = QUDA_CA_GCR_INVERTER;
	      if (i != mg_param.n_level - 1)
	      {
		QDPIO::cout << "ERROR: Right now CA_GCR is only allowed on the bottom level as a "
			       "coarse solver"
			    << std::endl;
		QDP_abort(1);
	      }
	      break;
	    default:
	      QDPIO::cout
		<< "Invalid coarse solver. Right now only GCR and CA_GCR solvers are allowed"
		<< std::endl;
	      break;
	    };
	  }
	  else
	  {
	    // Level  0 isolver is the outer solver. So this parameter is ignored. Set it to something sensible.
	    mg_param.coarse_solver[0] = QUDA_GCR_INVERTER;
	  }

	  // However maxiter and tol are used on level 0 to set precond_tol and precond_maxiter
	  mg_param.coarse_solver_tol[i] = toDouble(ip.tol[i]);
	  mg_param.coarse_solver_maxiter[i] = ip.maxIterations[i];

	  // Smoother Type is needed for every level as we may want to use a smoother
	  // as a preconditioner on the coarsest level too (even tho we don't recurse to yet another leve0
	  switch (ip.smootherType[i])
	  {
	  case MR:
	    mg_param.smoother[i] = QUDA_MR_INVERTER;
	    mg_param.smoother_tol[i] = toDouble(ip.smootherTol[i]);
	    mg_param.smoother_solve_type[i] = QUDA_DIRECT_PC_SOLVE;
	    mg_param.omega[i] = toDouble(ip.relaxationOmegaMG[i]);
	    mg_param.smoother_schwarz_type[i] =
	      theChromaToQudaSchwarzTypeMap::Instance()[ip.smootherSchwarzType[i]];
	    mg_param.smoother_schwarz_cycle[i] = ip.smootherSchwarzCycle[i];
	    break;
	  case CA_GCR:
	    mg_param.smoother[i] = QUDA_CA_GCR_INVERTER;
	    mg_param.smoother_tol[i] = toDouble(ip.smootherTol[i]);
	    mg_param.smoother_solve_type[i] = QUDA_DIRECT_PC_SOLVE;
	    mg_param.omega[i] = toDouble(ip.relaxationOmegaMG[i]);
	    mg_param.smoother_schwarz_type[i] =
	      theChromaToQudaSchwarzTypeMap::Instance()[ip.smootherSchwarzType[i]];
	    mg_param.smoother_schwarz_cycle[i] = ip.smootherSchwarzCycle[i];
	    break;
	  default:
	    QDPIO::cout
	      << "Unknown or no smother type specified, no smoothing inverter will be used."
	      << std::endl;
	    mg_param.smoother[i] = QUDA_INVALID_INVERTER;
	    QDP_abort(1);
	    break;
	  }

	  // if the value is DEFAULT -- leave the smoother halo precision unset.
	  if (ip.smootherHaloPrecision[i] != DEFAULT)
	  {
	    mg_param.smoother_halo_precision[i] =
	      (theChromaToQudaPrecisionTypeMap::Instance())[ip.smootherHaloPrecision[i]];
	  }
#if 0
					if( ip.smootherHaloPrecision[i] != DEFAULT ) {

						switch( ip.smootherHaloPrecision[i] ) {
						case QUARTER :
							mg_param.smoother_halo_precision[i] = QUDA_QUARTER_PRECISION;
							break;
						case HALF :
							mg_param.smoother_halo_precision[i] = QUDA_HALF_PRECISION;
							break;
						case SINGLE :
							mg_param.smoother_halo_precision[i] = QUDA_SINGLE_PRECISION;
							break;
						case DOUBLE :
							mg_param.smoother_halo_precision[i] = QUDA_DOUBLE_PRECISION;
							break;
						default:
							// Leave unset -- default behaviour -- should never be reached
							break;
						}
					}
					QDPIO::cout << "CKheckopoint 2" << std::endl << std::flush;
#endif
	  mg_param.global_reduction[i] = (mg_param.smoother_schwarz_type[i] == QUDA_INVALID_SCHWARZ)
					   ? QUDA_BOOLEAN_YES
					   : QUDA_BOOLEAN_NO;
	  mg_param.location[i] = QUDA_CUDA_FIELD_LOCATION;

	  if (ip.cycle_type == "MG_VCYCLE")
	  {
	    mg_param.cycle_type[i] = QUDA_MG_CYCLE_VCYCLE;
	  }
	  else if (ip.cycle_type == "MG_RECURSIVE")
	  {
	    mg_param.cycle_type[i] = QUDA_MG_CYCLE_RECURSIVE;
	  }
	  else
	  {
	    QDPIO::cout << "Unknown Cycle Type" << ip.cycle_type << std::endl;
	    QDP_abort(1);
	  }

	  switch (mg_param.cycle_type[i])
	  {
	  case QUDA_MG_CYCLE_RECURSIVE:
	    mg_param.coarse_grid_solution_type[i] = QUDA_MATPC_SOLUTION;
	    break;
	  case QUDA_MG_CYCLE_VCYCLE:
	    mg_param.coarse_grid_solution_type[i] = QUDA_MAT_SOLUTION;
	    break;
	  default:
	    QDPIO::cerr << "Should never get here" << std::endl;
	    QDP_abort(1);
	    break;
	  }
	}
	mg_param.setup_type = QUDA_NULL_VECTOR_SETUP;
	mg_param.pre_orthonormalize = QUDA_BOOLEAN_NO;
	mg_param.post_orthonormalize = QUDA_BOOLEAN_YES;

	// LEvel 0 must always be matpc
	mg_param.coarse_grid_solution_type[0] = QUDA_MATPC_SOLUTION;
	// only coarsen the spin on the first restriction
	mg_param.spin_block_size[0] = 2;
	// coarse grid solver is GCR
	mg_param.smoother[ip.mg_levels - 1] = QUDA_GCR_INVERTER;
	mg_param.compute_null_vector =
	  ip.generate_nullspace ? QUDA_COMPUTE_NULL_VECTOR_YES : QUDA_COMPUTE_NULL_VECTOR_NO;
	mg_param.generate_all_levels = ip.generate_all_levels ? QUDA_BOOLEAN_YES : QUDA_BOOLEAN_NO;
	for (int l = 0; l < ip.mg_levels; ++l)
	{
	  mg_param.vec_infile[l][0] = '\0';
	  mg_param.vec_outfile[l][0] = '\0';
	}
	QDPIO::cout << "Basic MULTIGRID params copied." << std::endl;
      }
      // setup the multigrid solver
      // this allocates memory
      QDPIO::cout << "About to Call newMultigridQuda" << std::endl;

      ret_val->preconditioner = newMultigridQuda(&mg_param);
      QDPIO::cout << "NewMultigridQuda state initialized." << std::endl;
      QDPIO::cout << "MULTIGRID preconditioner set." << std::endl;
      QDPIO::cout << "After call to newMultigridQuda" << std::endl;

      for (int i = 0; i < mg_param.n_level; ++i)
      {
	QDPIO::cout << "Set Level " << i << " blocking as: ";

	for (int j = 0; j < 4; ++j)
	{
	  QDPIO::cout << " " << mg_param.geo_block_size[i][j];
	}
	QDPIO::cout << std::endl;
      }

#if 1
      QDPIO::cout << "MULTIGrid Param Dump" << std::endl;
      printQudaMultigridParam(&mg_param);
#endif
      // We have just refreshed so not due for refresh.
      return ret_val;
    }

    inline void delete_subspace(const std::string SubspaceID)
    {
      // REcover pointer to subspace pointers array
      MGSubspacePointers* subspace_pointers =
	TheNamedObjMap::Instance().getData<MGSubspacePointers*>(SubspaceID);

      // Nuke the preconditioner
      destroyMultigridQuda(subspace_pointers->preconditioner);

      // Delete the storage (this automatically takes care the mg_params and mg_inv_param
      // Which are just member structs.
      delete subspace_pointers;

      // Erase the entry
      TheNamedObjMap::Instance().erase(SubspaceID);
    }

    inline size_t getCUDAFreeMem(void)
    {
      cudaError_t ret;
      size_t free, total;
      ret = cudaMemGetInfo(&free, &total);
      if (ret != cudaSuccess)
      {
	QDPIO::cout << "cudaMemGetInfo() returned unsuccesful status: " << ret << std::endl;
	QDP_abort(1);
      }
      return free;
    }

    struct AbstractQUDAClover
    {
      using T = LatticeFermion;
      using U = LatticeColorMatrix;
      using Q = multi1d<LatticeColorMatrix>;

      SysSolverQUDAMULTIGRIDCloverParams invParam;
      Handle<FermState<T, Q, Q>> state, fstate;

      U GFixMat;

      QudaGaugeParam& gauge_param()
      {
	static QudaGaugeParam params;
	return params;
      }

      QudaInvertParam quda_inv_param;

#ifndef BUILD_QUDA_DEVIFACE_CLOVER
      multi1d<QUDAPackedClovSite<REALT>> packed_clov;
      multi1d<QUDAPackedClovSite<REALT>> packed_invclov(all.siteTable().size());
#else
      CloverTermT<T, U> clov;
      CloverTermT<T, U> invclov;
#endif

      MGSubspacePointers* subspace_pointers;

      AbstractQUDAClover(Handle<FermState<T, Q, Q>> state,
			 const SysSolverQUDAMULTIGRIDCloverParams& invParam)
	: state(state), invParam(invParam), subspace_pointers(nullptr)
      {
	load();
      }

      void load()
      {
	// If this state was the last loaded don't load it again
	static Handle<FermState<T, Q, Q>> lastLoadedState;
	static bool lastLoadedAxialGaugeP;
	static QudaPrecisionType lastLoadedPrec;
	static QudaReconsType lastLoadedReconstruct;

	StopWatch init_swatch;
	init_swatch.reset();
	init_swatch.start();

	std::string solver_string = "AbstractQUDAClover";

	// FOLLOWING INITIALIZATION in test QUDA program

	// 1) work out cpu_prec, cuda_prec, cuda_prec_sloppy
	QudaPrecision_s cpu_prec;
	QudaPrecision_s gpu_prec;
	QudaPrecision_s gpu_half_prec;

	int s = sizeof(WordType<T>::Type_t);
	if (s == 4)
	{
	  cpu_prec = QUDA_SINGLE_PRECISION;
	}
	else
	{
	  cpu_prec = QUDA_DOUBLE_PRECISION;
	}

	// Work out GPU precision
	switch (invParam.cudaPrecision)
	{
	case HALF: gpu_prec = QUDA_HALF_PRECISION; break;
	case SINGLE: gpu_prec = QUDA_SINGLE_PRECISION; break;
	case DOUBLE: gpu_prec = QUDA_DOUBLE_PRECISION; break;
	default: gpu_prec = cpu_prec; break;
	}

	// Work out GPU Sloppy precision
	// Default: No Sloppy
	switch (invParam.cudaSloppyPrecision)
	{
	case HALF: gpu_half_prec = QUDA_HALF_PRECISION; break;
	case SINGLE: gpu_half_prec = QUDA_SINGLE_PRECISION; break;
	case DOUBLE: gpu_half_prec = QUDA_DOUBLE_PRECISION; break;
	default: gpu_half_prec = gpu_prec; break;
	}

	// 2) pull 'new; GAUGE and Invert params
	quda_inv_param = newQudaInvertParam();

	// 3) set lattice size
	const multi1d<int>& latdims = Layout::subgridLattSize();

	//
	// Load gauge
	//

	if (lastLoadedState.operator->() != state.operator->() ||
	    lastLoadedAxialGaugeP != invParam.axialGaugeP ||
	    (invParam.MULTIGRIDParams.operator->() &&
	     (lastLoadedPrec != invParam.MULTIGRIDParams->prec ||
	      lastLoadedReconstruct != invParam.MULTIGRIDParams->reconstruct)))
	{
	  QudaGaugeParam& q_gauge_param = gauge_param();
	  q_gauge_param = newQudaGaugeParam();
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
	  if (invParam.AntiPeriodicT)
	  {
	    q_gauge_param.t_boundary = QUDA_ANTI_PERIODIC_T;
	  }
	  else
	  {
	    q_gauge_param.t_boundary = QUDA_PERIODIC_T;
	  }

	  // Set cpu_prec, cuda_prec, reconstruct and sloppy versions
	  q_gauge_param.cpu_prec = cpu_prec;
	  q_gauge_param.cuda_prec = gpu_prec;

	  switch (invParam.cudaReconstruct)
	  {
	  case RECONS_NONE: q_gauge_param.reconstruct = QUDA_RECONSTRUCT_NO; break;
	  case RECONS_8: q_gauge_param.reconstruct = QUDA_RECONSTRUCT_8; break;
	  case RECONS_12: q_gauge_param.reconstruct = QUDA_RECONSTRUCT_12; break;
	  default: q_gauge_param.reconstruct = QUDA_RECONSTRUCT_12; break;
	  };

	  q_gauge_param.cuda_prec_sloppy = gpu_half_prec;

	  // Default. This may be overwritten by inner params
	  q_gauge_param.cuda_prec_precondition = gpu_half_prec;

	  switch (invParam.cudaSloppyReconstruct)
	  {
	  case RECONS_NONE: q_gauge_param.reconstruct_sloppy = QUDA_RECONSTRUCT_NO; break;
	  case RECONS_8: q_gauge_param.reconstruct_sloppy = QUDA_RECONSTRUCT_8; break;
	  case RECONS_12: q_gauge_param.reconstruct_sloppy = QUDA_RECONSTRUCT_12; break;
	  default: q_gauge_param.reconstruct_sloppy = QUDA_RECONSTRUCT_12; break;
	  };
	  q_gauge_param.reconstruct_precondition = q_gauge_param.reconstruct_sloppy;
	  // Gauge fixing:

	  // These are the links
	  // They may be smeared and the BC's may be applied
	  Q links_single(Nd);

	  // Now downcast to single prec fields.
	  for (int mu = 0; mu < Nd; mu++)
	  {
	    links_single[mu] = state->getLinks()[mu];
	  }

	  // GaugeFix
	  if (invParam.axialGaugeP)
	  {
	    temporalGauge(links_single, GFixMat, Nd - 1);
	    for (int mu = 0; mu < Nd; mu++)
	    {
	      links_single[mu] =
		GFixMat * (state->getLinks())[mu] * adj(shift(GFixMat, FORWARD, mu));
	    }
	    q_gauge_param.gauge_fix = QUDA_GAUGE_FIXED_YES;
	  }
	  else
	  {
	    // No GaugeFix
	    q_gauge_param.gauge_fix = QUDA_GAUGE_FIXED_NO; // No Gfix yet
	  }

	  // deferred 4) Gauge Anisotropy
	  const AnisoParam_t& aniso = invParam.CloverParams.anisoParam;
	  if (aniso.anisoP)
	  { // Anisotropic case
	    Real gamma_f = aniso.xi_0 / aniso.nu;
	    q_gauge_param.anisotropy = toDouble(gamma_f);
	  }
	  else
	  {
	    q_gauge_param.anisotropy = 1.0;
	  }

	  // MAKE FSTATE BEFORE RESCALING links_single
	  // Because the clover term expects the unrescaled links...
	  fstate = Handle<FermState<T, Q, Q>>(new PeriodicFermState<T, Q, Q>(links_single));

	  if (aniso.anisoP)
	  { // Anisotropic case
	    multi1d<Real> cf = makeFermCoeffs(aniso);
	    for (int mu = 0; mu < Nd; mu++)
	    {
	      links_single[mu] *= cf[mu];
	    }
	  }
	  // Setup padding

	  multi1d<int> face_size(4);
	  face_size[0] = latdims[1] * latdims[2] * latdims[3] / 2;
	  face_size[1] = latdims[0] * latdims[2] * latdims[3] / 2;
	  face_size[2] = latdims[0] * latdims[1] * latdims[3] / 2;
	  face_size[3] = latdims[0] * latdims[1] * latdims[2] / 2;

	  int max_face = face_size[0];
	  for (int i = 1; i <= 3; i++)
	  {
	    if (face_size[i] > max_face)
	    {
	      max_face = face_size[i];
	    }
	  }

	  q_gauge_param.ga_pad = max_face;

	  // Remember the last gauge loaded
	  lastLoadedState = state;
	  lastLoadedAxialGaugeP = invParam.axialGaugeP;

	  if (invParam.MULTIGRIDParams.operator->())
	  {
	    // Set preconditioner precision
	    switch (invParam.MULTIGRIDParams->prec)
	    {
	    case HALF: q_gauge_param.cuda_prec_precondition = QUDA_HALF_PRECISION; break;
	    case SINGLE: q_gauge_param.cuda_prec_precondition = QUDA_SINGLE_PRECISION; break;
	    case DOUBLE: q_gauge_param.cuda_prec_precondition = QUDA_DOUBLE_PRECISION; break;
	    default: q_gauge_param.cuda_prec_precondition = QUDA_HALF_PRECISION; break;
	    }

	    switch (invParam.MULTIGRIDParams->reconstruct)
	    {
	    case RECONS_NONE: q_gauge_param.reconstruct_precondition = QUDA_RECONSTRUCT_NO; break;
	    case RECONS_8: q_gauge_param.reconstruct_precondition = QUDA_RECONSTRUCT_8; break;
	    case RECONS_12: q_gauge_param.reconstruct_precondition = QUDA_RECONSTRUCT_12; break;
	    default: q_gauge_param.reconstruct_precondition = QUDA_RECONSTRUCT_12; break;
	    }

	    lastLoadedPrec = invParam.MULTIGRIDParams->prec;
	    lastLoadedReconstruct = invParam.MULTIGRIDParams->reconstruct;
	  }

	  // Set up the links
	  void* gauge[4];

#ifndef BUILD_QUDA_DEVIFACE_GAUGE
	  for (int mu = 0; mu < Nd; mu++)
	  {
	    gauge[mu] = (void*)&(links_single[mu].elem(all.start()).elem().elem(0, 0).real());
	  }
#else
	  //gauge[mu] = GetMemoryPtr( links_single[mu].getId() );
	  GetMemoryPtrGauge(gauge, links_single);
#endif

	  loadGaugeQuda((void*)gauge, &q_gauge_param);
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
	if (invParam.tuneDslashP)
	{
	  quda_inv_param.tune = QUDA_TUNE_YES;
	}
	else
	{
	  quda_inv_param.tune = QUDA_TUNE_NO;
	}

	// PADDING
	quda_inv_param.sp_pad = 0;
	quda_inv_param.cl_pad = 0;

	// Clover precision and order
	quda_inv_param.clover_cpu_prec = cpu_prec;
	quda_inv_param.clover_cuda_prec = gpu_prec;
	quda_inv_param.clover_cuda_prec_sloppy = gpu_half_prec;
	quda_inv_param.clover_cuda_prec_precondition = gpu_half_prec;
	if (invParam.MULTIGRIDParamsP)
	{
	  const MULTIGRIDSolverParams& ip = *(invParam.MULTIGRIDParams);

	  // Set preconditioner precision
	  switch (ip.prec)
	  {
	  case HALF:
	    quda_inv_param.cuda_prec_precondition = QUDA_HALF_PRECISION;
	    quda_inv_param.clover_cuda_prec_precondition = QUDA_HALF_PRECISION;
	    break;

	  case SINGLE:
	    quda_inv_param.cuda_prec_precondition = QUDA_SINGLE_PRECISION;
	    quda_inv_param.clover_cuda_prec_precondition = QUDA_SINGLE_PRECISION;
	    break;

	  case DOUBLE:
	    quda_inv_param.cuda_prec_precondition = QUDA_DOUBLE_PRECISION;
	    quda_inv_param.clover_cuda_prec_precondition = QUDA_DOUBLE_PRECISION;
	    break;
	  default:
	    quda_inv_param.cuda_prec_precondition = QUDA_HALF_PRECISION;
	    quda_inv_param.clover_cuda_prec_precondition = QUDA_HALF_PRECISION;
	    break;
	  }

	}

	MULTIGRIDSolverParams ip = *(invParam.MULTIGRIDParams);
	//
	quda_inv_param.tol_precondition = toDouble(ip.tol[0]);
	quda_inv_param.maxiter_precondition = ip.maxIterations[0];
	quda_inv_param.gcrNkrylov = ip.outer_gcr_nkrylov;
	quda_inv_param.residual_type = static_cast<QudaResidualType>(QUDA_L2_RELATIVE_RESIDUAL);

	//Replacing above with what's in the invert test.
	switch (ip.schwarzType)
	{
	case ADDITIVE_SCHWARZ: quda_inv_param.schwarz_type = QUDA_ADDITIVE_SCHWARZ; break;
	case MULTIPLICATIVE_SCHWARZ:
	  quda_inv_param.schwarz_type = QUDA_MULTIPLICATIVE_SCHWARZ;
	  break;
	default: quda_inv_param.schwarz_type = QUDA_ADDITIVE_SCHWARZ; break;
	}
	quda_inv_param.precondition_cycle = 1;
	//Invert test always sets this to 1.

	if (invParam.verboseP)
	{
	  quda_inv_param.verbosity = QUDA_VERBOSE;
	}
	else
	{
	  quda_inv_param.verbosity = QUDA_SUMMARIZE;
	}

	quda_inv_param.verbosity_precondition = QUDA_SILENT;

	quda_inv_param.inv_type_precondition = QUDA_MG_INVERTER;

	QDPIO::cout << solver_string << "Basic MULTIGRID params copied." << std::endl;

	//      Setup the clover term...
#ifndef BUILD_QUDA_DEVIFACE_CLOVER
	CloverTermT<T, U> clov;
	CloverTermT<T, U> invclov;
#endif
	QDPIO::cout << solver_string << "Creating CloverTerm" << std::endl;

	clov.create(fstate, invParam.CloverParams);

	// Don't recompute, just copy
	invclov.create(fstate, invParam.CloverParams);

	QDPIO::cout << solver_string << "Inverting CloverTerm" << std::endl;
	invclov.choles(0);
	invclov.choles(1);

#ifndef BUILD_QUDA_DEVIFACE_CLOVER
#  warning "NOT USING QUDA DEVICE IFACE"
	quda_inv_param.clover_order = QUDA_PACKED_CLOVER_ORDER;

	packed_clov.resize(all.siteTable().size());

	clov.packForQUDA(packed_clov, 0);
	clov.packForQUDA(packed_clov, 1);

	// Always need inverse
	packed_invclov.resize(all.siteTable().size());
	invclov.packForQUDA(packed_invclov, 0);
	invclov.packForQUDA(packed_invclov, 1);

	loadCloverQuda(&(packed_clov[0]), &(packed_invclov[0]), &quda_inv_param);

#else

#  warning "USING QUDA DEVICE IFACE"
	quda_inv_param.clover_location = QUDA_CUDA_FIELD_LOCATION;
	quda_inv_param.clover_order = QUDA_QDPJIT_CLOVER_ORDER;

	void* clover[2];
	void* cloverInv[2];

	GetMemoryPtrClover(clov.getOffId(), clov.getDiaId(), invclov.getOffId(),
			   invclov.getDiaId());

	loadCloverQuda((void*)(clover), (void*)(cloverInv), &quda_inv_param);

#endif

	quda_inv_param.omega = toDouble(ip.relaxationOmegaOuter);

	subspace_pointers = create_subspace<T>(invParam);
	quda_inv_param.preconditioner = subspace_pointers->preconditioner;
      }
    };
  } // MG Utils
} // Chroma

#endif /* LIB_ACTIONS_FERM_INVERT_QUDA_SOLVERS_QUDA_MG_UTILS_H_ */

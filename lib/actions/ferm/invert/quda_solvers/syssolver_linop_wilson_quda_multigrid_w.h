
// -*- C++ -*-
/*! \file
 *  \QUDA MULTIGRID Wilson solver.
 */

#ifndef __syssolver_linop_quda_multigrid_wilson_h__
#define __syssolver_linop_quda_multigrid_wilson_h__

#include "chroma_config.h"

#ifdef BUILD_QUDA
#include <quda.h>

#include "handle.h"
#include "state.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/fermbcs/simple_fermbc.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/invert/quda_solvers/syssolver_quda_multigrid_wilson_params.h"
#include "meas/gfix/temporal_gauge.h"
#include "io/aniso_io.h"
#include <string>

#include "util/gauge/reunit.h"

//#include <util_quda.h>

namespace Chroma
{

  //! Richardson system solver namespace
  namespace LinOpSysSolverQUDAMULTIGRIDWilsonEnv
  {
    //! Register the syssolver
    bool registerAll();
  }



  //! Solve a Wilson Fermion System using the QUDA inverter
  /*! \ingroup invert
 *** WARNING THIS SOLVER WORKS FOR Wilson FERMIONS ONLY ***
   */
 
  class LinOpSysSolverQUDAMULTIGRIDWilson : public LinOpSystemSolver<LatticeFermion>
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
    LinOpSysSolverQUDAMULTIGRIDWilson(Handle< LinearOperator<T> > A_,
					 Handle< FermState<T,Q,Q> > state_,
					 const SysSolverQUDAMULTIGRIDWilsonParams& invParam_) : 
      A(A_), invParam(invParam_)
    {
      QDPIO::cout << "LinOpSysSolverQUDAMULTIGRIDWilson:" << std::endl;

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
      mg_inv_param = newQudaInvertParam();
 

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
      const AnisoParam_t& aniso = invParam.WilsonParams.anisoParam;
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

      /****!!! FIXME: Before the final code remember to reset this to QUDA_CLOVER_WILSON_DSLASH */
      QDPIO::cout << "Remember for production to reset quda_inv_param.dslash_typeto QUDA_CLOVER_WILSON_DSLASH" << std::endl;
      quda_inv_param.dslash_type = QUDA_WILSON_DSLASH;
      mg_inv_param.dslash_type = QUDA_WILSON_DSLASH;

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
	QDPIO::cerr << "Unknown Solver type" << std::endl;
	QDP_abort(1);
	break;
      }
      //Params added for now to get this to initialize.
      mg_inv_param.inv_type = QUDA_GCR_INVERTER;
      mg_inv_param.tol = 1e-10;
      mg_inv_param.maxiter = 10000;
      mg_inv_param.reliable_delta = 1e-10;
      mg_inv_param.gcrNkrylov = 10;
      mg_inv_param.verbosity = QUDA_VERBOSE;
      mg_inv_param.verbosity_precondition = QUDA_VERBOSE;


      
      Real diag_mass;
      {
	// auto is C++11 so I don't have to remember all the silly typenames
	auto wlparams = invParam.WilsonParams;

	auto aniso = wlparams.anisoParam;
	
	Real ff = where(aniso.anisoP, aniso.nu / aniso.xi_0, Real(1));
	diag_mass = 1 + (Nd-1)*ff + wlparams.Mass;
      }


      quda_inv_param.kappa = static_cast<double>(1)/(static_cast<double>(2)*toDouble(diag_mass));
      /**** END FIXME XXX ***/
     
      quda_inv_param.tol = toDouble(invParam.RsdTarget);
      quda_inv_param.maxiter = invParam.MaxIter;
      quda_inv_param.reliable_delta = toDouble(invParam.Delta);

      // Solution type
      //quda_inv_param.solution_type = QUDA_MATPC_SOLUTION;
      //Taken from invert test.
      quda_inv_param.solution_type = QUDA_MATPC_SOLUTION;
      quda_inv_param.solve_type = QUDA_DIRECT_PC_SOLVE; 

#if 1
      if( invParam.asymmetricP ) { 
	QDPIO::cout << "Using Asymmetric Linop: A_oo - D A^{-1}_ee D" << std::endl;
	quda_inv_param.matpc_type = QUDA_MATPC_ODD_ODD_ASYMMETRIC;
      }
      else { 
	QDPIO::cout << "Using Symmetric Linop: 1 - A^{-1}_oo D A^{-1}_ee D" << std::endl;
	quda_inv_param.matpc_type = QUDA_MATPC_ODD_ODD;
      }
#endif

      quda_inv_param.dagger = QUDA_DAG_NO;
      quda_inv_param.mass_normalization = QUDA_KAPPA_NORMALIZATION;

      quda_inv_param.cpu_prec = cpu_prec;
      quda_inv_param.cuda_prec = gpu_prec;
      quda_inv_param.cuda_prec_sloppy = gpu_half_prec;
      //Add some lines for mg_inv_param.
      mg_inv_param.cpu_prec = cpu_prec;
      mg_inv_param.cuda_prec = gpu_prec;
      mg_inv_param.cuda_prec_sloppy = gpu_half_prec;

      //Clover stuff
      // This doesn't hurt to leave in, to make sure these are in a well defined state
      // But yucky!!!
      mg_inv_param.clover_cpu_prec = cpu_prec;
      mg_inv_param.clover_cuda_prec = gpu_prec;
      mg_inv_param.clover_cuda_prec_sloppy = gpu_half_prec;
      mg_inv_param.clover_cuda_prec_precondition = gpu_prec;
      mg_inv_param.clover_order = QUDA_PACKED_CLOVER_ORDER;
      //
      //Done...
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

      // Autotuning
      if( invParam.tuneDslashP ) { 
	QDPIO::cout << "Enabling Dslash Autotuning" << std::endl;

	quda_inv_param.tune = QUDA_TUNE_YES;
	mg_inv_param.tune = QUDA_TUNE_YES;
      }
      else { 
	QDPIO::cout << "Disabling Dslash Autotuning" << std::endl;
       
	quda_inv_param.tune = QUDA_TUNE_NO;
	mg_inv_param.tune = QUDA_TUNE_NO;
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

     if( invParam.MULTIGRIDParamsP ) {
	QDPIO::cout << "Setting MULTIGRID solver params" << std::endl;
	// Dereference handle
	MULTIGRIDSolverParams ip = *(invParam.MULTIGRIDParams);

	// Set preconditioner precision
	switch( ip.prec ) { 
	case HALF:
	  mg_inv_param.cuda_prec_precondition = QUDA_HALF_PRECISION;
	  quda_inv_param.cuda_prec_precondition = QUDA_HALF_PRECISION;
	  mg_inv_param.clover_cuda_prec_precondition = QUDA_HALF_PRECISION;
	  quda_inv_param.clover_cuda_prec_precondition = QUDA_HALF_PRECISION;
	  q_gauge_param.cuda_prec_precondition = QUDA_HALF_PRECISION;
	  break;

	case SINGLE:
	  mg_inv_param.cuda_prec_precondition = QUDA_SINGLE_PRECISION;
          mg_inv_param.clover_cuda_prec_precondition = QUDA_SINGLE_PRECISION;
          quda_inv_param.cuda_prec_precondition = QUDA_SINGLE_PRECISION;
          quda_inv_param.clover_cuda_prec_precondition = QUDA_SINGLE_PRECISION;
	  q_gauge_param.cuda_prec_precondition = QUDA_SINGLE_PRECISION;
	  break;

	case DOUBLE:
	  mg_inv_param.cuda_prec_precondition = QUDA_DOUBLE_PRECISION;
          mg_inv_param.clover_cuda_prec_precondition = QUDA_DOUBLE_PRECISION;
          quda_inv_param.cuda_prec_precondition = QUDA_DOUBLE_PRECISION;
          quda_inv_param.clover_cuda_prec_precondition = QUDA_DOUBLE_PRECISION;
	  q_gauge_param.cuda_prec_precondition = QUDA_DOUBLE_PRECISION;
	  break;
	default:
	  mg_inv_param.cuda_prec_precondition = QUDA_HALF_PRECISION;
          mg_inv_param.clover_cuda_prec_precondition = QUDA_HALF_PRECISION;
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
     
     for(int mu=0; mu < Nd; mu++) {
#ifndef BUILD_QUDA_DEVIFACE_GAUGE
        gauge[mu] = (void *)&(links_single[mu].elem(all.start()).elem().elem(0,0).real());
#else
       gauge[mu] = QDPCache::Instance().getDevicePtr( links_single[mu].getId() );
       QDPIO::cout << "MDAGM CUDA gauge[" << mu << "] in = " << gauge[mu] << "\n";
#endif
     }

    loadGaugeQuda((void *)gauge, &q_gauge_param);


            MULTIGRIDSolverParams ip = *(invParam.MULTIGRIDParams);
            //
	quda_inv_param.tol_precondition = toDouble(ip.tol);
	quda_inv_param.maxiter_precondition = ip.maxIterations;
	//quda_inv_param.gcrNkrylov = ip.gcrNkrylov;
	quda_inv_param.gcrNkrylov = 10;
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
	
	quda_inv_param.verbosity_precondition = QUDA_VERBOSE;
		
	//MG is the only option.
	quda_inv_param.inv_type_precondition = QUDA_MG_INVERTER;
	//New invert test changes here.
	
	mg_inv_param.sp_pad = 0;
  	mg_inv_param.cl_pad = 0;

  	mg_inv_param.preserve_source = QUDA_PRESERVE_SOURCE_NO;
  	mg_inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
  	mg_inv_param.dirac_order = QUDA_DIRAC_ORDER;

	mg_inv_param.input_location = QUDA_CPU_FIELD_LOCATION;
  	mg_inv_param.output_location = QUDA_CPU_FIELD_LOCATION;

	mg_inv_param.kappa = quda_inv_param.kappa;

	mg_inv_param.dagger = QUDA_DAG_NO;
  	mg_inv_param.mass_normalization = QUDA_KAPPA_NORMALIZATION;

  	mg_inv_param.matpc_type = QUDA_MATPC_ODD_ODD;
  	mg_inv_param.solution_type = QUDA_MAT_SOLUTION;
	mg_inv_param.solve_type = QUDA_DIRECT_SOLVE;

	mg_param.invert_param = &mg_inv_param;
	mg_param.n_level = ip.mg_levels;
	for (int i=0; i<mg_param.n_level; i++) {
    		for (int j=0; j<QUDA_MAX_DIM; j++) {
		  if( i < mg_param.n_level-1 ) {
      			mg_param.geo_block_size[i][j] = ip.blocking[i][j];
                  }
		  else {
		    mg_param.geo_block_size[i][j] = 4;
		  }
    		}
    		mg_param.spin_block_size[i] = 1;
    		mg_param.n_vec[i] = ip.nvec;
    		mg_param.nu_pre[i] = ip.nu_pre;
    		mg_param.nu_post[i] = ip.nu_post;

		switch( ip.smootherType ) { 
      		case MR: 
			mg_param.smoother[i] = QUDA_MR_INVERTER;
			mg_param.omega[i] = 0.85;
			break;
      		default:
			QDPIO::cout << "Unknown or no smother type specified, no smoothing inverter will be used." << std::endl;
			mg_param.smoother[i] = QUDA_INVALID_INVERTER;
			break;
      		}	
		mg_param.location[i] = QUDA_CUDA_FIELD_LOCATION;
		mg_param.smoother_solve_type[i] = QUDA_DIRECT_PC_SOLVE;
                mg_param.coarse_grid_solution_type[i] = QUDA_MATPC_SOLUTION;
    	}

  	// only coarsen the spin on the first restriction
  	mg_param.spin_block_size[0] = 2;

  	// coarse grid solver is GCR
  	mg_param.smoother[ip.mg_levels-1] = QUDA_GCR_INVERTER;

  	mg_param.compute_null_vector = ip.generate_nullspace ? QUDA_COMPUTE_NULL_VECTOR_YES
    	: QUDA_COMPUTE_NULL_VECTOR_NO;
	
	mg_param.vec_infile[0] = '\0';
        mg_param.vec_outfile[0] = '\0';

	QDPIO::cout<<"Basic MULTIGRID params copied."<<std::endl;
	quda_inv_param.verbosity = QUDA_VERBOSE;
	

   
     // setup the multigrid solver
     void *mg_preconditioner = newMultigridQuda(&mg_param);
     QDPIO::cout<<"NewMultigridQuda state initialized."<<std::endl;
     quda_inv_param.preconditioner = mg_preconditioner;
     QDPIO::cout<<"MULTIGRID preconditioner set."<<std::endl;
     //
   
      
    }
    

    //! Destructor is automatic
    ~LinOpSysSolverQUDAMULTIGRIDWilson() 
    {
      QDPIO::cout << "Destructing" << std::endl;
      freeGaugeQuda();

      /* FIXME: I Don't use it. Should I still do this? */
      // freeCloverQuda();
      destroyMultigridQuda(quda_inv_param.preconditioner);
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
	res = qudaInvert( g_chi,g_psi);      
	QDPIO::cout << "Untransforming solution." << std::endl;
	psi[ rb[1]]  = adj(GFixMat)*g_psi;

      }
      else { 
	res = qudaInvert(chi,psi);
      }      

      swatch.stop();
      double time = swatch.getTimeInSeconds();


      { 
	T r;
	r[A->subset()]=chi;
	T tmp;
	(*A)(tmp, psi, PLUS);
	r[A->subset()] -= tmp;
	res.resid = sqrt(norm2(r, A->subset()));
      }

      Double rel_resid = res.resid/sqrt(norm2(chi,A->subset()));

      QDPIO::cout << "QUDA_MULTIGRID_"<< solver_string <<"_WILSON_SOLVER: " << res.n_count << " iterations. Rsd = " << res.resid << " Relative Rsd = " << rel_resid << std::endl;
   
      // Convergence Check/Blow Up
      if ( ! invParam.SilentFailP ) { 
	      if (  toBool( rel_resid >  invParam.RsdToleranceFactor*invParam.RsdTarget) ) { 
        	QDPIO::cerr << "ERROR: QUDA MULTIGRID Solver residuum is outside tolerance: QUDA resid="<< rel_resid << " Desired =" << invParam.RsdTarget << " Max Tolerated = " << invParam.RsdToleranceFactor*invParam.RsdTarget << std::endl; 
        	QDP_abort(1);
      	      }
      }

      END_CODE();
      return res;
    }


  private:
    // Hide default constructor
    LinOpSysSolverQUDAMULTIGRIDWilson() {}
    
#if 1
    Q links_orig;
#endif

    U GFixMat;
    QudaPrecision_s cpu_prec;
    QudaPrecision_s gpu_prec;
    QudaPrecision_s gpu_half_prec;

    Handle< LinearOperator<T> > A;
    const SysSolverQUDAMULTIGRIDWilsonParams invParam;
    QudaGaugeParam q_gauge_param;
    QudaInvertParam quda_inv_param;
    QudaInvertParam mg_inv_param;	
    QudaMultigridParam mg_param;

    SystemSolverResults_t qudaInvert( const T& chi_s,
				      T& psi_s     
				      )const ;

    std::string solver_string;
  };


} // End namespace

#endif // BUILD_QUDA
#endif 


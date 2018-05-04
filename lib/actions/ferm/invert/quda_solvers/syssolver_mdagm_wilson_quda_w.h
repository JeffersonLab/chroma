// -*- C++ -*-
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by BiCGStab
 */

#ifndef __syssolver_mdagm_quda_wilson_h__
#define __syssolver_mdagm_quda_wilson_h__

#include "chroma_config.h"

#ifdef BUILD_QUDA

#include "handle.h"
#include "state.h"
#include "syssolver.h"
#include "linearop.h"
#include "lmdagm.h"
#include "actions/ferm/fermbcs/simple_fermbc.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/invert/quda_solvers/syssolver_quda_wilson_params.h"
#include "actions/ferm/invert/quda_solvers/quda_gcr_params.h"
#include "meas/gfix/temporal_gauge.h"
#include "io/aniso_io.h"
#include <string>

#include "util/gauge/reunit.h"

#include <quda.h>


namespace Chroma
{

  //! Richardson system solver namespace
  namespace MdagMSysSolverQUDAWilsonEnv
  {
    //! Register the syssolver
    bool registerAll();
  }



  //! Solve a Wilson Fermion System using the QUDA inverter
  /*! \ingroup invert
 *** WARNING THIS SOLVER WORKS FOR Wilson FERMIONS ONLY ***
   */
 
  class MdagMSysSolverQUDAWilson : public MdagMSystemSolver<LatticeFermion>
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
    MdagMSysSolverQUDAWilson(Handle< LinearOperator<T> > A_,
					 Handle< FermState<T,Q,Q> > state_,
					 const SysSolverQUDAWilsonParams& invParam_) : 
      A(A_), invParam(invParam_)
    {
      QDPIO::cout << "MdagMSysSolverQUDAWilson:" << std::endl;

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
      q_gauge_param.gauge_order = QUDA_QDP_GAUGE_ORDER; // gauge[mu], p

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
      quda_inv_param.dslash_type = QUDA_WILSON_DSLASH;

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


      Real massParam = Real(1) + Real(3)/Real(q_gauge_param.anisotropy) + invParam.WilsonParams.Mass;

      quda_inv_param.kappa = 1.0/(2*toDouble(massParam));

      // FIXME: If QUDA ever starts to compute our clover term we will need to fix this
      // Right now it is a dummy value, since we pass in the clover term
      quda_inv_param.clover_coeff = 1.0;

      quda_inv_param.mass_normalization = QUDA_ASYMMETRIC_MASS_NORMALIZATION;
      
      quda_inv_param.tol = toDouble(invParam.RsdTarget);
      quda_inv_param.maxiter = invParam.MaxIter;
      quda_inv_param.reliable_delta = toDouble(invParam.Delta);

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
      case MR:
	quda_inv_param.solve_type = QUDA_DIRECT_PC_SOLVE;
	break;

      default:
	quda_inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;   
	
	break;
      }


      quda_inv_param.matpc_type = QUDA_MATPC_ODD_ODD;


      quda_inv_param.dagger = QUDA_DAG_NO;


      quda_inv_param.cpu_prec = cpu_prec;
      quda_inv_param.cuda_prec = gpu_prec;
      quda_inv_param.cuda_prec_sloppy = gpu_half_prec;
      quda_inv_param.preserve_source = QUDA_PRESERVE_SOURCE_YES;
      quda_inv_param.use_init_guess = QUDA_USE_INIT_GUESS_NO;
      quda_inv_param.dirac_order = QUDA_DIRAC_ORDER;
      quda_inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;

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

	// Set preconditioner precision
	switch( ip.precPrecondition ) { 
	case HALF:
	  quda_inv_param.cuda_prec_precondition = QUDA_HALF_PRECISION;
	  q_gauge_param.cuda_prec_precondition = QUDA_HALF_PRECISION;
	  break;

	case SINGLE:
	  quda_inv_param.cuda_prec_precondition = QUDA_SINGLE_PRECISION;
	  q_gauge_param.cuda_prec_precondition = QUDA_SINGLE_PRECISION;
	  break;

	case DOUBLE:
	  quda_inv_param.cuda_prec_precondition = QUDA_DOUBLE_PRECISION;
	   q_gauge_param.cuda_prec_precondition = QUDA_DOUBLE_PRECISION;
	  break;
	default:
	  quda_inv_param.cuda_prec_precondition = QUDA_HALF_PRECISION;
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

      for(int mu=0; mu < Nd; mu++) { 
	gauge[mu] = (void *)&(links_single[mu].elem(all.start()).elem().elem(0,0).real());

      }
      loadGaugeQuda((void *)gauge, &q_gauge_param); 
      
      
    }
    

    //! Destructor is automatic
    ~MdagMSysSolverQUDAWilson() 
    {
      QDPIO::cout << "Destructing" << std::endl;
      freeGaugeQuda();
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

      START_CODE();
      StopWatch swatch;
      swatch.start();

      //    T MdagChi;

      // This is a CGNE. So create new RHS
      //      (*A)(MdagChi, chi, MINUS);
      // Handle< LinearOperator<T> > MM(new MdagMMdagM<T>(A));
      if ( invParam.axialGaugeP ) { 
	T g_chi,g_psi;

	// Gauge Fix source and initial guess
	QDPIO::cout << "Gauge Fixing source and initial guess" << std::endl;
        g_chi[ rb[1] ]  = GFixMat * chi;
	g_psi[ rb[1] ]  = GFixMat * psi;
	QDPIO::cout << "Solving" << std::endl;
	res = qudaInvert(g_chi,
			 g_psi);      
	QDPIO::cout << "Untransforming solution." << std::endl;
	psi[ rb[1]]  = adj(GFixMat)*g_psi;

      }
      else { 
	QDPIO::cout << "Calling QUDA Invert" << std::endl;
	res = qudaInvert(chi,
			 psi);      
      }      

      swatch.stop();
      double time = swatch.getTimeInSeconds();


      { 
	T r;
	r[A->subset()]=chi;
	T tmp,tmp2;
	(*A)(tmp, psi, PLUS);
	(*A)(tmp2, tmp, MINUS);
	r[A->subset()] -= tmp2;
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


   SystemSolverResults_t operator() (T& psi, const T& chi, Chroma::AbsChronologicalPredictor4D<T>& predictor ) const
    {
      SystemSolverResults_t res;

      START_CODE();
      StopWatch swatch;
      swatch.start();
      {
	Handle< LinearOperator<T> > MdagM( new MdagMLinOp<T>(A) );
	predictor(psi, (*MdagM), chi);
      }
      res = (*this)(psi, chi);
      predictor.newVector(psi);
      swatch.stop();
      double time = swatch.getTimeInSeconds();
      QDPIO::cout << "QUDA_"<< solver_string <<"_CLOVER_SOLVER: Total time (with prediction)=" << time << std::endl;
      END_CODE();
      return res;
    }


  private:
    // Hide default constructor
    MdagMSysSolverQUDAWilson() {}
    
#if 1
    Q links_orig;
#endif

    U GFixMat;
    QudaPrecision_s cpu_prec;
    QudaPrecision_s gpu_prec;
    QudaPrecision_s gpu_half_prec;

    Handle< LinearOperator<T> > A;
    const SysSolverQUDAWilsonParams invParam;
    QudaGaugeParam q_gauge_param;
    QudaInvertParam quda_inv_param;

    SystemSolverResults_t qudaInvert(const T& chi_s,
				     T& psi_s     
				     )const ;

    std::string solver_string;
  };


} // End namespace

#endif // BUILD_QUDA
#endif 


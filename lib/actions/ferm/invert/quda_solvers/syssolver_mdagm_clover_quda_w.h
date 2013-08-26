// -*- C++ -*-
// $Id: syssolver_linop_quda_clover.h,v 1.9 2009-10-16 13:37:39 bjoo Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by BiCGStab
 */

#ifndef __syssolver_mdagm_quda_clover_w_h__
#define __syssolver_mdagm_quda_clover_w_h__

#include "chroma_config.h"


#ifdef BUILD_QUDA

#include "handle.h"
#include "state.h"
#include "syssolver.h"
#include "linearop.h"
#include "lmdagm.h"
#include "actions/ferm/fermbcs/simple_fermbc.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/invert/quda_solvers/syssolver_quda_clover_params.h"
#include "actions/ferm/linop/clover_term_w.h"
#include "meas/gfix/temporal_gauge.h"
#include "io/aniso_io.h"
#include <string>
#include "update/molecdyn/predictor/null_predictor.h"

#include "util/gauge/reunit.h"

#include <quda.h>

//#undef BUILD_QUDA_DEVIFACE_GAUGE
//#undef BUILD_QUDA_DEVIFACE_SPINOR
//#undef BUILD_QUDA_DEVIFACE_CLOVER

//#include <util_quda.h>
using namespace std;

namespace Chroma
{

  //! Richardson system solver namespace
  namespace MdagMSysSolverQUDACloverEnv
  {
    //! Register the syssolver
    bool registerAll();
  }



  //! Solve a Clover Fermion System using the QUDA inverter
  /*! \ingroup invert
 *** WARNING THIS SOLVER WORKS FOR Clover FERMIONS ONLY ***
   */
 
  class MdagMSysSolverQUDAClover : public MdagMSystemSolver<LatticeFermion>
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
    MdagMSysSolverQUDAClover(Handle< LinearOperator<T> > A_,
					 Handle< FermState<T,Q,Q> > state_,
					 const SysSolverQUDACloverParams& invParam_) : 
      A(A_), state(state_), invParam(invParam_), clov(new CloverTermT<T, U>::Type_t()), invclov(new CloverTermT<T, U>::Type_t())
    {
      QDPIO::cout << "MdagMSysSolverQUDAClover:" << endl;

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
      QDPIO::cout << "MDAGM Using QDP-JIT gauge order" << endl;
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
	QDPIO::cout << "Fixing Temporal Gauge" << endl;
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
     quda_inv_param.inv_type = QUDA_CG_INVERTER;   
     solver_string = "CG";
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

      if( ! invParam.asymmetricP ) { 
	QDPIO::cout << "For MdagM we can only use asymmetric Linop: A_oo - D A^{-1}_ee D, overriding your choice" << endl;
	
      }
      // Only support Asymmetric linop
      quda_inv_param.matpc_type = QUDA_MATPC_ODD_ODD_ASYMMETRIC;


      quda_inv_param.dagger = QUDA_DAG_NO;
      quda_inv_param.mass_normalization = QUDA_KAPPA_NORMALIZATION;

      quda_inv_param.cpu_prec = cpu_prec;
      quda_inv_param.cuda_prec = gpu_prec;
      quda_inv_param.cuda_prec_sloppy = gpu_half_prec;
      quda_inv_param.preserve_source = QUDA_PRESERVE_SOURCE_NO;
      quda_inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;

#ifndef BUILD_QUDA_DEVIFACE_SPINOR
      quda_inv_param.dirac_order = QUDA_DIRAC_ORDER;
#else
      QDPIO::cout << "MDAGM Using QDP-JIT spinor order" << endl;
      quda_inv_param.dirac_order    = QUDA_QDPJIT_DIRAC_ORDER;
      quda_inv_param.input_location = QUDA_CUDA_FIELD_LOCATION;
      quda_inv_param.output_location = QUDA_CUDA_FIELD_LOCATION;
#endif


      // Autotuning
      if( invParam.tuneDslashP ) { 
	QDPIO::cout << "Enabling Dslash Autotuning" << endl;

	quda_inv_param.tune = QUDA_TUNE_YES;
      }
      else { 
	QDPIO::cout << "Disabling Dslash Autotuning" << endl;
       
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
	QDPIO::cout << "Setting inner solver params" << endl;
	// Dereference handle
	GCRInnerSolverParams ip = *(invParam.innerParams);

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
	case MR:
	  quda_inv_param.inv_type_precondition= QUDA_MR_INVERTER;
	  break;
	  
	default:
	  quda_inv_param.inv_type_precondition = QUDA_MR_INVERTER;   
	  break;
	}
      }
      else { 
	QDPIO::cout << "Setting Precondition stuff to defaults for not using" << endl;
	quda_inv_param.inv_type_precondition= QUDA_INVALID_INVERTER;
	quda_inv_param.tol_precondition = 1.0e-1;
	quda_inv_param.maxiter_precondition = 1000;
	quda_inv_param.verbosity_precondition = QUDA_SILENT;
        quda_inv_param.gcrNkrylov = 1;
      }
      
      // Clover precision and order
      quda_inv_param.clover_cpu_prec = cpu_prec;
      quda_inv_param.clover_cuda_prec = gpu_prec;
      quda_inv_param.clover_cuda_prec_sloppy = gpu_half_prec;

#ifndef BUILD_QUDA_DEVIFACE_CLOVER
      quda_inv_param.clover_order = QUDA_PACKED_CLOVER_ORDER;
#else      
      QDPIO::cout << "MDAGM Clover CUDA location\n";
      quda_inv_param.clover_location = QUDA_CUDA_FIELD_LOCATION;
      quda_inv_param.clover_order = QUDA_QDPJIT_CLOVER_ORDER;
#endif

      if( invParam.verboseP ) { 
	quda_inv_param.verbosity = QUDA_VERBOSE;
      }
      else { 
	quda_inv_param.verbosity = QUDA_SUMMARIZE;
      }
      
      // Set up the links     
      void* gauge[4]; 

      for(int mu=0; mu < Nd; mu++) { 
#ifndef BUILD_QUDA_DEVIFACE_GAUGE
	gauge[mu] = (void *)&(links_single[mu].elem(all.start()).elem().elem(0,0).real());
#else
	gauge[mu] = QDPCache::Instance().getDevicePtr( links_single[mu].getId() );
	//std::cout << "MDAGM CUDA gauge[" << mu << "] in = " << gauge[mu] << "\n";
#endif
      }


      loadGaugeQuda((void *)gauge, &q_gauge_param); 
      
      // Setup Clover Term
      QDPIO::cout << "Creating CloverTerm" << endl;
      clov->create(fstate, invParam_.CloverParams);
      // Don't recompute, just copy
      invclov->create(fstate, invParam_.CloverParams);
      
      QDPIO::cout << "Inverting CloverTerm" << endl;
      invclov->choles(0);
      invclov->choles(1);

#ifndef BUILD_QUDA_DEVIFACE_CLOVER
      multi1d<QUDAPackedClovSite<REALT> > packed_clov;

      // Only compute clover if we're using asymmetric preconditioner
      packed_clov.resize(all.siteTable().size());

      clov->packForQUDA(packed_clov, 0);
      clov->packForQUDA(packed_clov, 1);

      // Always need inverse
      multi1d<QUDAPackedClovSite<REALT> > packed_invclov(all.siteTable().size());
      invclov->packForQUDA(packed_invclov, 0);
      invclov->packForQUDA(packed_invclov, 1);

      loadCloverQuda(&(packed_clov[0]), &(packed_invclov[0]),&quda_inv_param);
#else
      void *clover[2];
      void *cloverInv[2];

      clover[0] = QDPCache::Instance().getDevicePtr( clov->getOffId() );
      clover[1] = QDPCache::Instance().getDevicePtr( clov->getDiaId() );

      cloverInv[0] = QDPCache::Instance().getDevicePtr( invclov->getOffId() );
      cloverInv[1] = QDPCache::Instance().getDevicePtr( invclov->getDiaId() );

      // std::cout << "MDAGM clover CUDA pointers: " 
      // 		<< clover[0] << " "
      // 		<< clover[1] << " "
      // 		<< cloverInv[0] << " "
      // 		<< cloverInv[1] << "\n";

      loadCloverQuda( (void*)(clover) , (void*)(cloverInv) ,&quda_inv_param);
#endif
      
    }
    

    //! Destructor is automatic
    ~MdagMSysSolverQUDAClover() 
    {
      QDPIO::cout << "Destructing" << endl;
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
    SystemSolverResults_t operator() (T& psi, const T& chi ) const
    {
      SystemSolverResults_t res;
      Null4DChronoPredictor not_predicting;
      (*this)(psi,chi, not_predicting);
      START_CODE();
      END_CODE();
      return res;
    }


   SystemSolverResults_t operator() (T& psi, const T& chi, Chroma::AbsChronologicalPredictor4D<T>& predictor ) const
    {
      SystemSolverResults_t res;
      SystemSolverResults_t res1;
      SystemSolverResults_t res2;

      START_CODE();
      // This is a solve with initial guess. So reset the policy  (quda_inv_param is mutable)
      QudaUseInitGuess old_guess_policy = quda_inv_param.use_init_guess;
      quda_inv_param.use_init_guess = QUDA_USE_INIT_GUESS_YES;


      StopWatch swatch;
      swatch.start();

      // Determine if 2 step solve is needed
      bool two_step=false;
      switch(  invParam.solverType  ) { 
      case BICGSTAB:
	two_step = true;
	break;
      case GCR:
	two_step = true;
	break;
      case MR:
	two_step = true;
	break;
      case CG:
	two_step = false;
	break;
      default:
	two_step = false;
	break;
      };

      // This is a solver whic will use an initial guess:

      // Create MdagM op
      Handle< LinearOperator<T> > MdagM( new MdagMLinOp<T>(A) );

      bool failP = false;

      if ( ! two_step ) { 

	// Single Step Solve
	QDPIO::cout << "Single Step Solve" << endl;
	predictor(psi, (*MdagM), chi);
	res = qudaInvert(*clov,
			 *invclov,
			 chi,
			 psi);      
	predictor.newVector(psi);

	
      }
      else { 

	// TWO STEP SOLVE
	try {
	  QDPIO::cout << "Two Step Solve" << endl;

	  T Y;
	  Y[ A->subset() ] = psi; // Y is initial guess

	  // Try to cast the predictor to a two step predictor 
	  AbsTwoStepChronologicalPredictor4D<T>& two_step_predictor = 
	     dynamic_cast<AbsTwoStepChronologicalPredictor4D<T>& >(predictor);

	  // Predict Y
	  two_step_predictor.predictY(Y,*A,chi);

	  // Now want to solve with M^\dagger on this.
	  quda_inv_param.solution_type = QUDA_MATPC_SOLUTION;
	  quda_inv_param.dagger = QUDA_DAG_YES;
	  res1 = qudaInvert(*clov,
			    *invclov,
			    chi,
			    Y);  
	  two_step_predictor.newYVector(Y);

	  // Step 2: Solve M X = Y
	  // Predict X 
	  two_step_predictor.predictX(psi,*MdagM, chi);
	  quda_inv_param.dagger = QUDA_DAG_NO; // Solve without dagger
	  res2 = qudaInvert(*clov,
			    *invclov,
			    Y,
			    psi);  
	  two_step_predictor.newXVector(psi);
	  
	}
	catch( std::bad_cast ) { 
	  // Failed to cast the predictor to a two step predictor
	  // In this case we do not predict Y directly, but 
	  // try to predict  X ~ (MdagM)^{-1} chi => X ~ M^{-1} M^\dag{-1} chi
	  // we get Y from M X ~ M^\dag{-1} chi 
	  // 
	  T Y;
	  Y[ A->subset() ] = psi;  // Set Y to what the user first gave (in case of NULL predictor)
	  predictor(psi, (*MdagM), chi); // We can only predict X
	  (*A)(Y, psi, PLUS); // and then Y = M X is a good guess for Y
	
	  quda_inv_param.solution_type = QUDA_MATPC_SOLUTION;
	  quda_inv_param.dagger = QUDA_DAG_YES;
	  res1 = qudaInvert(*clov,
			    *invclov,
			    chi,
			    Y);  

	  // Step 2: Solve M X = Y
	  // Predict X 
	  quda_inv_param.dagger = QUDA_DAG_NO; // Solve without dagger
	  res2 = qudaInvert(*clov,
			    *invclov,
			    Y,
			    psi);  
	  predictor.newVector(psi);
	}

	// reset params to their default value
	quda_inv_param.solution_type = QUDA_MATPCDAG_MATPC_SOLUTION;
	quda_inv_param.dagger = QUDA_DAG_NO; // Solve without dagger
	res.n_count = res1.n_count + res2.n_count;  // Two step solve so combine iteration count
      }
      swatch.stop();
      double time = swatch.getTimeInSeconds();

      // reset init guess policy
      quda_inv_param.use_init_guess = old_guess_policy;

      // Check solution 
      { 
	T r;
	r[A->subset()]=chi;
	T tmp;
	(*MdagM)(tmp, psi, PLUS);
	r[A->subset()] -= tmp;
	res.resid = sqrt(norm2(r, A->subset()));
      }

      Double rel_resid = res.resid/sqrt(norm2(chi,A->subset()));

      QDPIO::cout << "QUDA_"<< solver_string <<"_CLOVER_SOLVER: " << res.n_count << " iterations. Rsd = " << res.resid << " Relative Rsd = " << rel_resid << endl;
      QDPIO::cout << "QUDA_"<< solver_string <<"_CLOVER_SOLVER: Total time (with prediction)=" << time << endl;


      // Check for failure
      if (  toBool( rel_resid >  invParam.RsdToleranceFactor*invParam.RsdTarget) ) { 
	if ( ! invParam.SilentFailP ) { 

	  // If we are here we are meant to be verbose about it. 
	  QDPIO::cerr << "ERROR: QUDA Solver residuum is outside tolerance: QUDA resid="<< rel_resid << " Desired =" << invParam.RsdTarget << " Max Tolerated = " << invParam.RsdToleranceFactor*invParam.RsdTarget << endl; 

	  // Check if we need to dump...
	  if ( invParam.dump_on_failP ) { 

	    // Dump config
	    int my_random_int = rand();
	    XMLBufferWriter rhs_filebuf;
	    XMLBufferWriter rhs_recbuf;
	    XMLBufferWriter gauge_filebuf;
	    XMLBufferWriter gauge_recbuf;

	    // Complete nonsense for the metadata
	    int foo=5;
	    push(rhs_filebuf, "RHSFile");
	    write(rhs_filebuf, "dummyInt", foo++);
	    pop(rhs_filebuf);

	    push(rhs_recbuf, "RHSRecord");
	    write(rhs_recbuf, "dummyInt", foo++);
	    pop(rhs_recbuf);

	    push(gauge_filebuf, "GAUGEFile");
	    write(gauge_filebuf, "dummyInt", foo++);
	    pop(gauge_filebuf);

	    push(gauge_recbuf, "GAUGERecord");
	    write(gauge_recbuf, "dummyInt", foo++);
	    pop(gauge_recbuf);

	    // Filenames
	    std::ostringstream rhs_filename;
	    rhs_filename << "./rhs_diagnostic_" << my_random_int <<".lime";

	    std::ostringstream gauge_filename;
	    gauge_filename << "./gauge_diagnostic_" << my_random_int <<".lime";

	    // Hardwired QDPIO_PARALLEL now, as it seems to work even in a single file write :)
	    QDPFileWriter rhs_writer(rhs_filebuf, rhs_filename.str(), QDPIO_SINGLEFILE, QDPIO_PARALLEL);
	    QDPFileWriter gauge_writer(gauge_filebuf,gauge_filename.str(), QDPIO_SINGLEFILE, QDPIO_PARALLEL);

	    // DUMP RHS (will anyone ever read it?)
	    QDPIO::cout << "DUMPING RHS to " << rhs_filename.str() << endl;
	    write(rhs_writer, rhs_recbuf, chi);

	    // DUMP Gauge (fill anyone ever read it?)
	    QDPIO::cout << "DUMPING GAUGE FIELD to " << gauge_filename.str() << endl;
	    write(gauge_writer, gauge_recbuf, state->getLinks());
	    
	    rhs_writer.close();
	    gauge_writer.close();
	  }

	  if( invParam.backup_invP) { 
	    // Create the Backup Solver...
	    QDPIO::cout << "CREATING BACKUP SOLVER" << endl;
	    std::istringstream is( invParam.backup_inv_param.xml );
	    XMLReader backup_solver_reader( is );
	    
	    Handle< MdagMSystemSolver< T > > backup_solver;

	    try { 
	      backup_solver = TheMdagMFermSystemSolverFactory::Instance().createObject(
										       invParam.backup_inv_param.id,  
										       backup_solver_reader,
										       "/BackupSolverParam", 
										       state,
										       A );

	    }
	    catch(const std::string e) {
	      backup_solver_reader.print(std::cout);
	      QDPIO::cout << "Caught exception: " << e << endl;
	      QDP_abort(1);
	    }

	    QDPIO::cout << "PERFORMING BACKUP SOLVE" << endl;

	    // Backup solver will re-predict
	    res =  (*backup_solver)(psi, chi, predictor);

	    // Check solution 
	    { 
	      T r;
	      r[A->subset()]=chi;
	      T tmp;
	      (*MdagM)(tmp, psi, PLUS);
	      r[A->subset()] -= tmp;
	      res.resid = sqrt(norm2(r, A->subset()));
	    }

	    rel_resid = res.resid/sqrt(norm2(chi,A->subset()));
	    if (  toBool( rel_resid >  invParam.RsdToleranceFactor*invParam.RsdTarget) ) { 
	      QDPIO::cout << "ERROR: BACKUP SOLVE FAILED" << endl;
	      QDP_abort(1);
	    }
	  }
	  else { 
	    // No backup solve 
	    QDPIO::cout << "SOLVER FAILED: Aborting" << endl;
	    QDP_abort(1);
	  }
	}
	else { 
	  // (not so)SILENT FAILURE:
	  QDPIO::cout << "WARNING: Solver failed but SILENT FAILURE is ENABLED. Continuing" << endl;
	}

      }

      END_CODE();
      return res;
    }


  private:
    // Hide default constructor
    MdagMSysSolverQUDAClover() {}
    
    Q links_orig;

    U GFixMat;
    QudaPrecision_s cpu_prec;
    QudaPrecision_s gpu_prec;
    QudaPrecision_s gpu_half_prec;

    Handle< LinearOperator<T> > A;
    Handle< FermState<T,Q,Q> > state;
    const SysSolverQUDACloverParams invParam;
    QudaGaugeParam q_gauge_param;
    mutable QudaInvertParam quda_inv_param;

    Handle< typename CloverTermT<T, U>::Type_t > clov;
    Handle< typename CloverTermT<T, U>::Type_t > invclov;

    SystemSolverResults_t qudaInvert(const CloverTermT<T, U>::Type_t& clover,
				     const CloverTermT<T, U>::Type_t& inv_clov,
				     const T& chi_s,
				     T& psi_s     
				     )const ;

    std::string solver_string;
  };


} // End namespace

#endif // BUILD_QUDA
#endif 


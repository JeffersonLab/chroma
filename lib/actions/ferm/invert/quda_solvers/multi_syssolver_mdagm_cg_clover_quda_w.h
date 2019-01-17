// -*- C++ -*-
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by CG2 using CG
 */

#ifndef __multi_syssolver_mdagm_cg_clover_quda_w_h__
#define __multi_syssolver_mdagm_cg_clover_quda_w_h__

#include <cfloat>
#include <cstdio>
#include "chroma_config.h"

#ifdef BUILD_QUDA

#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "lmdagm.h"
#include "actions/ferm/fermbcs/simple_fermbc.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/invert/quda_solvers/multi_syssolver_quda_clover_params.h"

#include "actions/ferm/linop/clover_term_w.h"
#include "meas/gfix/temporal_gauge.h"
#include "io/aniso_io.h"
#include <string>

#include "util/gauge/reunit.h"

#include <quda.h>



#ifdef QDP_IS_QDPJIT
#include "actions/ferm/invert/quda_solvers/qdpjit_memory_wrapper.h"
#endif

namespace Chroma
{

  //! CG system solver namespace
  namespace MdagMMultiSysSolverCGQudaCloverEnv
  {
    //! Register the syssolver
    bool registerAll();
  }

  //! Solve a CG2 system. Here, the operator is NOT assumed to be hermitian
  /*! \ingroup invert
   */
  class MdagMMultiSysSolverCGQudaClover : public MdagMMultiSystemSolver<LatticeFermion>
  {
  public:
    typedef LatticeFermion T;
    typedef LatticeColorMatrix U;
    typedef multi1d<LatticeColorMatrix> Q;
    typedef multi1d<LatticeColorMatrix> P;
 
    typedef LatticeFermionF TF;
    typedef LatticeColorMatrixF UF;
    typedef multi1d<LatticeColorMatrixF> QF;
    typedef multi1d<LatticeColorMatrixF> PF;

    typedef LatticeFermionD TD;
    typedef LatticeColorMatrixD UD;   
    typedef multi1d<LatticeColorMatrixD> QD;
    typedef multi1d<LatticeColorMatrixD> PD;

    typedef WordType<T>::Type_t REALT;
    //! Constructor
    /*!
     * \param M_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    MdagMMultiSysSolverCGQudaClover(Handle< LinearOperator<T> > M_,
				      Handle< FermState<T,P,Q> > state_,
				      const MultiSysSolverQUDACloverParams& invParam_) : 
      A(M_), invParam(invParam_), clov(new CloverTermT<T,U>()), invclov(new CloverTermT<T,U>())

    {
      QDPIO::cout << "MdagMMultiSysSolverCGQUDAClover: " << std::endl;
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

      // Work out GPU Sloppy precision
      // Default: No Sloppy
      switch( invParam.cudaRefinementPrecision ) {
      case HALF:
    	  gpu_ref_prec = QUDA_HALF_PRECISION;
    	  break;
      case SINGLE:
    	  gpu_ref_prec = QUDA_SINGLE_PRECISION;
    	  break;
      case DOUBLE:
    	  gpu_ref_prec = QUDA_DOUBLE_PRECISION;
    	  break;
      default:
    	  gpu_ref_prec = gpu_prec;
    	  break;
      }

      // 2) pull 'new; GAUGE and Invert params
      // 
      QDPIO::cout << " Calling new QUDA Invert Param" << std::endl;
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
      //QDPIO::cout << "MULTI MDAGM Using QDP-JIT gauge order" << std::endl;
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
      q_gauge_param.cuda_prec_precondition = gpu_half_prec;
      q_gauge_param.cuda_prec_refinement_sloppy = gpu_ref_prec;

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

      // Default
      q_gauge_param.reconstruct_precondition = q_gauge_param.reconstruct_sloppy;

      //  Mathias's Hardwired version:
      //    q_gauge_param.reconstruct_refinement_sloppy = q_gauge_param.reconstruct_sloppy;

      // Parameter file based version
      switch( invParam.cudaRefinementReconstruct ) {
      case RECONS_NONE:
    	  q_gauge_param.reconstruct_refinement_sloppy = QUDA_RECONSTRUCT_NO;
    	  break;
      case RECONS_8:
    	  q_gauge_param.reconstruct_refinement_sloppy = QUDA_RECONSTRUCT_8;
    	  break;
      case RECONS_12:
    	  q_gauge_param.reconstruct_refinement_sloppy = QUDA_RECONSTRUCT_12;
    	  break;
      default:
    	  q_gauge_param.reconstruct_refinement_sloppy = QUDA_RECONSTRUCT_12;
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
      solver_string = "MULTI_CG";
      quda_inv_param.inv_type = QUDA_CG_INVERTER;
      quda_inv_param.use_alternative_reliable = 1;
      // Mass

      // Fiendish idea from Ron. Set the kappa=1/2 and use 
      // unmodified clover term, and ask for Kappa normalization
      // This should give us A - (1/2)D as the unpreconditioned operator
      // and probabl 1 - {1/4} A^{-1} D A^{-1} D as the preconditioned
      // op. Apart from the A_oo stuff on the antisymmetric we have
      // nothing to do...
      if( invParam.asymmetricP ) {
	quda_inv_param.mass_normalization = QUDA_KAPPA_NORMALIZATION;
	quda_inv_param.kappa = 0.5;
      } 
      else { 
      	quda_inv_param.mass_normalization = QUDA_KAPPA_NORMALIZATION;
      	quda_inv_param.kappa = 0.5;
      } 

      // FIXME: If we ever use QUDA to compute our clover term we need to fix this
      // and QUDA will have to deal with anisotropy. Right now this is a hack to make
      // the check_param pass. We pass in our clover term so this outghtnt matter
      quda_inv_param.clover_coeff = 1.0; // dummy vaue

      quda_inv_param.tol = toDouble(invParam.RsdTarget[0]);
      quda_inv_param.maxiter = invParam.MaxIter;
      quda_inv_param.reliable_delta = toDouble(invParam.Delta);
      quda_inv_param.pipeline = invParam.Pipeline;

      // Solution type
      quda_inv_param.solution_type = QUDA_MATPCDAG_MATPC_SOLUTION;

      // Solve type
      switch( invParam.solverType ) { 
      case CG: 
	quda_inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;
	break;
      default:
	QDPIO::cerr << "Only CG Is currently implemented for multi-shift" << std::endl;
	QDP_abort(1);
	
	break;
      }

      // Only suppfkort Asymmetric linop
      if ( invParam.asymmetricP ) { 
	QDPIO::cout << "Working with asymmetric solution" << std::endl;
      	quda_inv_param.matpc_type = QUDA_MATPC_ODD_ODD_ASYMMETRIC;
      }
      else {
	QDPIO::cout << "Working with symmetric solution" << std::endl; 
	quda_inv_param.matpc_type = QUDA_MATPC_ODD_ODD;
      }

      quda_inv_param.dagger = QUDA_DAG_NO;

      quda_inv_param.cpu_prec = cpu_prec;
      quda_inv_param.cuda_prec = gpu_prec;
      quda_inv_param.cuda_prec_sloppy = gpu_half_prec;
      quda_inv_param.cuda_prec_precondition = gpu_half_prec;
      quda_inv_param.cuda_prec_refinement_sloppy = gpu_ref_prec;
      quda_inv_param.preserve_source = QUDA_PRESERVE_SOURCE_NO;
      quda_inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;

#ifndef BUILD_QUDA_DEVIFACE_SPINOR
      quda_inv_param.dirac_order = QUDA_DIRAC_ORDER;
#else
      //QDPIO::cout << "MULTI MDAGM Using QDP-JIT spinor order" << std::endl;
      quda_inv_param.dirac_order    = QUDA_QDPJIT_DIRAC_ORDER;
      quda_inv_param.input_location = QUDA_CUDA_FIELD_LOCATION;
      quda_inv_param.output_location = QUDA_CUDA_FIELD_LOCATION;
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


      // PADDING
      
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


      // Setting GCR Preconditioner to defaults, as we don't use it..
      // This is kinda yucky.

      QDPIO::cout << "Setting Precondition stuff to defaults for not using" << std::endl;
      quda_inv_param.inv_type_precondition= QUDA_INVALID_INVERTER;
      quda_inv_param.tol_precondition = 1.0e-1;
      quda_inv_param.maxiter_precondition = 1000;
      quda_inv_param.verbosity_precondition = QUDA_SILENT;
      quda_inv_param.gcrNkrylov = 1;

      // Clover precision and order
      quda_inv_param.clover_cpu_prec = cpu_prec;
      quda_inv_param.clover_cuda_prec = gpu_prec;
      quda_inv_param.clover_cuda_prec_sloppy = gpu_half_prec;
      quda_inv_param.clover_cuda_prec_precondition = gpu_half_prec;
      quda_inv_param.clover_cuda_prec_refinement_sloppy = gpu_ref_prec;

#ifndef BUILD_QUDA_DEVIFACE_CLOVER
      quda_inv_param.clover_order = QUDA_PACKED_CLOVER_ORDER;
#else      
      QDPIO::cout << "MULTI MDAGM clover CUDA location\n";
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

#ifndef BUILD_QUDA_DEVIFACE_GAUGE
      for(int mu=0; mu < Nd; mu++) { 
	gauge[mu] = (void *)&(links_single[mu].elem(all.start()).elem().elem(0,0).real());
      }
#else
      GetMemoryPtrGauge(gauge,links_single);
      //gauge[mu] = GetMemoryPtr( links_single[mu].getId() );
#endif


      loadGaugeQuda((void *)gauge, &q_gauge_param); 

      //      Setup the clover term...
      QDPIO::cout << "Creating CloverTerm" << std::endl;
      clov->create(fstate, invParam_.CloverParams);
      invclov->create(fstate, invParam_.CloverParams);
      
      QDPIO::cout << "Inverting CloverTerm" << std::endl;
      invclov->choles(0);
      invclov->choles(1);


#ifndef BUILD_QUDA_DEVIFACE_CLOVER
      multi1d<QUDAPackedClovSite<REALT> > packed_clov;

      packed_clov.resize(all.siteTable().size());

      clov->packForQUDA(packed_clov, 0);
      clov->packForQUDA(packed_clov, 1);
 
      // Always need inverse
      multi1d<QUDAPackedClovSite<REALT> > packed_invclov(all.siteTable().size());
      invclov->packForQUDA(packed_invclov, 0);
      invclov->packForQUDA(packed_invclov, 1);

      loadCloverQuda(&(packed_clov[0]), &(packed_invclov[0]),&quda_inv_param);
#else

	quda_inv_param.clover_location = QUDA_CUDA_FIELD_LOCATION;
	quda_inv_param.clover_order = QUDA_QDPJIT_CLOVER_ORDER;

      void *clover[2];
      void *cloverInv[2];

      GetMemoryPtrClover(clov->getOffId(),clov->getDiaId(),invclov->getOffId(),invclov->getDiaId());

   
      loadCloverQuda( (void*)(clover) , (void*)(cloverInv) ,&quda_inv_param);
#endif

    }

    //! Destructor is automatic
    ~MdagMMultiSysSolverCGQudaClover() {
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
    SystemSolverResults_t operator() (multi1d<T>& psi, const multi1d<Real>& shifts, const T& chi) const
    {
      START_CODE();
      StopWatch swatch;
      swatch.reset();
      swatch.start();
      SystemSolverResults_t res;
      res.n_count = 0; 

      if( invParam.RsdTarget.size() > 1 ) {
	if( invParam.RsdTarget.size() != shifts.size()) { 
	  QDPIO::cerr << "There are: " << invParam.RsdTarget.size() << " values for RsdTarget but " << shifts.size() << " shifts" << std::endl;
	  QDPIO::cerr << "This doesnt match. Aborting" << std::endl;
	}
      }
      
      if ( invParam.axialGaugeP ) { 
	T g_chi;
	multi1d<T> g_psi(psi.size());

	// Gauge Fix source and initial guess
	QDPIO::cout << "Gauge Fixing source and initial guess" << std::endl;
        g_chi[ rb[1] ]  = GFixMat * chi;
	for(int s=0; s < psi.size(); s++) {
	  g_psi[s][ rb[1] ]  = zero; // All initial guesses are zero
	}

	QDPIO::cout << "Solving" << std::endl;
	res = qudaInvertMulti(
			 g_chi,
			 g_psi,
			 shifts);      
	QDPIO::cout << "Untransforming solution." << std::endl;
	for(int s=0; s< psi.size(); s++) { 
	  psi[s][ rb[1]]  = adj(GFixMat)*g_psi[s];
	}

      }
      else { 
	res = qudaInvertMulti(chi,
			 psi,
			 shifts);      
      }      

      swatch.stop();
      double time = swatch.getTimeInSeconds();

      bool abortFound = false;
      if (invParam.checkShiftsP )  {
        Double chinorm=norm2(chi, A->subset());
        multi1d<Double> r_rel(shifts.size());

#ifdef QUDA_DEBUG
        for(int i=0; i < shifts.size(); i++) {
	  char normpsi_subset[256];
	  char normpsi_full[256];
	  std::sprintf( normpsi_subset, "%.*e", DECIMAL_DIG, toDouble(norm2(psi[i],A->subset())) );
          std::sprintf( normpsi_full, "%.*e", DECIMAL_DIG, toDouble(norm2(psi[i])) );
	  QDPIO::cout << "psi[ " << i << " ] : norm( A->subset() ) = " << normpsi_subset << " norm(total) = " << normpsi_full << std::endl;
        }
#endif
        for(int i=0; i < shifts.size(); i++) { 
          T tmp1,tmp2;
          (*A)(tmp1, psi[i], PLUS);
          (*A)(tmp2, tmp1, MINUS);  // tmp2 = A^\dagger A psi
          tmp2[ A->subset() ] +=  shifts[i]* psi[i]; // tmp2 = ( A^\dagger A + shift_i ) psi
          T r;
          r[ A->subset() ] = chi - tmp2;
          r_rel[i] = sqrt(norm2(r, A->subset())/chinorm );
#if 1
          QDPIO::cout << "r[" <<i <<"] = " << r_rel[i] << std::endl;
#endif
   	  if ( toBool( r_rel[i]  >  invParam.RsdTarget[i]*invParam.RsdToleranceFactor  ) ) { 
	    QDPIO::cout << "Shift " << i << " has rel. residuum " << r_rel[i] <<  " exceeding target " 
			<< invParam.RsdTarget[i] << " . Aborting! " << std::endl;
	    abortFound = true;
          } 
        }
      }
      QDPIO::cout << "MULTI_CG_QUDA_CLOVER_SOLVER: " << res.n_count << " iterations. Rsd = " << res.resid << std::endl;
      QDPIO::cout << "MULTI_CG_QUDA_CLOVER_SOLVER: "<<time<< " sec" << std::endl;
      if( abortFound ) QDP_abort(1);


      END_CODE();
      
      return res;
    }

    
  private:
    // Hide default constructor
    MdagMMultiSysSolverCGQudaClover() {}
    U GFixMat;
    QudaPrecision_s cpu_prec;
    QudaPrecision_s gpu_prec;
    QudaPrecision_s gpu_half_prec;
    QudaPrecision_s gpu_ref_prec;

    Handle< LinearOperator<T> > A;
    const MultiSysSolverQUDACloverParams invParam;
    QudaGaugeParam q_gauge_param;
    mutable QudaInvertParam quda_inv_param;

    Handle< CloverTermT<T, U> > clov;
    Handle< CloverTermT<T, U> > invclov;

    SystemSolverResults_t qudaInvertMulti(const T& chi_s,
				     multi1d<T>& psi_s,
				     multi1d<Real> shifts
				     )const ;

    std::string solver_string;
    
  };


} // End namespace

#endif
#endif 


// -*- C++ -*-
// $Id: syssolver_linop_quda_clover.h,v 1.9 2009-10-16 13:37:39 bjoo Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by BiCGStab
 */

#ifndef __syssolver_linop_quda_clover_h__
#define __syssolver_linop_quda_clover_h__

#include "chroma_config.h"

#ifdef BUILD_QUDA_0_3

#include "handle.h"
#include "state.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/fermbcs/simple_fermbc.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/invert/quda_solvers/syssolver_quda_clover_params.h"
#include "actions/ferm/linop/clover_term_qdp_w.h"
#include "meas/gfix/temporal_gauge.h"
#include "io/aniso_io.h"
#include <string>

#include "util/gauge/reunit.h"

#include <quda.h>
//#include <util_quda.h>
using namespace std;

namespace Chroma
{

  //! Richardson system solver namespace
  namespace LinOpSysSolverQUDACloverEnv
  {
    //! Register the syssolver
    bool registerAll();
  }



  //! Solve a Clover Fermion System using the QUDA inverter
  /*! \ingroup invert
 *** WARNING THIS SOLVER WORKS FOR Clover FERMIONS ONLY ***
   */
 
  class LinOpSysSolverQUDAClover : public LinOpSystemSolver<LatticeFermion>
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
    LinOpSysSolverQUDAClover(Handle< LinearOperator<T> > A_,
					 Handle< FermState<T,Q,Q> > state_,
					 const SysSolverQUDACloverParams& invParam_) : 
      A(A_), invParam(invParam_), clov(new QDPCloverTermT<T, U>()), invclov(new QDPCloverTermT<T, U>())
    {
      QDPIO::cout << "LinOpSysSolverQUDAClover:" << endl;

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
      quda_inv_param.solution_type = QUDA_MATPC_SOLUTION;

      // Solve type
      switch( invParam.solverType ) { 
      case CG: 
	quda_inv_param.solve_type = QUDA_NORMEQ_PC_SOLVE;
	break;
      case BICGSTAB:
	quda_inv_param.solve_type = QUDA_DIRECT_PC_SOLVE;
	break;
      default:
	quda_inv_param.solve_type = QUDA_NORMEQ_PC_SOLVE;   
	
	break;
      }

      if( invParam.asymmetricP ) { 
	QDPIO::cout << "Using Asymmetric Linop: A_oo - D A^{-1}_ee D" << endl;
	quda_inv_param.matpc_type = QUDA_MATPC_ODD_ODD_ASYMMETRIC;
      }
      else { 
	QDPIO::cout << "Using Symmetric Linop: 1 - A^{-1}_oo D A^{-1}_ee D" << endl;
	quda_inv_param.matpc_type = QUDA_MATPC_ODD_ODD;
      }

      quda_inv_param.dagger = QUDA_DAG_NO;
      quda_inv_param.mass_normalization = QUDA_KAPPA_NORMALIZATION;

      quda_inv_param.cpu_prec = cpu_prec;
      quda_inv_param.cuda_prec = gpu_prec;
      quda_inv_param.cuda_prec_sloppy = gpu_half_prec;
      quda_inv_param.preserve_source = QUDA_PRESERVE_SOURCE_NO;
      quda_inv_param.dirac_order = QUDA_DIRAC_ORDER;

      // Autotuning
      if( invParam.tuneDslashP ) { 
	QDPIO::cout << "Enabling Dslash Autotuning" << endl;

	quda_inv_param.dirac_tune = QUDA_TUNE_YES;
      }
      else { 
	QDPIO::cout << "Disabling Dslash Autotuning" << endl;
       
	quda_inv_param.dirac_tune = QUDA_TUNE_NO;
      }

      if( invParam.cacheDslashTuningP) { 
	// Retune for every solve
	QDPIO::cout << "Will cache Dslash tuning params accross solves" << endl;
	quda_inv_param.preserve_dirac = QUDA_PRESERVE_DIRAC_YES;
      }
      else { 
	
	quda_inv_param.preserve_dirac = QUDA_PRESERVE_DIRAC_NO;
      }

      // PADDING
      q_gauge_param.ga_pad = (latdims[0]*latdims[1]*latdims[2])/2;
      quda_inv_param.sp_pad = 0;
      quda_inv_param.cl_pad = 0;

      // Clover precision and order
      quda_inv_param.clover_cpu_prec = cpu_prec;
      quda_inv_param.clover_cuda_prec = gpu_prec;
      quda_inv_param.clover_cuda_prec_sloppy = gpu_half_prec;
      quda_inv_param.clover_order = QUDA_PACKED_CLOVER_ORDER;
    
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

      //      Setup the clover term...
      QDPIO::cout << "Creating CloverTerm" << endl;
      clov->create(fstate, invParam_.CloverParams);
      // Don't recompute, just copy
      invclov->create(fstate, invParam_.CloverParams);
      
      QDPIO::cout << "Inverting CloverTerm" << endl;
      invclov->choles(0);
      invclov->choles(1);

      multi1d<QUDAPackedClovSite<REALT> > packed_clov;

      // Only compute clover if we're using asymmetric preconditioner
      if( invParam.asymmetricP ) { 
	packed_clov.resize(all.siteTable().size());

	clov->packForQUDA(packed_clov, 0);
	clov->packForQUDA(packed_clov, 1);
      }

      // Always need inverse
      multi1d<QUDAPackedClovSite<REALT> > packed_invclov(all.siteTable().size());
      invclov->packForQUDA(packed_invclov, 0);
      invclov->packForQUDA(packed_invclov, 1);





      
      if( invParam.asymmetricP ) { 
	loadCloverQuda(&(packed_clov[0]), &(packed_invclov[0]),&quda_inv_param);
      }
      else { 
	loadCloverQuda(NULL, &(packed_invclov[0]), &quda_inv_param);
      }
   
      
    }
    

    //! Destructor is automatic
    ~LinOpSysSolverQUDAClover() 
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
	QDPIO::cout << "Gauge Fixing source and initial guess" << endl;
        g_chi[ rb[1] ]  = GFixMat * chi;
	g_psi[ rb[1] ]  = GFixMat * psi;
	QDPIO::cout << "Solving" << endl;
	res = qudaInvert(*clov,
			 *invclov,
			 g_chi,
			 g_psi);      
	QDPIO::cout << "Untransforming solution." << endl;
	psi[ rb[1]]  = adj(GFixMat)*g_psi;

      }
      else { 
	res = qudaInvert(*clov,
			 *invclov,
			 chi,
			 psi);      
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

      QDPIO::cout << "QUDA_"<< solver_string <<"_CLOVER_SOLVER: " << res.n_count << " iterations. Rsd = " << res.resid << " Relative Rsd = " << rel_resid << endl;
   
      // Convergence Check/Blow Up
      if ( ! invParam.SilentFailP ) { 
	      if (  toBool( rel_resid >  invParam.RsdToleranceFactor*invParam.RsdTarget) ) { 
        	QDPIO::cerr << "ERROR: QUDA Solver residuum is outside tolerance: QUDA resid="<< rel_resid << " Desired =" << invParam.RsdTarget << " Max Tolerated = " << invParam.RsdToleranceFactor*invParam.RsdTarget << endl; 
        	QDP_abort(1);
      	      }
      }

      END_CODE();
      return res;
    }


  private:
    // Hide default constructor
    LinOpSysSolverQUDAClover() {}
    
#if 1
    Q links_orig;
#endif

    U GFixMat;
    QudaPrecision_s cpu_prec;
    QudaPrecision_s gpu_prec;
    QudaPrecision_s gpu_half_prec;

    Handle< LinearOperator<T> > A;
    const SysSolverQUDACloverParams invParam;
    QudaGaugeParam q_gauge_param;
    QudaInvertParam quda_inv_param;

    Handle< QDPCloverTermT<T, U> > clov;
    Handle< QDPCloverTermT<T, U> > invclov;

    SystemSolverResults_t qudaInvert(const QDPCloverTermT<T, U>& clover,
				     const QDPCloverTermT<T, U>& inv_clov,
				     const T& chi_s,
				     T& psi_s     
				     )const ;

    std::string solver_string;
  };


} // End namespace

#endif // BUILD_QUDA
#endif 


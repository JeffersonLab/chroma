// -*- C++ -*-
// $Id: syssolver_linop_quda_nef.h,v 1.9 2009-10-16 13:37:39 bjoo Exp $
/*! \file
*  \brief Solve a MdagM*psi=chi linear system by BiCGStab
*/

#ifndef __syssolver_linop_quda_nef_h__
#define __syssolver_linop_quda_nef_h__

#include "chroma_config.h"

#ifdef BUILD_QUDA
#include <quda.h>

#include "handle.h"
#include "state.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/fermbcs/simple_fermbc.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/invert/quda_solvers/syssolver_quda_nef_params.h"
#include "actions/ferm/linop/eoprec_nef_linop_array_w.h"
#include "meas/gfix/temporal_gauge.h"
#include <string>

#include "util/gauge/reunit.h"

//#include <util_quda.h>

namespace Chroma
{
  
  //! Richardson system solver namespace
  namespace LinOpSysSolverQUDANEFEnv
  {
    //! Register the syssolver
    bool registerAll();
  }
  
  
  
  //! Solve a Clover Fermion System using the QUDA inverter
  /*! \ingroup invert
   * WARNING THIS SOLVER WORKS FOR MOEBIUS/DWF ONLY ***
   */
 
  class LinOpSysSolverQUDANEF : public LinOpSystemSolverArray<LatticeFermion>
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
    LinOpSysSolverQUDANEF(Handle< LinearOperatorArray<T> > A_,
			  Handle< FermState<T,Q,Q> > state_,
			  const SysSolverQUDANEFParams& invParam_) : 
      A(A_), invParam(invParam_) 
    {
      START_CODE();
      
      QDPIO::cout << "LinOpSysSolverQUDANEF:" << std::endl;
      
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
      
      // Don't support anisotorpy for Moebius
      q_gauge_param.anisotropy = 1.0;
      
      // Now onto the inv param:
      // Dslash type
      quda_inv_param.dslash_type = QUDA_MOBIUS_DWF_DSLASH;
      // Invert type:
      switch( invParam.solverType ) { 
      case CG: 
	quda_inv_param.inv_type = QUDA_CG_INVERTER;
	solver_string = "CG";
	break;
      case BICGSTAB:
	QDPIO::cerr << "Solver BICGSTAB not supported for MDWF" << std::endl;
	QDP_abort(1);
	break;
      case GCR:
	QDPIO::cerr << "Solver GCR not supported for MDWF" << std::endl;
	QDP_abort(1);
	break;
      default:
	QDPIO::cerr << "Unknown Solver type" << std::endl;
	QDP_abort(1);
	break;
      }

      // Mass
      quda_inv_param.mass = toDouble(invParam.NEFParams.Mass);
      quda_inv_param.m5=toDouble(-invParam.NEFParams.OverMass);
      quda_inv_param.Ls=invParam.NEFParams.N5;
      // Mike made these static so no need to alloc them. 
      if ( invParam.NEFParams.N5 >= QUDA_MAX_DWF_LS ) { 
	QDPIO::cerr << "LS can be at most " << QUDA_MAX_DWF_LS << std::endl;
	QDP_abort(1);
      }
      
      // Copy b5 and c5 into static array
      QDPIO::cout << "Ls from matrix: " << A->size() << std::endl;
      QDPIO::cout << "Ls from params: " << invParam.NEFParams.N5 << std::endl;
      QDPIO::cout << "Ls from quda: " << quda_inv_param.Ls << std::endl;
      for(unsigned int s = 0; s < quda_inv_param.Ls; s++){
	quda_inv_param.b_5[s] = toDouble(invParam.NEFParams.b5[s]);
	quda_inv_param.c_5[s] = toDouble(invParam.NEFParams.c5[s]);
	QDPIO::cout << " b5[" <<s<<"] = " << quda_inv_param.b_5[s] << "   c5[" << s << "] = " << quda_inv_param.c_5[s] << std::endl;
      }
	  
      quda_inv_param.tol = toDouble(invParam.RsdTarget);
      quda_inv_param.maxiter = invParam.MaxIter;
      quda_inv_param.reliable_delta = toDouble(invParam.Delta);
      
      // Solution type
      quda_inv_param.solution_type = QUDA_MATPC_SOLUTION;

      // Solve type, only CG-types supported so far:
      switch( invParam.solverType ) { 
      case CG: 
	quda_inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;
	break;
	
      default:
	quda_inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;
	break;
      }

      //only symmetric DWF supported at the moment:
      QDPIO::cout << "Using Symmetric Linop: A_oo - D_oe A^{-1}_ee D_eo" << std::endl;
      quda_inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN_ASYMMETRIC;
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
      QDPIO::cout << "MDAGM Using QDP-JIT spinor order" << std::endl;
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

      
      // Setup padding
      multi1d<int> face_size(4);
      face_size[0] = latdims[1]*latdims[2]*latdims[3]/2;
      face_size[1] = latdims[0]*latdims[2]*latdims[3]/2;
      face_size[2] = latdims[0]*latdims[1]*latdims[3]/2;
      face_size[3] = latdims[0]*latdims[1]*latdims[2]/2;
      
      int max_face = face_size[0];
      for(int i=1; i < 4; i++) { 
	if ( face_size[i] > max_face ) { 
	  max_face = face_size[i]; 
	}
      }
      
      q_gauge_param.ga_pad = max_face;
      // PADDING
      quda_inv_param.sp_pad = 0;
      quda_inv_param.cl_pad = 0;
      
      QDPIO::cout << "Setting Precondition stuff to defaults for not using" << std::endl;
      quda_inv_param.inv_type_precondition= QUDA_INVALID_INVERTER;
      quda_inv_param.tol_precondition = 1.0e-1;
      quda_inv_param.maxiter_precondition = 1000;
      quda_inv_param.verbosity_precondition = QUDA_SILENT;
      quda_inv_param.gcrNkrylov = 1;
      
 
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
	std::cout << "MDAGM CUDA gauge[" << mu << "] in = " << gauge[mu] << "\n";
#endif
      }

      loadGaugeQuda((void *)gauge, &q_gauge_param); 
      
      END_CODE();
    }
    

    //! Destructor is automatic
    ~LinOpSysSolverQUDANEF() 
    {
      QDPIO::cout << "Destructing" << std::endl;
      freeGaugeQuda();
    }

    //! Return the subset on which the operator acts
    const Subset& subset() const {return A->subset();}
    
    //return the size:
    int size() const {return invParam.NEFParams.N5;}

    //! Solver the linear system
    /*!
     * \param psi      solution ( Modify )
     * \param chi      source ( Read )
     * \return syssolver results
     */
    SystemSolverResults_t operator() (multi1d<T>& psi, const multi1d<T>& chi) const
    {
      SystemSolverResults_t res;
      
      START_CODE();
      StopWatch swatch;
      swatch.start();



      // 1/( 2 kappa_b ) =  b5 * ( Nd - M5 ) + 1
      Real invTwoKappaB = invParam.NEFParams.b5[0]*( Nd - invParam.NEFParams.OverMass) + Real(1);
      Real twoKappaB = Real(1)/invTwoKappaB;

#if 0
      // Test code....
      {
	const multi1d<int>& latdims = Layout::subgridLattSize();
	int halfsize=latdims[0]*latdims[1]*latdims[2]*latdims[3]/2;
	int fermsize=halfsize*Nc*Ns*2;
	
	// In1 is input to Chroma, Out1 is result
	multi1d<T> in1( this->size() );
	multi1d<T> out1(this->size() );

	// In2 is input to QUDA, Out2 is result
	multi1d<T> in2( this->size() );
	multi1d<T> out2(this->size() );

	
	for(int s=0; s < this->size(); s++ ) { 
	  gaussian(in1[s]);  // Gaussian into in1
	  in2[s] = in1[s];   // copy to in2
	  
	}

	for(int d=0; d < 2; d++) { 
	  for(int s=0; s < this->size(); s++ ) { 
	    out1[s]=zero;   // zero both out1 and out2
	    out2[s]=zero;
	  }
	  
	  if ( d==0 ) { 
	    // Apply A to in2
	    QDPIO::cout << "DOing Mat" << std::endl;
	    (*A)(out2, in2, PLUS);
	  }
	  else {
	    QDPIO::cout << "Doing MatDag" << std::endl;
	    (*A)(out2, in2, MINUS);
	  }
	    
	// Copy in1 into QUDA
	  REAL* spinorIn = new REAL[quda_inv_param.Ls*fermsize];
	  REAL* spinorOut = new REAL[quda_inv_param.Ls*fermsize];
	  memset((spinorIn), 0, fermsize*quda_inv_param.Ls*sizeof(REAL));
	  memset((spinorOut), 0, fermsize*quda_inv_param.Ls*sizeof(REAL));
	  
	  for(unsigned int s=0; s<quda_inv_param.Ls; s++){
	    memcpy((&spinorIn[fermsize*s]),&(in1[s].elem(rb[1].start()).elem(0).elem(0).real()),fermsize*sizeof(REAL));
	  }
	  
	  // Apply QUDA
	  if( d==0 ) { 
	    quda_inv_param.dagger = QUDA_DAG_NO;
	    MatQuda((void *)spinorOut, (void *)spinorIn, (QudaInvertParam*)&quda_inv_param);
	  }
	  else {
	    quda_inv_param.dagger = QUDA_DAG_YES;
	    MatQuda((void *)spinorOut, (void *)spinorIn, (QudaInvertParam*)&quda_inv_param);
	  }

	  for(unsigned int s=0; s<quda_inv_param.Ls; s++){
	    //      memset(reinterpret_cast<char*>(&(psi_s[s].elem(all.start()).elem(0).elem(0).real())),0,fermsize*2*sizeof(REAL));
	    memcpy((&out1[s].elem(rb[1].start()).elem(0).elem(0).real()),(&spinorOut[fermsize*s]),fermsize*sizeof(REAL));
	  }
	  
	  // Now compare out1 and out2
	  for(int s=0; s < this->size();s++) { 
	    out1[s] *= invTwoKappaB;
	    
	    QDPIO::cout << "s=" << s << "  diff=" << norm2(out2[s]-out1[s]) << std::endl;
	  }


	  delete [] spinorIn;
	  delete [] spinorOut;
	}


      }
#endif

      QDPIO::cout << "Norm of chi = " << norm2(chi,rb[1]) << std::endl;
      QDPIO::cout << "TwoKappa = " << twoKappaB << "   invTwoKappa_b = " << invTwoKappaB << std::endl;
    
      if ( invParam.axialGaugeP ) { 
	multi1d<T> g_chi( this->size());
	multi1d<T> g_psi( this->size());

	// Gauge Fix source and initial guess
	QDPIO::cout << "Gauge Fixing source and initial guess" << std::endl;
	for(unsigned int s=0; s<invParam.NEFParams.N5; s++){
	  g_chi[s][ rb[1] ]  = GFixMat * chi[s];
	  g_psi[s][ rb[1] ]  = GFixMat * psi[s];
	}
	QDPIO::cout << "Solving" << std::endl;
	res = qudaInvert(g_chi,g_psi);      
	for(int s=0; s < this->size(); s++) {
	  g_psi[s][rb[1]] *= twoKappaB;
	}

	QDPIO::cout << "Untransforming solution." << std::endl;
	for(unsigned int s=0; s<invParam.NEFParams.N5; s++){
	  psi[s][ rb[1] ]  = adj(GFixMat)*g_psi[s];
	}
	
      }
      else { 
	QDPIO::cout << "Calling qudaInvert" << std::endl;
	res = qudaInvert(chi,psi);      
	for(int s=0; s < this->size(); s++) {
	  psi[s][rb[1]] *= twoKappaB;
	}
      }      
      
      

      swatch.stop();
      double time = swatch.getTimeInSeconds();
      
      // Check Solution
      {

	multi1d<T> r(A->size());
	multi1d<T> Ax(A->size());
	r=zero;
	Ax=zero;
	Double r_norm(zero);
	Double b_norm(zero);
	(*A)(Ax, psi, PLUS);
	for(int s=0; s < A->size(); s++) { 
	  r[s][rb[1]] = chi[s] - Ax[s];
	  r_norm += norm2(r[s], rb[1]);
	  b_norm += norm2(chi[s], rb[1]);
	}

	Double resid = sqrt(r_norm);
	Double rel_resid = sqrt(r_norm/b_norm);

	res.resid = resid;
	QDPIO::cout << "QUDA_"<< solver_string <<"_NEF_SOLVER: " << res.n_count << " iterations. Max. Rsd = " << res.resid << " Max. Relative Rsd = " << rel_resid << std::endl;


	// Convergence Check/Blow Up
	if ( ! invParam.SilentFailP ) { 
	  if (  toBool( rel_resid >  invParam.RsdToleranceFactor*invParam.RsdTarget) ) { 
	    QDPIO::cerr << "ERROR: QUDA Solver residuum is outside tolerance: QUDA resid="<< rel_resid << " Desired =" << invParam.RsdTarget << " Max Tolerated = " << invParam.RsdToleranceFactor*invParam.RsdTarget << std::endl; 
	    QDP_abort(1);
	  }

	}
      }

      
      END_CODE();
      return res;
    }
    
  private:
    // Hide default constructor
    LinOpSysSolverQUDANEF() {}
    
#if 1
		Q links_orig;
#endif

		U GFixMat;
		QudaPrecision_s cpu_prec;
		QudaPrecision_s gpu_prec;
		QudaPrecision_s gpu_half_prec;

		Handle< LinearOperatorArray<T> > A;
		const SysSolverQUDANEFParams invParam;
		QudaGaugeParam q_gauge_param;
		mutable QudaInvertParam quda_inv_param;

		SystemSolverResults_t qudaInvert(const multi1d<T>& chi_s, multi1d<T>& psi_s)const;

		std::string solver_string;
	};


} // End namespace

#endif // BUILD_QUDA
#endif 


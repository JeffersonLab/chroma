// $Id: syssolver_linop_quda_clover.cc,v 1.3 2009-10-05 20:19:13 bjoo Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by CG2
 */

#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"
#include "actions/ferm/invert/quda_solvers/syssolver_quda_clover_params.h"
#include "actions/ferm/invert/quda_solvers/syssolver_linop_quda_clover.h"
#include "io/aniso_io.h"


#include "handle.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/linop/lwldslash_w.h"
#include "meas/glue/mesplq.h"
// QUDA Headers
#include <quda.h>
#include <util_quda.h>

namespace Chroma
{
  namespace LinOpSysSolverQUDACloverEnv
  {

    //! Anonymous namespace
    namespace
    {
      //! Name to be used
      const std::string name("QUDA_CLOVER_INVERTER");

      //! Local registration flag
      bool registered = false;
    }



    LinOpSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,	
						  const std::string& path,
						  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state, 
						  
						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new LinOpSysSolverQUDAClover(A, state,SysSolverQUDACloverParams(xml_in, path));
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
  }

  SystemSolverResults_t 
  LinOpSysSolverQUDAClover::qudaInvert(const QDPCloverTermT<T, U>& clover,
				       const QDPCloverTermT<T, U>& invclov,
				       const T& chi_s,
				       T& psi_s) const{

    SystemSolverResults_t ret;
    QudaInvertParam inv_param;
    // Definitely no clover here...

    inv_param.dslash_type = QUDA_CLOVER_WILSON_DSLASH; // Sets Clover Matrix
    //  (1-k^2 Doe Deo)
    inv_param.matpc_type = QUDA_MATPC_ODD_ODD; 
    // inv_param.matpc_type = QUDA_MATPC_ODD_ODD_ASYMMETRIC;

    // Fiendish idea from Ron. Set the kappa=1/2 and use 
    // unmodified clover term, and ask for Kappa normalization
    // This should give us A - (1/2)D as the unpreconditioned operator
    // and probabl 1 - {1/4} A^{-1} D A^{-1} D as the preconditioned
    // op. Apart from the A_oo stuff on the antisymmetric we have
    // nothing to do...
    inv_param.kappa = 0.5;

    inv_param.tol = toDouble(invParam.RsdTarget);
    inv_param.maxiter = invParam.MaxIter;
    inv_param.reliable_delta = toDouble(invParam.Delta);


    // Solve the preconditioned matrix (rather than the prop
    inv_param.solution_type = QUDA_MATPC_SOLUTION;
    inv_param.mass_normalization = QUDA_KAPPA_NORMALIZATION;

    inv_param.cpu_prec = cpu_prec;
    inv_param.cuda_prec = cpu_prec;
    inv_param.cuda_prec_sloppy = half_prec;
    inv_param.preserve_source = QUDA_PRESERVE_SOURCE_NO;

    // Even-odd colour inside spin
    inv_param.dirac_order = QUDA_DIRAC_ORDER;

    if( invParam.verboseP ) { 
      inv_param.verbosity = QUDA_VERBOSE;
    }
    else { 
      inv_param.verbosity = QUDA_SUMMARIZE;
    }



    if ( invParam.solverType == "CG" ) { 
       inv_param.inv_type = QUDA_CG_INVERTER;   
    }
    else { 
      if( invParam.solverType == "BICGSTAB" ) { 
	inv_param.inv_type = QUDA_BICGSTAB_INVERTER;   
      }
      else { 
	QDPIO::cout << "LINOP_QUDA_CLOVER_SOLVER: Unknown solver type: " << invParam.solverType << endl;
	QDP_abort(1);
      }
    }

    // Solving  A_oo ( 1 - A^{-1}_oo D A^{-1}_ee D ) psi = chi
    // so            ( 1 - A^{-1}_oo D A^{-1}_ee D ) psi = A^{-1}_oo chi
    // So set up A^{-1}_oo chi


    T mod_chi;
    if ( inv_param.matpc_type == QUDA_MATPC_ODD_ODD_ASYMMETRIC ) {
      // asymmetric 
      mod_chi = chi_s;
    }
    else if( inv_param.matpc_type == QUDA_MATPC_ODD_ODD) { 
	invclov.apply(mod_chi, chi_s, PLUS, 1);
    }
    else { 
      QDPIO::cout << "MATPC Type not allowed." << endl;
      QDPIO::cout << " Allowed are: QUDA_MATPC_ODD_ODD_ASYMMETRIC or QUDA_MATPC_ODD_ODD" << endl;
      QDP_abort(1);
    }

    StopWatch swatch1; 
    swatch1.reset();
    swatch1.start();
 
// DEAD Test Code
#if 0
    void* spinorIn =(void *)&(mod_chi.elem(rb[1].start()).elem(0).elem(0).real());
    void* spinorOut =(void *)&(psi_s.elem(rb[0].start()).elem(0).elem(0).real());
    
    // OK Here I have a chance to test directly 
    // Even Target Checkerboard, No Dagger
    psi_s = zero;
    inv_param.dslash_type = QUDA_CLOVER_WILSON_DSLASH; // Sets Clover Matrix    
    dslashQuda(spinorOut, spinorIn, &inv_param, 0, 0);



    // Need to create a simple ferm state from the links_single...
    Handle< FermState<T, Q, Q> > pstate(new PeriodicFermState<TF,QF,QF>(links_orig));
    const AnisoParam_t& aniso = invParam.CloverParams.anisoParam;
    QDPWilsonDslashT<T,Q,Q>  qdp_dslash(pstate, aniso);
    
    T tmp,psi2;
    tmp=zero;
    psi2=zero;
    // qdp_dslash.apply(psi2, mod_chi, PLUS, 0);
    qdp_dslash.apply(tmp, mod_chi, PLUS, 0);
    invclov.apply(psi2,tmp, PLUS, 0);


    T r=zero;
    r = psi2 - psi_s;
    
    QDPIO::cout << "CB=0" << endl;
    QDPIO::cout << "Dslash Test: || r || = " << sqrt(norm2(r,rb[0])) << endl;
    //    QDPIO::cout << "Dslash Test: || r ||/|| psi || = " << sqrt(norm2(r,rb[0])/norm2(psi_s, rb[0])) << endl;

    QDPIO::cout << "CB=1: Should be zero" << endl;
    QDPIO::cout << "Dslash Test: || r || = " << sqrt(norm2(r,rb[1])) << endl;
    //QDPIO::cout << "Dslash Test: || r ||/|| psi || = " << sqrt(norm2(r,rb[1])/norm2(psi_s, rb[1])) << endl;
    
    const int* tab = rb[0].siteTable().slice();
    for(int i=0; i < rb[0].numSiteTable(); i++) { 
      int j = tab[i];
      bool printSite=false;

      for(int spin=0; spin < 4; spin++) {
	for(int col=0; col < 3; col++) { 
	  if( (fabs(r.elem(j).elem(spin).elem(col).real()) > 1.0e-5 )
	      || (fabs(r.elem(j).elem(spin).elem(col).imag()) > 1.0e-5 )) {
	    printSite=true;
	  }
	}
      }
      if( printSite ) { 
	  
	for(int spin=0; spin < 4; spin++) { 
	  for(int col=0; col < 3; col++) { 
	    QDPIO::cout << "Site= " << j << " Spin= "<< spin << " Col= " << col << " spinor = ( " 
			<< psi2.elem(j).elem(spin).elem(col).real()  << " , " 
			<< psi2.elem(j).elem(spin).elem(col).imag()  << " )" << endl;
	  }
	}
	QDPIO::cout << endl;
      }
    }
    QDP_abort(1);
#else

   void* spinorIn =(void *)&(mod_chi.elem(rb[1].start()).elem(0).elem(0).real());
    void* spinorOut =(void *)&(psi_s.elem(rb[1].start()).elem(0).elem(0).real());

    invertQuda(spinorOut, spinorIn, &inv_param);
#endif
    // Take care of mass normalization
    //psi_s *= (invMassParam);
    swatch1.stop();


    QDPIO::cout << "Cuda Space Required" << endl;
    QDPIO::cout << "\t Spinor:" << inv_param.spinorGiB << " GiB" << endl;
    QDPIO::cout << "\t Gauge :" << q_gauge_param.gaugeGiB << " GiB" << endl;
    QDPIO::cout << "\t InvClover :" << inv_param.cloverGiB << " GiB" << endl;
    QDPIO::cout << "QUDA_"<<invParam.solverType<<"_CLOVER_SOLVER: time="<< inv_param.secs <<" s" ;
    QDPIO::cout << "\tPerformance="<<  inv_param.gflops/inv_param.secs<<" GFLOPS" ; 
    QDPIO::cout << "\tTotal Time (incl. load gauge)=" << swatch1.getTimeInSeconds() <<" s"<<endl;

    ret.n_count =inv_param.iter;

    return ret;

  }
  

}



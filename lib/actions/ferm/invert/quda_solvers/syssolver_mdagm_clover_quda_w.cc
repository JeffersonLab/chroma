// $Id: syssolver_linop_quda_clover.cc,v 1.6 2009-10-09 13:59:46 bjoo Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by CG2
 */

#include "actions/ferm/invert/syssolver_mdagm_factory.h"
#include "actions/ferm/invert/syssolver_mdagm_aggregate.h"
#include "actions/ferm/invert/quda_solvers/syssolver_quda_clover_params.h"
#include "actions/ferm/invert/quda_solvers/syssolver_mdagm_clover_quda_w.h"
#include "io/aniso_io.h"


#include "handle.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/linop/lwldslash_w.h"
#include "meas/glue/mesplq.h"
// QUDA Headers
#include <quda.h>
// #include <util_quda.h>


//#undef BUILD_QUDA_DEVIFACE_GAUGE
//#undef BUILD_QUDA_DEVIFACE_SPINOR
//#undef BUILD_QUDA_DEVIFACE_CLOVER


namespace Chroma
{
  namespace MdagMSysSolverQUDACloverEnv
  {

    //! Anonymous namespace
    namespace
    {
      //! Name to be used
      const std::string name("QUDA_CLOVER_INVERTER");

      //! Local registration flag
      bool registered = false;
    }



    MdagMSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,	
						  const std::string& path,
						  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state, 
						  
						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new MdagMSysSolverQUDAClover(A, state,SysSolverQUDACloverParams(xml_in, path));
    }

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheMdagMFermSystemSolverFactory::Instance().registerObject(name, createFerm);
	registered = true;
      }
      return success;
    }
  }

  SystemSolverResults_t 
  MdagMSysSolverQUDAClover::qudaInvert(const typename CloverTermT<T, U>::Type_t& clover,
				       const typename CloverTermT<T, U>::Type_t& invclov,
				       const T& chi_s,
				       T& psi_s) const{

    SystemSolverResults_t ret;

    void *spinorIn;

    T mod_chi;
    if ( quda_inv_param.matpc_type == QUDA_MATPC_ODD_ODD_ASYMMETRIC ) {

      // Because of the vaguaries of our HMC Formulation we need to 
      // solve the Asymmetric system. Symmetric won't work unless we change
      // preconditioning strategy
      // asymmetric 
      //
      // Solve A_oo - D A^{-1}_ee D -- chroma conventions.
      // No need to transform source
#ifndef BUILD_QUDA_DEVIFACE_SPINOR
      spinorIn =(void *)&(chi_s.elem(rb[1].start()).elem(0).elem(0).real());
#else
      spinorIn = QDPCache::Instance().getDevicePtr( chi_s.getId() );
#endif
    }
    else { 
      QDPIO::cout << "MATPC Type not allowed." << endl;
      QDP_abort(1);
    }

#ifndef BUILD_QUDA_DEVIFACE_SPINOR
    void* spinorOut =(void *)&(psi_s.elem(rb[1].start()).elem(0).elem(0).real());
#else
    void* spinorOut = QDPCache::Instance().getDevicePtr( psi_s.getId() );
#endif

    // Do the solve here 
    StopWatch swatch1; 
    swatch1.reset();
    swatch1.start();
    invertQuda(spinorOut, spinorIn, (QudaInvertParam*)&quda_inv_param);
    swatch1.stop();


    QDPIO::cout << "Cuda Space Required" << endl;
    QDPIO::cout << "\t Spinor:" << quda_inv_param.spinorGiB << " GiB" << endl;
    QDPIO::cout << "\t Gauge :" << q_gauge_param.gaugeGiB << " GiB" << endl;
    QDPIO::cout << "\t InvClover :" << quda_inv_param.cloverGiB << " GiB" << endl;
    QDPIO::cout << "QUDA_"<<solver_string<<"_CLOVER_SOLVER: time="<< quda_inv_param.secs <<" s" ;
    QDPIO::cout << "\tPerformance="<<  quda_inv_param.gflops/quda_inv_param.secs<<" GFLOPS" ; 
    QDPIO::cout << "\tTotal Time (incl. load gauge)=" << swatch1.getTimeInSeconds() <<" s"<<endl;

    ret.n_count =quda_inv_param.iter;

    return ret;

  }
  

}


 
// DEAD Test Code
#if 0
    void* spinorIn =(void *)&(mod_chi.elem(rb[1].start()).elem(0).elem(0).real());
    void* spinorOut =(void *)&(psi_s.elem(rb[0].start()).elem(0).elem(0).real());
    
    // OK Here I have a chance to test directly 
    // Even Target Checkerboard, No Dagger
    psi_s = zero;
    quda_inv_param.dslash_type = QUDA_CLOVER_WILSON_DSLASH; // Sets Clover Matrix    
    dslashQuda(spinorOut, spinorIn, &quda_inv_param, 0, 0);



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
#endif

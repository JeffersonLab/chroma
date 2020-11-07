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
  MdagMSysSolverQUDAClover::qudaInvert(const CloverTermT<T, U>& clover,
				       const CloverTermT<T, U>& invclov,
				       const T& chi_s,
				       T& psi_s) const{

    SystemSolverResults_t ret;

    void *spinorIn;
    void *spinorOut;

#ifdef BUILD_QUDA_DEVIFACE_SPINOR
    std::vector<QDPCache::ArgKey> ids;
#endif
  
      // Solve A_oo - D A^{-1}_ee D -- chroma conventions.
      // No need to transform source
#ifndef BUILD_QUDA_DEVIFACE_SPINOR
      spinorIn =(void *)&(chi_s.elem(rb[1].start()).elem(0).elem(0).real());
#else
      //spinorIn = GetMemoryPtr( chi_s.getId() );
      ids.push_back(chi_s.getId());
#endif
  
#ifndef BUILD_QUDA_DEVIFACE_SPINOR
    spinorOut =(void *)&(psi_s.elem(rb[1].start()).elem(0).elem(0).real());
#else
    ids.push_back(psi_s.getId());
    auto dev_ptr = GetMemoryPtr(ids);
    spinorIn  = dev_ptr[0];
    spinorOut = dev_ptr[1];
    //void* spinorOut = GetMemoryPtr( psi_s.getId() );
#endif

    // Do the solve here 
    StopWatch swatch1; 
    swatch1.reset();
    swatch1.start();
    invertQuda(spinorOut, spinorIn, (QudaInvertParam*)&quda_inv_param);
    swatch1.stop();



    QDPIO::cout << "QUDA_"<<solver_string<<"_CLOVER_SOLVER: time="<< quda_inv_param.secs <<" s" ;
    QDPIO::cout << "\tPerformance="<<  quda_inv_param.gflops/quda_inv_param.secs<<" GFLOPS" ; 
    QDPIO::cout << "\tTotal Time (incl. load gauge)=" << swatch1.getTimeInSeconds() <<" s"<<std::endl;

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
    
    QDPIO::cout << "CB=0" << std::endl;
    QDPIO::cout << "Dslash Test: || r || = " << sqrt(norm2(r,rb[0])) << std::endl;
    //    QDPIO::cout << "Dslash Test: || r ||/|| psi || = " << sqrt(norm2(r,rb[0])/norm2(psi_s, rb[0])) << std::endl;

    QDPIO::cout << "CB=1: Should be zero" << std::endl;
    QDPIO::cout << "Dslash Test: || r || = " << sqrt(norm2(r,rb[1])) << std::endl;
    //QDPIO::cout << "Dslash Test: || r ||/|| psi || = " << sqrt(norm2(r,rb[1])/norm2(psi_s, rb[1])) << std::endl;
    
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
			<< psi2.elem(j).elem(spin).elem(col).imag()  << " )" << std::endl;
	  }
	}
	QDPIO::cout << std::endl;
      }
    }
    QDP_abort(1);
#endif

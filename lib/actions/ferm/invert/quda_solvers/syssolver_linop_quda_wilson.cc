// $Id: syssolver_linop_quda_wilson.cc,v 1.4 2009-10-05 20:19:13 bjoo Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by CG2
 */

#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"
#include "actions/ferm/invert/quda_solvers/syssolver_quda_wilson_params.h"
#include "actions/ferm/invert/quda_solvers/syssolver_linop_quda_wilson.h"
#include "io/aniso_io.h"


#include "handle.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/linop/lwldslash_w.h"
#include "meas/glue/mesplq.h"
// QUDA Headers
#include <quda.h>
#include <util_quda.h>
#include <invert_quda.h>

namespace Chroma
{
  namespace LinOpSysSolverQUDAWilsonEnv
  {

    //! Anonymous namespace
    namespace
    {
      //! Name to be used
      const std::string name("QUDA_WILSON_INVERTER");

      //! Local registration flag
      bool registered = false;
    }



    LinOpSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,	
						  const std::string& path,
						  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state, 
						  
						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new LinOpSysSolverQUDAWilson(A, state,SysSolverQUDAWilsonParams(xml_in, path));
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
  LinOpSysSolverQUDAWilson::qudaInvert(const T& chi_s,
				       T& psi_s) const{

    SystemSolverResults_t ret;
    
    QudaInvertParam inv_param;
    

    // Definitely no clover here...
    inv_param.dslash_type = QUDA_WILSON_DSLASH; // Sets Wilson Matrix

    float massParam = 1.0 + 3.0/q_gauge_param.anisotropy+ toDouble(invParam.WilsonParams.Mass); 
    float invMassParam = 1.0/massParam;

    inv_param.kappa = 1.0/(2*massParam);
    inv_param.tol = toDouble(invParam.RsdTarget);
    inv_param.maxiter = invParam.MaxIter;
    inv_param.reliable_delta = toDouble(invParam.Delta);

    //  (1-k^2 Doe Deo)
    inv_param.matpc_type = QUDA_MATPC_ODD_ODD; 

    // Solve the preconditioned matrix (rather than the prop
    inv_param.solution_type = QUDA_MATPC_SOLUTION;
    inv_param.mass_normalization = QUDA_KAPPA_NORMALIZATION;



    inv_param.cpu_prec = cpu_prec;
    inv_param.cuda_prec = cpu_prec;
    inv_param.cuda_prec_sloppy = half_prec;
    inv_param.preserve_source = QUDA_PRESERVE_SOURCE_YES;

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
	QDPIO::cout << "LINOP_QUDA_WILSON_SOLVER: Unknown solver type: " << invParam.solverType << endl;
	QDP_abort(1);
      }
    }


    StopWatch swatch1; 
    swatch1.reset();
    swatch1.start();

    void* spinorIn =(void *)&(chi_s.elem(rb[1].start()).elem(0).elem(0).real());
    void* spinorOut =(void *)&(psi_s.elem(rb[1].start()).elem(0).elem(0).real());

    invertQuda(spinorOut, spinorIn, &inv_param);

    // Take care of mass normalization
    psi_s *= (invMassParam);
    swatch1.stop();


    QDPIO::cout << "Cuda Space Required" << endl;
    QDPIO::cout << "\t Spinor:" << inv_param.spinorGiB << " GiB" << endl;
    QDPIO::cout << "\t Gauge :" << q_gauge_param.gaugeGiB << " GiB" << endl;
    QDPIO::cout << "QUDA_" << invParam.solverType << "_WILSON_SOLVER: time="<< inv_param.secs <<" s" ;
    QDPIO::cout << "\tPerformance="<<  inv_param.gflops/inv_param.secs<<" GFLOPS" ; 
    QDPIO::cout << "\tTotal Time (incl. load gauge)=" << swatch1.getTimeInSeconds() <<" s"<<endl;

    ret.n_count =inv_param.iter;

    return ret;

  }
  

}


// DEAD Test Code
#if 0
    // OK Here I have a chance to test directly 
    // Even Target Checkerboard, No Dagger
    dslashQuda(spinorOut, spinorIn, &inv_param, 0, 0);

    // Need to create a simple ferm state from the links...
    Handle< FermState<T, Q, Q> > pstate(new PeriodicFermState<T,Q,Q>(links_single));
    QDPWilsonDslashT<T,Q,Q>  qdp_dslash(pstate, aniso);
    


    qdp_dslash.apply(psi2, chi_s, PLUS, 0);

    T r=zero;
    r = psi2 - psi_s;
    
    QDPIO::cout << "CB=0" << endl;
    QDPIO::cout << "Dslash Test: || r || = " << sqrt(norm2(r,rb[0])) << endl;
    QDPIO::cout << "Dslash Test: || r ||/|| psi || = " << sqrt(norm2(r,rb[0])/norm2(psi_s, rb[0])) << endl;

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
			<< r.elem(j).elem(spin).elem(col).real()  << " , " 
			<< r.elem(j).elem(spin).elem(col).imag()  << " )" << endl;
	  }
	}
	QDPIO::cout << endl;
      }
    }
#endif

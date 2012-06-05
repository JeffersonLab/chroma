// $Id: multi_syssolver_mdagm_cg.cc,v 3.2 2006-09-20 20:28:00 edwards Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by CG2
 */

#include "actions/ferm/invert/multi_syssolver_mdagm_factory.h"
#include "actions/ferm/invert/multi_syssolver_mdagm_aggregate.h"
#include "actions/ferm/invert/quda_solvers/multi_syssolver_mdagm_cg_clover_quda_w.h"
#include "actions/ferm/invert/quda_solvers/syssolver_quda_clover_params.h"
#include "quda.h"

#include <cstdlib>
using namespace std;

namespace Chroma
{

  //! CG2 system solver namespace
  namespace MdagMMultiSysSolverCGQudaCloverEnv
  {
    //! Callback function
    MdagMMultiSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						       const std::string& path,
						       Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state, 
						       Handle< LinearOperator<LatticeFermion> > A)
    {
      return new MdagMMultiSysSolverCGQudaClover(A, state,SysSolverQUDACloverParams(xml_in, path));
    }

    //! Name to be used
    const std::string name("MULTI_CG_QUDA_CLOVER_INVERTER");

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheMdagMFermMultiSystemSolverFactory::Instance().registerObject(name, createFerm);
	registered = true;
      }
      return success;
    }
  }

  SystemSolverResults_t 
  MdagMMultiSysSolverCGQudaClover::qudaInvertMulti(const T& chi_s,
				       multi1d<T>& psi_s,
				       const multi1d<Real> shifts) const{

    SystemSolverResults_t ret;

    void *spinorIn;

    if ( quda_inv_param.matpc_type == QUDA_MATPC_ODD_ODD_ASYMMETRIC ) {

      // Because of the vaguaries of our HMC Formulation we need to 
      // solve the Asymmetric system. Symmetric won't work unless we change
      // preconditioning strategy
      // asymmetric 
      //
      // Solve A_oo - D A^{-1}_ee D -- chroma conventions.
      // No need to transform source
      spinorIn =(void *)&(chi_s.elem(rb[1].start()).elem(0).elem(0).real());
    }
    else { 
      QDPIO::cout << "MATPC Type not allowed." << endl;
      QDP_abort(1);
    }

    void** spinorOut;
    spinorOut = (void **)malloc(shifts.size()*sizeof(void *));
    if (spinorOut == NULL) { 
      QDPIO::cerr << "Couldn't allocate spinorOut" << endl;
      QDP_abort(1);
    }

    double* shift_d = new double [ shifts.size() ];
    double resid_sq;
    
    psi_s.resize( shifts.size());

    for(int s=0; s < shifts.size(); s++) {
      psi_s[s][ rb[1] ] = zero;
      spinorOut[s] = (void *)&(psi_s[s].elem(rb[1].start()).elem(0).elem(0).real());
      shift_d[s] = toDouble(shifts[s]);
    }
   // Do the solve here 
    StopWatch swatch1; 
    swatch1.reset();
    swatch1.start();
    QDPIO::cout << "CALLING QUDA SOLVER" << endl << flush ; 
    invertMultiShiftQuda(spinorOut, spinorIn, (QudaInvertParam*)&quda_inv_param, shift_d, shifts.size(), &resid_sq);
    swatch1.stop();

    // Tidy Up
    delete [] spinorOut;
    delete [] shift_d;

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


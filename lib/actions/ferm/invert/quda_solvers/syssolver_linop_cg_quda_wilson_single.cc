// $Id: syssolver_linop_cg_quda_wilson_single.cc,v 1.1 2009-09-25 12:41:23 bjoo Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by CG2
 */

#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"
#include "actions/ferm/invert/quda_solvers/syssolver_cg_quda_wilson_params.h"
#include "actions/ferm/invert/quda_solvers/syssolver_linop_cg_quda_wilson_single.h"
#include "io/aniso_io.h"
// QUDA Headers
#include <quda.h>
#include <util_quda.h>

namespace Chroma
{
  namespace LinOpSysSolverCGQUDAWilsonEnv
  {

    //! Anonymous namespace
    namespace
    {
      //! Name to be used
      const std::string name("QUDA_CG_WILSON_INVERTER");

      //! Local registration flag
      bool registered = false;
    }



    LinOpSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,	
						  const std::string& path,
						  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state, 
						  
						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new LinOpSysSolverCGQUDAWilson(A, state,SysSolverCGQUDAWilsonParams(xml_in, path));
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
  LinOpSysSolverCGQUDAWilson::qudaInvert(const QF& links, 
					 const TF& chi_s,
					 TF& psi_s,      
					 const SysSolverCGQUDAWilsonParams& invParam) const{

    int device = 0;

    SystemSolverResults_t ret;
    QudaGaugeParam q_gauge_param;
    QudaInvertParam inv_param;
    const multi1d<int>& latdims = Layout::lattSize();

    QF links_trans;

    // Kappa norm chi

    q_gauge_param.X[0] = latdims[0];
    q_gauge_param.X[1] = latdims[1];
    q_gauge_param.X[2] = latdims[2];
    q_gauge_param.X[3] = latdims[3];

    QDPIO::cout << "Gauge Param X[0] = " <<q_gauge_param.X[0] << endl;
    QDPIO::cout << "Gauge Param X[1] = " <<q_gauge_param.X[1] << endl;
    QDPIO::cout << "Gauge Param X[2] = " <<q_gauge_param.X[2] << endl;
    QDPIO::cout << "Gauge Param X[3] = " <<q_gauge_param.X[3] << endl;

    const AnisoParam_t& aniso = invParam.WilsonParams.anisoParam;

    if( aniso.anisoP ) {                     // Anisotropic case
      Real gamma_f = aniso.nu / aniso.xi_0; 
      q_gauge_param.anisotropy = toDouble(gamma_f);
    }
    else {
      q_gauge_param.anisotropy = 1.0;
    }

    if( invParam.AntiPeriodicT ) { 
      q_gauge_param.t_boundary = QUDA_ANTI_PERIODIC_T;
    }
    else { 
      q_gauge_param.t_boundary = QUDA_PERIODIC_T;
    }

    q_gauge_param.gauge_order = QUDA_QDP_GAUGE_ORDER; // Could be QDP...
    q_gauge_param.cpu_prec = QUDA_SINGLE_PRECISION;  // Single Prec G-field
    q_gauge_param.cuda_prec = QUDA_SINGLE_PRECISION; 
    q_gauge_param.reconstruct = QUDA_RECONSTRUCT_12;
    q_gauge_param.cuda_prec_sloppy = QUDA_SINGLE_PRECISION; 
    q_gauge_param.reconstruct_sloppy = QUDA_RECONSTRUCT_12;
    q_gauge_param.gauge_fix = QUDA_GAUGE_FIXED_NO;  // No Gfix yet

    q_gauge_param.blockDim = 64;         // I copy these 
    q_gauge_param.blockDim_sloppy = 64;
    
    // OK! This is ugly: gauge_param is an 'extern' in dslash_quda.h
    gauge_param = &q_gauge_param;


    // Definitely no clover here...
    inv_param.dslash_type = QUDA_WILSON_DSLASH; // Sets Wilson Matrix
    inv_param.inv_type = QUDA_CG_INVERTER;      // CG Inverter

    float massParam = 1.0 + 3.0/q_gauge_param.anisotropy+ toDouble(invParam.WilsonParams.Mass); 
    float invMassParam = 1.0/massParam;

    
    inv_param.kappa = 1.0/(2*massParam);

    QDPIO::cout << "New Kappa is " << inv_param.kappa << endl;
    inv_param.tol = toDouble(invParam.RsdTarget);
    inv_param.maxiter = invParam.MaxIter;
    inv_param.reliable_delta = toDouble(invParam.Delta);

    //  (1-k^2 Doe Deo)
    inv_param.matpc_type = QUDA_MATPC_ODD_ODD; 

    // Solve the preconditioned matrix (rather than the prop
    inv_param.solution_type = QUDA_MATPC_SOLUTION;
    inv_param.mass_normalization = QUDA_KAPPA_NORMALIZATION;

    inv_param.cpu_prec = QUDA_SINGLE_PRECISION;
    inv_param.cuda_prec = QUDA_SINGLE_PRECISION;
    inv_param.cuda_prec_sloppy = QUDA_SINGLE_PRECISION;
    inv_param.preserve_source = QUDA_PRESERVE_SOURCE_NO;

    // Even-odd colour inside spin
    inv_param.dirac_order = QUDA_QDP_DIRAC_ORDER;
    inv_param.verbosity = QUDA_VERBOSE;

    initQuda(device);
    void* gauge[4];
    for(int mu=0; mu < Nd; mu++) { 
      gauge[mu] = (void *)&(links[mu].elem(all.start()).elem().elem(0,0).real());
    }

    loadGaugeQuda((void *)gauge, &q_gauge_param);
    void* spinorIn =(void *)&(chi_s.elem(rb[1].start()).elem(0).elem(0).real());
    void* spinorOut =(void *)&(psi_s.elem(rb[1].start()).elem(0).elem(0).real());
    



    invertQuda(spinorOut, spinorIn, &inv_param);

    //psi_s *= (invMassParam);

    return ret;

  }
  

}

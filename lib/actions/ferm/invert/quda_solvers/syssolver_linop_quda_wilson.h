// -*- C++ -*-
// $Id: syssolver_linop_quda_wilson.h,v 1.3 2009-10-05 20:19:13 bjoo Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by BiCGStab
 */

#ifndef __syssolver_linop_cg_cuda_wilson_single_h__
#define __syssolver_linop_cg_coda_wilson_single_h__

#include "chroma_config.h"

#ifdef BUILD_QUDA

#include "handle.h"
#include "state.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/fermbcs/simple_fermbc.h"
#include "actions/ferm/invert/quda_solvers/syssolver_quda_wilson_params.h"
#include "io/aniso_io.h"
#include <string>

#include <quda.h>
#include <util_quda.h>
using namespace std;

namespace Chroma
{

  //! Richardson system solver namespace
  namespace LinOpSysSolverQUDAWilsonEnv
  {
    //! Register the syssolver
    bool registerAll();
  }



  //! Solve a Wilson Fermion System using the QUDA inverter
  /*! \ingroup invert
 *** WARNING THIS SOLVER WORKS FOR Wilson FERMIONS ONLY ***
   */
 
  class LinOpSysSolverQUDAWilson : public LinOpSystemSolver<LatticeFermion>
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


    //! Constructor
    /*!
     * \param M_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    LinOpSysSolverQUDAWilson(Handle< LinearOperator<T> > A_,
					 Handle< FermState<T,Q,Q> > state_,
					 const SysSolverQUDAWilsonParams& invParam_) : 
      A(A_), invParam(invParam_) 
    {
      QDPIO::cout << "LinOpSysSolverQUDAWilson:" << endl;
      const AnisoParam_t& aniso = invParam.WilsonParams.anisoParam;
      
     

      // These are the links
      // They may be smeared and the BC's may be applied
      Q links_single(Nd);
      
      // Now downcast to single prec fields.
      for(int mu=0; mu < Nd; mu++) {
	links_single[mu] = (state_->getLinks())[mu];
      }
      if( aniso.anisoP ) {                     // Anisotropic case
	multi1d<Real> cf=makeFermCoeffs(aniso);
	for(int mu=0; mu < Nd; mu++) { 
	  links_single[mu] *= cf[mu];
	}
      }

      const multi1d<int>& latdims = Layout::lattSize();
      
      q_gauge_param.X[0] = latdims[0];
      q_gauge_param.X[1] = latdims[1];
      q_gauge_param.X[2] = latdims[2];
      q_gauge_param.X[3] = latdims[3];
      
      if( aniso.anisoP ) {                     // Anisotropic case
	Real gamma_f = aniso.xi_0 / aniso.nu; 
	q_gauge_param.anisotropy = toDouble(gamma_f);
      }
      else {
	q_gauge_param.anisotropy = 1.0;
      }
      
      // Convention: BC has to be applied already
      // This flag just tells QUDA that this is so,
      // so that QUDA can take care in the reconstruct
      if( invParam.AntiPeriodicT ) { 
	q_gauge_param.t_boundary = QUDA_ANTI_PERIODIC_T;
      }
      else { 
	q_gauge_param.t_boundary = QUDA_PERIODIC_T;
      }


      q_gauge_param.gauge_order = QUDA_QDP_GAUGE_ORDER; // gauge[mu], p, col col

      int s = sizeof( WordType<T>::Type_t );
      if (s == 4 ) { 
	cpu_prec = QUDA_SINGLE_PRECISION;
	half_prec = QUDA_SINGLE_PRECISION; 
      }
      else { 
	cpu_prec = QUDA_DOUBLE_PRECISION;
	half_prec = QUDA_SINGLE_PRECISION;
      }

      q_gauge_param.cpu_prec = cpu_prec;  // Single Prec G-field
      q_gauge_param.cuda_prec = cpu_prec; 
      q_gauge_param.reconstruct = QUDA_RECONSTRUCT_12;
      q_gauge_param.cuda_prec_sloppy = half_prec;
      q_gauge_param.reconstruct_sloppy = QUDA_RECONSTRUCT_12; // No Sloppy
      
      // Do I want to Gauge Fix? -- Not yet
      q_gauge_param.gauge_fix = QUDA_GAUGE_FIXED_NO;  // No Gfix yet
      
      q_gauge_param.blockDim = 64;         // I copy these from invert test
      q_gauge_param.blockDim_sloppy = 64;
      
      // OK! This is ugly: gauge_param is an 'extern' in dslash_quda.h
      gauge_param = &q_gauge_param;
      
      // Set up the links
      void* gauge[4];
      for(int mu=0; mu < Nd; mu++) { 
	gauge[mu] = (void *)&(links_single[mu].elem(all.start()).elem().elem(0,0).real());
      }
      loadGaugeQuda((void *)gauge, &q_gauge_param);
      

    }

    //! Destructor is automatic
    ~LinOpSysSolverQUDAWilson() 
    {
      
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

      //   TD psi_d = psi;
      // TD chi_d = chi;


      // Call the QUDA Thingie here
      res = qudaInvert(chi,
		       psi);      


      //psi = psi_d;

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
      QDPIO::cout << "QUDA_" << invParam.solverType << "_WILSON_SOLVER: " << res.n_count << " iterations. Rsd = " << res.resid << " Relative Rsd = " << res.resid/sqrt(norm2(chi,A->subset())) << endl;
   
      
      END_CODE();
      return res;
    }


  private:
    // Hide default constructor
    LinOpSysSolverQUDAWilson() {}
    QudaPrecision_s cpu_prec;
    QudaPrecision_s half_prec;

    Handle< LinearOperator<T> > A;
    const SysSolverQUDAWilsonParams invParam;
    QudaGaugeParam q_gauge_param;
    SystemSolverResults_t qudaInvert(const TD& chi_d,
				     TD& psi_d     
				     )const ;


  };


} // End namespace

#endif // BUILD_QUDA
#endif 


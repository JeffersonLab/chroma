// -*- C++ -*-
// $Id: syssolver_linop_quda_clover.h,v 1.2 2009-10-01 20:21:53 bjoo Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by BiCGStab
 */

#ifndef __syssolver_linop_quda_clover_h__
#define __syssolver_linop_quda_clover_h__

#include "chroma_config.h"

#ifdef BUILD_QUDA

#include "handle.h"
#include "state.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/fermbcs/simple_fermbc.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/invert/quda_solvers/syssolver_quda_clover_params.h"
#include "actions/ferm/linop/clover_term_qdp_w.h"
#include "io/aniso_io.h"
#include <string>

#include <quda.h>
#include <util_quda.h>
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


    //! Constructor
    /*!
     * \param M_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    LinOpSysSolverQUDAClover(Handle< LinearOperator<T> > A_,
					 Handle< FermState<T,Q,Q> > state_,
					 const SysSolverQUDACloverParams& invParam_) : 
      A(A_), invParam(invParam_), clov(new QDPCloverTermT<TF, UF>()), invclov(new QDPCloverTermT<TF, UF>())
    {
      QDPIO::cout << "LinOpSysSolverQUDAClover:" << endl;
      const AnisoParam_t& aniso = invParam.CloverParams.anisoParam;
      
      // These are the links
      // They may be smeared and the BC's may be applied
      links_single.resize(Nd);
      
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
      q_gauge_param.cpu_prec = QUDA_SINGLE_PRECISION;  // Single Prec G-field
      q_gauge_param.cuda_prec = QUDA_SINGLE_PRECISION; 
      q_gauge_param.reconstruct = QUDA_RECONSTRUCT_12;
      q_gauge_param.cuda_prec_sloppy = QUDA_SINGLE_PRECISION; // No Sloppy
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
      

      
      // These are the links
      // They may be smeared and the BC's may be applied
      links_orig.resize(Nd);
      
      // Now downcast to single prec fields.
      for(int mu=0; mu < Nd; mu++) {
	links_orig[mu] = (state_->getLinks())[mu];
      }

      Handle<FermState<TF,QF,QF> > fstate( new PeriodicFermState<TF,QF,QF>(links_orig));

      QDPIO::cout << "Creating CloverTerm" << endl;
      clov->create(fstate, invParam_.CloverParams);
      // Don't recompute, just copy
      invclov->create(fstate, invParam_.CloverParams);
      
      QDPIO::cout << "Inverting CloverTerm" << endl;
      invclov->choles(0);
      invclov->choles(1);
      multi1d<QUDAPackedClovSite<REAL> > packed_invclov(all.siteTable().size());
      invclov->packForQUDA(packed_invclov, 0);
      invclov->packForQUDA(packed_invclov, 1);


      inv_param.clover_cpu_prec = QUDA_SINGLE_PRECISION;
      inv_param.clover_cuda_prec = QUDA_SINGLE_PRECISION;
      inv_param.clover_cuda_prec_sloppy = QUDA_SINGLE_PRECISION;
      inv_param.clover_order = QUDA_PACKED_CLOVER_ORDER;

      loadCloverQuda(NULL, &(packed_invclov[0]), &inv_param);

    }


    //! Destructor is automatic
    ~LinOpSysSolverQUDAClover() {
      discardCloverQuda(&inv_param);


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

      TF psi_s = psi;
      TF chi_s = chi;


      // Call the QUDA Thingie here
      res = qudaInvert(*clov,
		       *invclov,
		       chi_s,
		       psi_s);      


      psi = psi_s;

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
      QDPIO::cout << "QUDA_CG_CLOVER_SOLVER: " << res.n_count << " iterations. Rsd = " << res.resid << " Relative Rsd = " << res.resid/sqrt(norm2(chi,A->subset())) << endl;
   
      
      END_CODE();
      return res;
    }


  private:
    // Hide default constructor
    LinOpSysSolverQUDAClover() {}
    
    QF links_single;
    QF links_orig;

    Handle< LinearOperator<T> > A;
    const SysSolverQUDACloverParams invParam;
    QudaGaugeParam q_gauge_param;
    QudaInvertParam inv_param;

    Handle< QDPCloverTermT<TF, UF> > clov;
    Handle< QDPCloverTermT<TF, UF> > invclov;

    SystemSolverResults_t qudaInvert(const QDPCloverTermT<TF, UF>& clover,
				     const QDPCloverTermT<TF, UF>& inv_clov,
				     const TF& chi_s,
				     TF& psi_s     
				     )const ;


  };


} // End namespace

#endif // BUILD_QUDA
#endif 


// -*- C++ -*-
// $Id: syssolver_linop_quda_clover.h,v 1.1 2009-09-29 23:10:30 bjoo Exp $
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
      links_single; links_single.resize(Nd);
      
      // Now downcast to single prec fields.
      for(int mu=0; mu < Nd; mu++) {
	links_single[mu] = (state_->getLinks())[mu];
      }

      Handle<FermState<TF,QF,QF> > fstate( new PeriodicFermState<TF,QF,QF>(links_single));


      QDPIO::cout << "Creating CloverTerm" << endl;
      clov->create(fstate, invParam_.CloverParams);
      invclov->create(fstate, invParam_.CloverParams);
      
      QDPIO::cout << "Inverting CloverTerm" << endl;
      invclov->choles(0);
      invclov->choles(1);

      if( aniso.anisoP ) {                     // Anisotropic case
	multi1d<Real> cf=makeFermCoeffs(aniso);
	for(int mu=0; mu < Nd; mu++) { 
	  links_single[mu] *= cf[mu];
	}
      }

    }


    //! Destructor is automatic
    ~LinOpSysSolverQUDAClover() {}

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
      res = qudaInvert(links_single, 
		       *clov,
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
    Q links_floating;
    Handle< LinearOperator<T> > A;
    const SysSolverQUDACloverParams invParam;
    Handle< QDPCloverTermT<TF, UF> > clov;
    Handle< QDPCloverTermT<TF, UF> > invclov;

    SystemSolverResults_t qudaInvert(const QF& links, 
				     const QDPCloverTermT<TF, UF>& clover,
				     const QDPCloverTermT<TF, UF>& inv_clov,
				     const TF& chi_s,
				     TF& psi_s     
				     )const ;


  };


} // End namespace

#endif // BUILD_QUDA
#endif 


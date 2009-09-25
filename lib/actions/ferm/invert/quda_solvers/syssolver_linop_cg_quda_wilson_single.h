// -*- C++ -*-
// $Id: syssolver_linop_cg_quda_wilson_single.h,v 1.3 2009-09-25 19:00:44 bjoo Exp $
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
#include "actions/ferm/invert/quda_solvers/syssolver_cg_quda_wilson_params.h"
#include "io/aniso_io.h"
#include <string>
using namespace std;

namespace Chroma
{

  //! Richardson system solver namespace
  namespace LinOpSysSolverCGQUDAWilsonEnv
  {
    //! Register the syssolver
    bool registerAll();
  }



  //! Solve a Wilson Fermion System using the QUDA inverter
  /*! \ingroup invert
 *** WARNING THIS SOLVER WORKS FOR Wilson FERMIONS ONLY ***
   */
 
  class LinOpSysSolverCGQUDAWilson : public LinOpSystemSolver<LatticeFermion>
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
    LinOpSysSolverCGQUDAWilson(Handle< LinearOperator<T> > A_,
					 Handle< FermState<T,Q,Q> > state_,
					 const SysSolverCGQUDAWilsonParams& invParam_) : 
      A(A_), invParam(invParam_) 
    {
      QDPIO::cout << "LinOpSysSolverCGQUDAWilson:" << endl;
      const AnisoParam_t& aniso = invParam.WilsonParams.anisoParam;
      
     

      // These are the links
      // They may be smeared and the BC's may be applied
      links_single; links_single.resize(Nd);
      
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
      
    }

    //! Destructor is automatic
    ~LinOpSysSolverCGQUDAWilson() {}

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
      QDPIO::cout << "QUDA_CG_WILSON_SOLVER: " << res.n_count << " iterations. Rsd = " << res.resid << " Relative Rsd = " << res.resid/sqrt(norm2(chi,A->subset())) << endl;
   
      
      END_CODE();
      return res;
    }


  private:
    // Hide default constructor
    LinOpSysSolverCGQUDAWilson() {}
    
    QF links_single;
    Q links_floating;
    Handle< LinearOperator<T> > A;
    const SysSolverCGQUDAWilsonParams invParam;
    SystemSolverResults_t qudaInvert(const QF& links, 
				     const TF& chi_s,
				     TF& psi_s     
				     )const ;


  };


} // End namespace

#endif // BUILD_QUDA
#endif 


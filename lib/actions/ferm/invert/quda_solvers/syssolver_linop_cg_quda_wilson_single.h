// -*- C++ -*-
// $Id: syssolver_linop_cg_quda_wilson_single.h,v 1.1 2009-09-25 12:41:23 bjoo Exp $
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

      // These are the links
      // They may be smeared and the BC's may be applied
      links_floating; links_floating.resize(Nd);
      links_single; links_single.resize(Nd);
     
      // This whole thing doesnt work if the BC's are more sophisticated
      // Than periodic/antiperiodic.
      // If so, I can apply the boundaries a second time, to recreate the original links...
      Handle< FermBC<T,Q,Q> > fbc = state_->getFermBC();

      try { 
	SimpleFermBC<T,Q,Q>& downcast = dynamic_cast< SimpleFermBC<T,Q,Q>& >(*fbc); 

	// OK The user entered either periodic/antiperiodic/dirichlet BCs
	// Ideally I'd like to go through the boundary array to see
	// if anything other than the time is not periodic.
	// Then I'd have to Bomb
	// Currently the 'boundary' array cannot be exposed so 
	// just hope the user was not stupid.


	QDPIO::cout << "WARNING: You can only use periodic BCs in space and" << endl;
	QDPIO::cout << "  either periodic or antiperiodic BCs in time" << endl;
	QDPIO::cout << "  if you did something else, your answer may be incorrect" << endl;

	QDPIO::cout << "BC Downcast succeeded. Undoing BCs" << endl;

	// Do the BC undo in the base precision
	for(int mu=0; mu < Nd; mu++) {
	  links_floating[mu] = (state_->getLinks())[mu];
	}
	fbc->modify(links_floating);

	// Now downcast to single prec fields.
	for(int mu=0; mu < Nd; mu++) {
	  links_single[mu] = links_floating[mu];
	}
      }
      catch( std::bad_cast ) {
	QDPIO::cout << "The Boundaries selected are not periodic/antiperiodic" << endl;
	QDPIO::cout << "You cant use this solver" << endl;
	QDP_abort(1);
      }

      // OK. Links Single contains the Single Prec Links Sans BC's
      
					     
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
		       psi_s,      
		       invParam);

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
      QDPIO::cout << "CG_QUDA_WILON_SOLVER: " << res.n_count << " iterations. Rsd = " << res.resid << " Relative Rsd = " << res.resid/sqrt(norm2(chi,A->subset())) << endl;
   
      
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
				     TF& psi_s,      
				     const SysSolverCGQUDAWilsonParams& invParam) const;


  };


} // End namespace

#endif // BUILD_QUDA
#endif 


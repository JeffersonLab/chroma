// -*- C++ -*-
// $Id: syssolver_mdagm_cg.h,v 3.7 2009-06-02 15:56:40 bjoo Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by CG2
 */

#ifndef __syssolver_mdagm_cg_h__
#define __syssolver_mdagm_cg_h__
#include "chroma_config.h"

#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "lmdagm.h"
#include "actions/ferm/invert/syssolver_mdagm.h"
#include "actions/ferm/invert/syssolver_cg_params.h"
#include "actions/ferm/invert/invcg2.h"


namespace Chroma
{

  //! CG2 system solver namespace
  namespace MdagMSysSolverCGEnv
  {
    //! Register the syssolver
    bool registerAll();
  }


  //! Solve a CG2 system. Here, the operator is NOT assumed to be hermitian
  /*! \ingroup invert
   */
  template<typename T>
  class MdagMSysSolverCG : public MdagMSystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param M_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    MdagMSysSolverCG(Handle< LinearOperator<T> > A_,
		     const SysSolverCGParams& invParam_) : 
      A(A_), invParam(invParam_) 
      {}

    //! Destructor is automatic
    ~MdagMSysSolverCG() {}

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
	START_CODE();
	StopWatch swatch;
	swatch.reset(); swatch.start();

	SystemSolverResults_t res;  // initialized by a constructor
	res = InvCG2(*A, chi, psi, invParam.RsdCG, invParam.MaxCG);

#ifdef CHROMA_DO_ONE_CG_RESTART
	int n_count = res.n_count;
	res = InvCG2(*A, chi, psi, invParam.RsdCGRestart, invParam.MaxCGRestart);
	res.n_count += n_count;
#endif 
	{ // Find true residuum
	  T tmp=zero;
	  T r=zero;
	  (*A)(tmp,psi, PLUS);
	  (*A)(r,tmp, MINUS);
	  r[A->subset()] -= chi;
	  res.resid = sqrt(norm2(r,A->subset()));
	}
	
	swatch.stop();
	QDPIO::cout << "CG_SOLVER: " << res.n_count 
		    << " iterations. Rsd = " << res.resid 
		    << " Relative Rsd = " << res.resid/sqrt(norm2(chi,A->subset())) << endl;
	
	double time = swatch.getTimeInSeconds();
	QDPIO::cout << "CG_SOLVER_TIME: "<<time<< " sec" << endl;
	

	END_CODE();

	return res;
      }


    //! Solve the linear system starting with a chrono guess 
    /*! 
     * \param psi solution (Write)
     * \param chi source   (Read)
     * \param predictor   a chronological predictor (Read)
     * \return syssolver results
     */

    SystemSolverResults_t operator()(T& psi, const T& chi, 
				     AbsChronologicalPredictor4D<T>& predictor) const 
    {
      
      START_CODE();

      // This solver uses InvCG2, so A is just the matrix.
      // I need to predict with A^\dagger A
      {
	Handle< LinearOperator<T> > MdagM( new MdagMLinOp<T>(A) );
	predictor(psi, (*MdagM), chi);
      }
      // Do solve
      SystemSolverResults_t res=(*this)(psi,chi);

      // Store result
      predictor.newVector(psi);
      END_CODE();
      return res;
    }

  private:
    // Hide default constructor
    MdagMSysSolverCG() {}

    Handle< LinearOperator<T> > A;
    SysSolverCGParams invParam;
  };


} // End namespace

#endif 


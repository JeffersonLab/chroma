// -*- C++ -*-
// $Id: syssolver_mdagm_mr.h,v 1.3 2009-06-02 15:56:40 bjoo Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by MR
 */

#ifndef __syssolver_mdagm_mr_h__
#define __syssolver_mdagm_mr_h__

#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "lmdagm.h"
#include "actions/ferm/invert/syssolver_mdagm.h"
#include "actions/ferm/invert/syssolver_mr_params.h"
#include "actions/ferm/invert/invmr.h"

namespace Chroma
{

  //! MR system solver namespace
  namespace MdagMSysSolverMREnv
  {
    //! Register the syssolver
    bool registerAll();
  }


  //! Solve a MR system. Here, the operator is NOT assumed to be hermitian
  /*! \ingroup invert
   */
  template<typename T>
  class MdagMSysSolverMR : public MdagMSystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param M_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    MdagMSysSolverMR(Handle< LinearOperator<T> > A_,
		     const SysSolverMRParams& invParam_) : 
      A(A_), invParam(invParam_) 
      {}

    //! Destructor is automatic
    ~MdagMSysSolverMR() {}

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
	SystemSolverResults_t res1,res2;  // initialized by a constructor
	T tmp;

	swatch.reset(); swatch.start();
	res1 = InvMR(*A, chi, tmp, invParam.RsdMR, invParam.MaxMR,MINUS);
	res2 = InvMR(*A, tmp, psi, invParam.RsdMR, invParam.MaxMR,PLUS);

	res2.n_count += res1.n_count;
	{ // Find true residuum
	  tmp=zero;
	  T re=zero;
	  (*A)(tmp, psi, PLUS);
	  (*A)(re,tmp, MINUS);
	  re[A->subset()] -= chi;
	  res2.resid = sqrt(norm2(re,A->subset()));
	}

	swatch.stop();
	QDPIO::cout << "MR_SOLVER: " << res2.n_count 
		    << " iterations. Rsd = " << res2.resid 
		    << " Relative Rsd = " << res2.resid/sqrt(norm2(chi,A->subset())) << endl;

	double time = swatch.getTimeInSeconds();
	QDPIO::cout << "MR_SOLVER_TIME: "<<time<< " sec" << endl;
	
	END_CODE();

	return res;
      }


    //! Solver the linear system
    /*!
     * \param psi      solution ( Modify )
     * \param chi      source ( Read )
     * \return syssolver results
     */
    SystemSolverResults_t operator() (T& psi, const T& chi,
				      AbsChronologicalPredictor4D<T>& predictor) const
      {
	START_CODE();
	StopWatch swatch;
	SystemSolverResults_t res1,res2;  // initialized by a constructor
	T tmp=zero;
	
	swatch.reset(); swatch.start();

	try { 
	  AbsTwoStepChronologicalPredictor4D<T>& two_step_predictor = 
	    dynamic_cast<AbsTwoStepChronologicalPredictor4D<T>&>(predictor);

	  two_step_predictor.predictY(tmp, *A, chi);

	  res1 = InvMR(*A, chi, tmp, invParam.RsdMR, invParam.MaxMR,MINUS);

	  two_step_predictor.newYVector(tmp);

	  two_step_predictor.predictX(psi, *A, tmp);

	  res2 = InvMR(*A, tmp, psi, invParam.RsdMR, invParam.MaxMR,PLUS);
	  
	  two_step_predictor.newXVector(psi);


	}
	catch(std::bad_cast) {
	  T Y ;
	  Handle< LinearOperator<T> > MdagM( new MdagMLinOp<T>(A) );
	  predictor(psi, (*MdagM), chi);
	  (*A)(Y, psi, PLUS); // Y = M X
	  res1 = InvMR(*A, chi, tmp, invParam.RsdMR, invParam.MaxMR,MINUS);
	  res2 = InvMR(*A, tmp, psi, invParam.RsdMR, invParam.MaxMR,PLUS);
	  predictor.newVector(psi);
	}

	res2.n_count += res1.n_count;

	{ // Find true residuum
	  tmp=zero;
	  T re=zero;
	  (*A)(tmp, psi, PLUS);
	  (*A)(re,tmp, MINUS);
	  re[A->subset()] -= chi;
	  res2.resid = sqrt(norm2(re,A->subset()));
	}
	  
	swatch.stop();
	QDPIO::cout << "MR_SOLVER: " << res2.n_count 
		    << " iterations. Rsd = " << res2.resid 
		    << " Relative Rsd = " << res2.resid/sqrt(norm2(chi,A->subset())) << endl;
	
	double time = swatch.getTimeInSeconds();
	QDPIO::cout << "MR_SOLVER_TIME: "<<time<< " sec" << endl;
	
	END_CODE();
	
	return res;
      }


  private:
    // Hide default constructor
    MdagMSysSolverMR() {}

    Handle< LinearOperator<T> > A;
    SysSolverMRParams invParam;
  };


} // End namespace

#endif 


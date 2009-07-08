// -*- C++ -*-
// $Id: syssolver_mdagm_ibicgstab.h,v 3.2 2009-07-08 18:46:47 bjoo Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by BiCGStab
 */

#ifndef __syssolver_mdagm_ibicgstab_h__
#define __syssolver_mdagm_ibicgstab_h__
#include "chroma_config.h"

#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/invert/syssolver_mdagm.h"
#include "actions/ferm/invert/syssolver_bicgstab_params.h"
#include "actions/ferm/invert/invibicgstab.h"
#ifdef CHROMA_DO_ONE_CG_RESTART
#include "actions/ferm/invert/invcg2.h"
#endif

#include "lmdagm.h"
#include "update/molecdyn/predictor/chrono_predictor.h"
#include "update/molecdyn/predictor/zero_guess_predictor.h"
namespace Chroma
{

  //! IBiCGStab system solver namespace
  namespace MdagMSysSolverIBiCGStabEnv
  {
    //! Register the syssolver
    bool registerAll();
  }


  //! Solve a IBiCGStab system. Here, the operator is NOT assumed to be hermitian
  /*! \ingroup invert
   */
  template<typename T>
  class MdagMSysSolverIBiCGStab : public MdagMSystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param M_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    MdagMSysSolverIBiCGStab(Handle< LinearOperator<T> > A_,
			   const SysSolverBiCGStabParams& invParam_) : 
      A(A_), invParam(invParam_) 
    {}

    //! Destructor is automatic
    ~MdagMSysSolverIBiCGStab() {}

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
      SystemSolverResults_t res1,res2,res3;  // initialized by a constructo
      swatch.reset(); swatch.start();

      T Y;
      Handle< LinearOperator<T> > MdagM( new MdagMLinOp<T>(A) );
      (*A)(Y, psi, PLUS); // Y = M X
      
      res1 = InvIBiCGStab(*A, chi, Y, invParam.RsdBiCGStab, invParam.MaxBiCGStab, MINUS );
	
      // Step 2: Solve M X = Y
      res2 = InvIBiCGStab(*A, Y, psi, invParam.RsdBiCGStab, invParam.MaxBiCGStab, PLUS );

      res3.n_count = 0;
      // Potential safety polishup
#ifdef CHROMA_DO_ONE_CG_RESTART
      // CG Polish - should be very quick
      res3 = InvCG2(*A, chi, psi, invParam.RsdBiCGStab, invParam.MaxBiCGStab);
#endif
      res3.n_count += res2.n_count + res1.n_count;

      { // Find true residuum
	Y=zero;
	T re=zero;
	(*A)(Y, psi, PLUS);
	(*A)(re,Y, MINUS);
	re[A->subset()] -= chi;
	res3.resid = sqrt(norm2(re,A->subset()));
      }

      swatch.stop();
      QDPIO::cout << "IBICGSTAB_SOLVER: " << res3.n_count 
		  << " iterations. Rsd = " << res3.resid 
		  << " Relative Rsd = " << res3.resid/sqrt(norm2(chi,A->subset())) << endl;
      
      double time = swatch.getTimeInSeconds();
      QDPIO::cout << "IBICGSTAB_SOLVER_TIME: "<<time<< " sec" << endl;
	
      
      END_CODE();

      return res3;

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
	SystemSolverResults_t res1,res2,res3;  // initialized by a constructo
	swatch.reset(); swatch.start();

	T Y = psi;

	try { 
	  // Get a two step solution plan
	  AbsTwoStepChronologicalPredictor4D<T>& two_step_predictor
	    = dynamic_cast<AbsTwoStepChronologicalPredictor4D<T>& >(predictor);

	  // Hooray , we succeeded.
	  // Step 1: Solve M^\dagger Y = chi

	  two_step_predictor.predictY(Y,*A,chi);
	  res1 = InvIBiCGStab(*A, chi, Y, invParam.RsdBiCGStab, invParam.MaxBiCGStab, MINUS );
	  two_step_predictor.newYVector(Y);

	  // Step 2: Solve M X = Y
	  Handle<LinearOperator<T> > MdagM(new MdagMLinOp<T>(A));
	  two_step_predictor.predictX(psi,*MdagM, chi);
	  res2 = InvIBiCGStab(*A, Y, psi, invParam.RsdBiCGStab, invParam.MaxBiCGStab, PLUS );
	  two_step_predictor.newXVector(psi);
	}
	catch(std::bad_cast) {

	  // Boo Hiss, we can't downcast to a two step predictor.
	  // We rely on the fact that we can predict 
	  //    X ~ (M^\dagger M)^{-1} chi
	  // and then 
	  //    X ~  M^{-1} M^{-\dagger} chi
	  //
	  //  Then MX ~ M^{-\dagger} chi ~ Y

	  T Y ;
	  Handle< LinearOperator<T> > MdagM( new MdagMLinOp<T>(A) );
	  predictor(psi, (*MdagM), chi);
	  (*A)(Y, psi, PLUS); // Y = M X
	  res1 = InvIBiCGStab(*A, chi, Y, invParam.RsdBiCGStab, invParam.MaxBiCGStab, MINUS );
	  // Step 2: Solve M X = Y
	  res2 = InvIBiCGStab(*A, Y, psi, invParam.RsdBiCGStab, invParam.MaxBiCGStab, PLUS );

	  predictor.newVector(psi);

	}

	res3.n_count = 0;
	// Potential safety polishup
#ifdef CHROMA_DO_ONE_CG_RESTART
	// CG Polish - should be very quick
	res3 = InvCG2(*A, chi, psi, invParam.RsdBiCGStab, invParam.MaxBiCGStab);
#endif

	res3.n_count += res2.n_count + res1.n_count;

	{ // Find true residuum
	  Y=zero;
	  T re=zero;
	  (*A)(Y, psi, PLUS);
	  (*A)(re,Y, MINUS);
	  re[A->subset()] -= chi;
	  res3.resid = sqrt(norm2(re,A->subset()));
	}

	swatch.stop();
	QDPIO::cout << "IBICGSTAB_SOLVER: " << res3.n_count 
		    << " iterations. Rsd = " << res3.resid 
		    << " Relative Rsd = " << res3.resid/sqrt(norm2(chi,A->subset())) << endl;

	double time = swatch.getTimeInSeconds();
	QDPIO::cout << "IBICGSTAB_SOLVER_TIME: "<<time<< " sec" << endl;
	

	END_CODE();

	return res3;
    }



  private:
    // Hide default constructor
    MdagMSysSolverIBiCGStab() {}

    Handle< LinearOperator<T> > A;
    SysSolverBiCGStabParams invParam;
  };


} // End namespace

#endif 


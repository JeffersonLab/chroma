// -*- C++ -*-
// $Id: syssolver_mdagm_rel_ibicgstab_clover.h,v 3.1 2009-07-06 19:02:34 bjoo Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by IBiCGStab
 */

#ifndef __syssolver_mdagm_rel_ibicgstab_multiprec_h__
#define __syssolver_mdagm_rel_ibicgstab_multiprec_h__
#include "chroma_config.h"

#include "handle.h"
#include "state.h"
#include "syssolver.h"
#include "linearop.h"
#include "lmdagm.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/invert/syssolver_mdagm.h"
#include "actions/ferm/invert/reliable_ibicgstab.h"
#include "actions/ferm/invert/syssolver_mdagm_factory.h"
#include "actions/ferm/invert/syssolver_rel_bicgstab_clover_params.h"
#include "actions/ferm/linop/eoprec_clover_dumb_linop_w.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"

#include <string>
using namespace std;

namespace Chroma
{

  //! Richardson system solver namespace
  namespace MdagMSysSolverReliableIBiCGStabCloverEnv
  {
    //! Register the syssolver
    bool registerAll();
  }



  //! Solve a system using Richardson iteration.
  /*! \ingroup invert
 *** WARNING THIS SOLVER WORKS FOR CLOVER FERMIONS ONLY ***
   */
 
  class MdagMSysSolverReliableIBiCGStabClover : public MdagMSystemSolver<LatticeFermion>
  {
  public:
    typedef LatticeFermion T;
    typedef LatticeColorMatrix U;
    typedef multi1d<LatticeColorMatrix> Q;
 
    typedef LatticeFermionF TF;
    typedef LatticeColorMatrixF UF;
    typedef multi1d<LatticeColorMatrixF> QF;

    typedef LatticeFermionD TD;
    typedef LatticeColorMatrixD UD;
    typedef multi1d<LatticeColorMatrixD> QD;

    //! Constructor
    /*!
     * \param M_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    MdagMSysSolverReliableIBiCGStabClover(Handle< LinearOperator<T> > A_,
					 Handle< FermState<T,Q,Q> > state_,
					 const SysSolverReliableBiCGStabCloverParams& invParam_) : 
      A(A_), invParam(invParam_) 
    {

      // Get the Links out of the state and convertto single.
      QF links_single; links_single.resize(Nd);
      QD links_double; links_double.resize(Nd);

      const Q& links = state_->getLinks();
      for(int mu=0; mu < Nd; mu++) { 
	links_single[mu] = links[mu];
	links_double[mu] = links[mu];
      }

      
      // Links single hold the possibly stouted links
      // with gaugeBCs applied... 
      // Now I need to create a single prec state...
      fstate_single = new PeriodicFermState<TF,QF,QF>(links_single);
      fstate_double = new PeriodicFermState<TD,QD,QD>(links_double);

      // Make single precision M
      M_single= new EvenOddPrecDumbCloverFLinOp( fstate_single, invParam_.clovParams );
      M_double= new EvenOddPrecDumbCloverDLinOp( fstate_double, invParam_.clovParams );

      
					     
    }

    //! Destructor is automatic
    ~MdagMSysSolverReliableIBiCGStabClover() {}

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
      SystemSolverResults_t res,res1,res2;

      START_CODE();
      StopWatch swatch;
      swatch.start();
      const Subset& s = (*M_double).subset();
      
      
      // Boo Hiss, we can't downcast to a two step predictor.
      // We rely on the fact that we can predict 
      //    X ~ (M^\dagger M)^{-1} chi
      // and then 
      //    X ~  M^{-1} M^{-\dagger} chi
      //
      //  Then MX ~ M^{-\dagger} chi ~ Y
	
      T Y ;
      Handle< LinearOperator<T> > MdagM( new MdagMLinOp<T>(A) );
      (*A)(Y, psi, PLUS); // Y = M X
	  
      TD psi_d;
      TD chi_d;
      psi_d[s] = Y;
      chi_d[s] = chi;
      // Two Step IBiCGStab:
      // Step 1:  M^\dagger Y = chi;
      
      res1=InvIBiCGStabReliable(*M_double,
			       *M_single,
			       chi_d,
			       psi_d,
			       invParam.RsdTarget,
			       invParam.Delta,
			       invParam.MaxIter,
			       MINUS);

      chi_d[s] = psi;

      res2=InvIBiCGStabReliable(*M_double,
			       *M_single,
			       psi_d,
			       chi_d,
			       invParam.RsdTarget,
			       invParam.Delta,
			       invParam.MaxIter,
			       PLUS);
      
      psi[s] = chi_d;
	
      // Check solution
      { 
	T r;
	r[A->subset()]=chi;
	T tmp,tmp2;
	// M^\dag M
	(*A)(tmp, psi, PLUS);
	(*A)(tmp2,tmp, MINUS);
	r[A->subset()] -= tmp2;
	res.n_count = res1.n_count + res2.n_count;
	res.resid = sqrt(norm2(r, A->subset()));
      }
      QDPIO::cout << "RELIABLE_IBICGSTAB_SOLVER: " << res.n_count 
		  << " iterations. Rsd = " << res.resid 
		  << " Relative Rsd = " << res.resid/sqrt(norm2(chi,A->subset())) << endl;
      
      swatch.stop();
      double time = swatch.getTimeInSeconds();
      QDPIO::cout << "RELIABLE_IBICGSTAB_SOLVER_TIME: "<<time<< " sec" << endl;
      
      
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
      SystemSolverResults_t res,res1,res2;

      START_CODE();
      StopWatch swatch;
      swatch.start();
      const Subset& s = (*M_double).subset();

      TD psi_d; TD chi_d;
      try { 
	// Get a two step solution plan
	AbsTwoStepChronologicalPredictor4D<T>& two_step_predictor
	  = dynamic_cast<AbsTwoStepChronologicalPredictor4D<T>& >(predictor);
      
        // Hooray , we succeeded.
        // Step 1: Solve M^\dagger Y = chi
        T Y = zero;
	two_step_predictor.predictY(Y,*A,chi);
	psi_d[s] = Y;
	chi_d[s] = chi;
	// Two Step IBiCGStab:
	// Step 1:  M^\dagger Y = chi;
	  
	res1=InvIBiCGStabReliable(*M_double,
				 *M_single,
				 chi_d,
				 psi_d,
				 invParam.RsdTarget,
				 invParam.Delta,
				 invParam.MaxIter,
				 MINUS);

	Y[s] = psi_d;
	two_step_predictor.newYVector(Y);
	
	Handle<LinearOperator<T> > MdagM(new MdagMLinOp<T>(A));

	two_step_predictor.predictX(psi,*MdagM, chi);
	chi_d[s] = psi;

	res2=InvIBiCGStabReliable(*M_double,
				 *M_single,
				 psi_d,
				 chi_d,
				 invParam.RsdTarget,
				 invParam.Delta,
				 invParam.MaxIter,
				 PLUS);
	

	psi[s] = chi_d;
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
	  
	psi_d[s] = Y;
	chi_d[s] = chi;
	// Two Step IBiCGStab:
	// Step 1:  M^\dagger Y = chi;
	  
	res1=InvIBiCGStabReliable(*M_double,
				 *M_single,
				 chi_d,
				 psi_d,
				 invParam.RsdTarget,
				 invParam.Delta,
				 invParam.MaxIter,
				 MINUS);

	chi_d[s] = psi;

	res2=InvIBiCGStabReliable(*M_double,
				 *M_single,
				 psi_d,
				 chi_d,
				 invParam.RsdTarget,
				 invParam.Delta,
				 invParam.MaxIter,
				 PLUS);

	psi[s] = chi_d;
	predictor.newVector(psi);

	
      }
	
      // Check solution
      { 
	T r;
	r[A->subset()]=chi;
	T tmp,tmp2;
	// M^\dag M
	(*A)(tmp, psi, PLUS);
	(*A)(tmp2,tmp, MINUS);
	r[A->subset()] -= tmp2;
	res.n_count = res1.n_count + res2.n_count;
	res.resid = sqrt(norm2(r, A->subset()));
      }
      QDPIO::cout << "RELIABLE_IBICGSTAB_SOLVER: " << res.n_count 
		  << " iterations. Rsd = " << res.resid 
		  << " Relative Rsd = " << res.resid/sqrt(norm2(chi,A->subset())) << endl;
      
      swatch.stop();
      double time = swatch.getTimeInSeconds();
      QDPIO::cout << "RELIABLE_IBICGSTAB_SOLVER_TIME: "<<time<< " sec" << endl;
      
      
      END_CODE();
      return res;
    }
    

  private:
    // Hide default constructor
    MdagMSysSolverReliableIBiCGStabClover() {}
    Handle< LinearOperator<T> > A;
    const SysSolverReliableBiCGStabCloverParams invParam;

    // Created and initialized here.
    Handle< FermState<TF, QF, QF> > fstate_single;
    Handle< FermState<TD, QD, QD> > fstate_double;
    Handle< LinearOperator<TF> > M_single;
    Handle< LinearOperator<TD> > M_double;
  };


} // End namespace

#endif 


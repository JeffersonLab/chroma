// -*- C++ -*-
// $Id: syssolver_mdagm_rel_cg_clover.h,v 3.4 2009-06-02 15:56:40 bjoo Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by BiCGStab
 */

#ifndef __syssolver_mdagm_rel_cg_multiprec_h__
#define __syssolver_mdagm_rel_cg_multiprec_h__
#include "chroma_config.h"

#include "handle.h"
#include "state.h"
#include "syssolver.h"
#include "linearop.h"
#include "lmdagm.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/invert/syssolver_mdagm.h"
#include "actions/ferm/invert/reliable_cg.h"
#include "actions/ferm/invert/syssolver_mdagm_factory.h"
#include "actions/ferm/invert/syssolver_rel_bicgstab_clover_params.h"
#include "actions/ferm/linop/eoprec_clover_dumb_linop_w.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"

#include <string>
using namespace std;

namespace Chroma
{

  //! Richardson system solver namespace
  namespace MdagMSysSolverReliableCGCloverEnv
  {
    //! Register the syssolver
    bool registerAll();
  }



  //! Solve a system using Richardson iteration.
  /*! \ingroup invert
 *** WARNING THIS SOLVER WORKS FOR CLOVER FERMIONS ONLY ***
   */
 
  class MdagMSysSolverReliableCGClover : public MdagMSystemSolver<LatticeFermion>
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
    MdagMSysSolverReliableCGClover(Handle< LinearOperator<T> > A_,
					 Handle< FermState<T,Q,Q> > state_,
					 const SysSolverReliableCGCloverParams& invParam_) : 
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
    ~MdagMSysSolverReliableCGClover() {}

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
      const Subset& s = (*M_double).subset();
      //    T MdagChi;

      // This is a CGNE. So create new RHS
      //      (*A)(MdagChi, chi, MINUS);
      // Handle< LinearOperator<T> > MM(new MdagMLinOp<T>(A));

      TD psi_d; psi_d[s] = psi;
      TD chi_d; chi_d[s] = chi;
      

      // Two Step CG:
      // Step 1:  M^\dagger Y = chi;

      res=InvCGReliable(*M_double,
			*M_single,
			chi_d,
			psi_d,
			invParam.RsdTarget,
			invParam.Delta,
			invParam.MaxIter);
      psi = psi_d;
      
      { 
	T r;
	r[A->subset()]=chi;
	T tmp,tmp2;
	// M^\dag M
	(*A)(tmp, psi, PLUS);
	(*A)(tmp2,tmp, MINUS);
	r[A->subset()] -= tmp2;
	res.resid = sqrt(norm2(r, A->subset()));
      }
      QDPIO::cout << "RELIABLE_CG_SOLVER: " << res.n_count << " iterations. Rsd = " << res.resid << " Relative Rsd = " << res.resid/sqrt(norm2(chi,A->subset())) << endl;
   
      
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
    MdagMSysSolverReliableCGClover() {}
    Handle< LinearOperator<T> > A;
    const SysSolverReliableCGCloverParams invParam;

    // Created and initialized here.
    Handle< FermState<TF, QF, QF> > fstate_single;
    Handle< FermState<TD, QD, QD> > fstate_double;
    Handle< LinearOperator<TF> > M_single;
    Handle< LinearOperator<TD> > M_double;
  };


} // End namespace

#endif 


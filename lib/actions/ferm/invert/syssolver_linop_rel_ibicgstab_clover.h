// -*- C++ -*-
// $Id: syssolver_linop_rel_ibicgstab_clover.h,v 3.1 2009-07-06 19:02:34 bjoo Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by BiCGStab
 */

#ifndef __syssolver_linop_rel_ibicgstab_multiprec_h__
#define __syssolver_linop_rel_ibicgstab_multiprec_h__
#include "chroma_config.h"

#include "handle.h"
#include "state.h"
#include "syssolver.h"
#include "linearop.h"
#include "lmdagm.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/invert/syssolver_mdagm.h"
#include "actions/ferm/invert/reliable_ibicgstab.h"
#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_rel_bicgstab_clover_params.h"
#include "actions/ferm/linop/eoprec_clover_dumb_linop_w.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"

#include <string>
using namespace std;

namespace Chroma
{

  //! Richardson system solver namespace
  namespace LinOpSysSolverReliableIBiCGStabCloverEnv
  {
    //! Register the syssolver
    bool registerAll();
  }



  //! Solve a system using Richardson iteration.
  /*! \ingroup invert
 *** WARNING THIS SOLVER WORKS FOR CLOVER FERMIONS ONLY ***
   */
 
  class LinOpSysSolverReliableIBiCGStabClover : public LinOpSystemSolver<LatticeFermion>
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
    LinOpSysSolverReliableIBiCGStabClover(Handle< LinearOperator<T> > A_,
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
    ~LinOpSysSolverReliableIBiCGStabClover() {}

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

      TD psi_d = psi;
      TD chi_d = chi;

      res=InvIBiCGStabReliable(*M_double,
			  *M_single,
			  chi_d,
			  psi_d,
			  invParam.RsdTarget,
			  invParam.Delta,
			  invParam.MaxIter,
			  PLUS);
      
      psi = psi_d;

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
      QDPIO::cout << "RELIABLE_IBICGSTAB_SOLVER: " << res.n_count << " iterations. Rsd = " << res.resid << " Relative Rsd = " << res.resid/sqrt(norm2(chi,A->subset())) << endl;
      QDPIO::cout << "RELIABLE_IBICGSTAB_SOLVER_TIME: "<<time<< " sec" << endl;
   
      
      END_CODE();
      return res;
    }


  private:
    // Hide default constructor
    LinOpSysSolverReliableIBiCGStabClover() {}
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


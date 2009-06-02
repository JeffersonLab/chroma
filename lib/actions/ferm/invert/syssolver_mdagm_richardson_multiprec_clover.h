// -*- C++ -*-
// $Id: syssolver_mdagm_richardson_multiprec_clover.h,v 3.2 2009-06-02 15:56:40 bjoo Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by BiCGStab
 */

#ifndef __syssolver_mdagm_richardson_multiprec_h__
#define __syssolver_mdadm_richardson_multiprec_h__
#include "chroma_config.h"

#include "handle.h"
#include "state.h"
#include "syssolver.h"
#include "linearop.h"
#include "lmdagm.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/invert/syssolver_mdagm.h"
#include "actions/ferm/invert/inv_multiprec_richardson.h"
#include "actions/ferm/invert/syssolver_mdagm_factory.h"
#include "actions/ferm/invert/syssolver_richardson_clover_params.h"
#include "actions/ferm/linop/eoprec_clover_dumb_linop_w.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"

#include <string>
using namespace std;

namespace Chroma
{

  //! Richardson system solver namespace
  namespace MdagMSysSolverRichardsonCloverEnv
  {
    //! Register the syssolver
    bool registerAll();
  }



  //! Solve a system using Richardson iteration.
  /*! \ingroup invert
   */
 
  class MdagMSysSolverRichardsonClover : public MdagMSystemSolver<LatticeFermion>
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
    MdagMSysSolverRichardsonClover(Handle< LinearOperator<T> > A_,
			     Handle< FermState<T,Q,Q> > state_,
			   const SysSolverRichardsonCloverParams& invParam_) : 
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

      std::istringstream is( invParam_.innerSolverParams.xml );
      XMLReader paramtop(is);

      DInv = TheMdagMFermFSystemSolverFactory::Instance().createObject( invParam_.innerSolverParams.id, paramtop, 
									invParam_.innerSolverParams.path,
									fstate_single, 
									M_single);
      
					     
    }

    //! Destructor is automatic
    ~MdagMSysSolverRichardsonClover() {}

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
      Handle< LinearOperator<TD> > MM_double( new MdagMLinOp<TD>(M_double));

      TD psi_d = psi;
      TD chi_d = chi;

      InvMultiPrecRichardson(*DInv, 
			     *MM_double,
			     chi_d,
			     psi_d,
			     invParam.MaxIter,
			     invParam.RsdTarget,
			     res);
      
      psi = psi_d;
      swatch.stop();
      double time = swatch.getTimeInSeconds();
      { 
	T r;
	r[A->subset()]=chi;
	T tmp, tmp1;
	(*A)(tmp, psi, PLUS);
	(*A)(tmp1, tmp, MINUS);
	r[A->subset()] -= tmp1;
	res.resid = sqrt(norm2(r, A->subset())/norm2(chi, A->subset()));
      }

      QDPIO::cout << "MULTIPREC_RICHARDSON_SOLVER: " << res.n_count << " iterations. Rsd = " << res.resid << " Relative Rsd = " << res.resid/sqrt(norm2(chi,A->subset())) << endl;
      QDPIO::cout << "MULTIPREC_RICHARDSON_SOLVER_TIME: "<<time<< " sec" << endl;

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
    MdagMSysSolverRichardsonClover() {}
    Handle< LinearOperator<T> > A;
    const SysSolverRichardsonCloverParams invParam;

    // Created and initialized here.
    Handle< FermState<TF, QF, QF> > fstate_single;
    Handle< FermState<TD, QD, QD> > fstate_double;
    Handle< LinearOperator<TF> > M_single;
    Handle< LinearOperator<TD> > M_double;
    Handle< MdagMSystemSolver< TF > > DInv;


  };


} // End namespace

#endif 


// $Id: syssolver_MdagM_qdp_mg.cc, v1.0 2013-06-20 22:12 sdcohen $
/*! \file
 *  \brief Make contact with the QDP clover multigrid solver, transfer
 *         the gauge field, generate the coarse grids, solve systems
 */
#include "state.h"
#include "meas/inline/io/named_objmap.h"
#include "actions/ferm/invert/syssolver_mdagm_factory.h"
#include "actions/ferm/invert/syssolver_mdagm_aggregate.h"


#include "actions/ferm/invert/qop_mg/syssolver_linop_qop_mg_w.h"

#include "meas/glue/mesplq.h"

namespace Chroma
{
  
  template<typename T> // T is the Lattice Fermion type
  MdagMSysSolverQOPMG<T>::
  MdagMSysSolverQOPMG(Handle< LinearOperator<T> > A_,
		      Handle< FermState<T,Q,Q> > state_, 
		      const SysSolverQOPMGParams& invParam_) : 
    A(A_), state(state_), linop_solver(new LinOpSysSolverQOPMG<T>(A_,state_,invParam_)), invParam(invParam_) {}
				    
				       
  
  template<typename T> // T is the Lattice Fermion type
  MdagMSysSolverQOPMG<T>::~MdagMSysSolverQOPMG(){}
  
  //! Solve the linear system
  /*!
   * \param psi      solution ( Modify )
   * \param chi      source ( Read )
   * \return syssolver results
   */
  template<typename T> // T is the Lattice Fermion type
  SystemSolverResults_t
  MdagMSysSolverQOPMG<T>::operator() (T& psi, const T& chi) const
  {
    START_CODE();
    
    StopWatch swatch;
    SystemSolverResults_t res1,res2,res3;  // initialized by a constructo
    swatch.reset(); swatch.start();
    
    T Y;
    // psi is some initial guess for X, so then 
    // we turn that into an initial guess for Y
    (*A)(Y, psi, PLUS); // Y = M X
    res1=(*linop_solver)(Y,chi,MINUS);  // solve Y = M^{-dagger} chi
    res2=(*linop_solver)(psi,Y,PLUS);   // solve psi = M^{-1} Y = M^{-1} M^{-dagger} chi
    { // Find true residuum
      Y=zero;
      T re=zero;
      (*A)(Y, psi, PLUS);
      (*A)(re,Y, MINUS);
      re[A->subset()] -= chi;
      res3.resid = sqrt(norm2(re,A->subset()));
    }
    res3.n_count = res2.n_count + res1.n_count;

    swatch.stop();
    QDPIO::cout << "QOPMG_SOLVER: " << res3.n_count 
		<< " iterations. Rsd = " << res3.resid 
		<< " Relative Rsd = " << res3.resid/sqrt(norm2(chi,A->subset())) << endl;
     
    double time = swatch.getTimeInSeconds();
    QDPIO::cout << "QOMPG_SOLVER_TIME: "<<time<< " sec" << endl;
     
    END_CODE();
    return res3;
  }


  //! Solve the linear system
  /*!
   * \param psi      solution ( Modify )
   * \param chi      source ( Read )
   * \param isign    solve with dagger or not
   * \return syssolver results
   */
  template<typename T> // T is the Lattice Fermion type
  SystemSolverResults_t
  MdagMSysSolverQOPMG<T>::operator() (T& psi, const T& chi, AbsChronologicalPredictor4D<T>& predictor) const
  {
    // Ignore chrono for now.
    (*this)(psi,chi);
  }


  //! QDP multigrid system solver namespace
  namespace MdagMSysSolverQOPMGEnv
  {
    //! Anonymous namespace
    namespace
    {
      //! Name to be used
      const std::string name("QOP_CLOVER_MULTIGRID_INVERTER");
      
      //! Local registration flag
      bool registered = false;
    }
    

    //! Callback function for standard precision
    MdagMSystemSolver<LatticeFermion>*
    createFerm( XMLReader& xml_in,
		const std::string& path,
		Handle< FermState< LatticeFermion, 
		multi1d<LatticeColorMatrix>,
		multi1d<LatticeColorMatrix> > > state, 
    Handle< LinearOperator<LatticeFermion> >           A)
    {
      return new MdagMSysSolverQOPMG<LatticeFermion>
	(A, state, SysSolverQOPMGParams(xml_in, path));
    }

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
	{
	  success &= Chroma::TheMdagMFermSystemSolverFactory::Instance().registerObject(name, createFerm);
	  registered = true;
	}
      return success;
    }
  }
}

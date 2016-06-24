/*! \file
 *  \brief Make contact with the QDP clover multigrid solver, transfer
 *         the gauge field, generate the coarse grids, solve systems
 */
#include "state.h"
#include "meas/inline/io/named_objmap.h"
#include "actions/ferm/invert/syssolver_mdagm_factory.h"
#include "actions/ferm/invert/syssolver_mdagm_aggregate.h"
#include <cstdio>
#include <ostream>

#include "actions/ferm/invert/qop_mg/syssolver_linop_qop_mg_w.h"
#include "actions/ferm/invert/qop_mg/syssolver_mdagm_qop_mg_w.h"
//Added support for MG predictor.
#include "update/molecdyn/predictor/null_predictor.h"
#include "actions/ferm/linop/lunprec_w.h"

#include "meas/glue/mesplq.h"

#if BASE_PRECISION == 32
#define QDP_Precision 'F'
#define QLA_Precision 'F'
#define toReal toFloat
#elif BASE_PRECISION == 64
#define QDP_Precision 'D'
#define QLA_Precision 'D'
#define toReal toDouble
#endif

extern "C" {
  // This should be placed on the include path.
#include "wilsonmg-interface.h"

}

namespace Chroma
{

  namespace MGMdagMInternal { 
    
    /*! This will remap the MdagM params for the linop. For example,
     * we will turn off the TerminateOnFail feature, so that this outer solver
     * can control termination
     */
    SysSolverQOPMGParams remapParams(const SysSolverQOPMGParams& invParam_)
    {
      SysSolverQOPMGParams ret_val = invParam_;
      ret_val.TerminateOnFail = false; // Turn off terminate on fail as we 
				       // will retry from MdagM

      ret_val.MaxIter = invParam_.MaxIter/2; // MdagM is 2 solves, so I want a maxIter that is for 1 solve, which hopefully should be about 1/2 of MdagM
      return ret_val;
    }
  };

  // This will come from syssolver_linop_qop_mg_w.h
  //  static multi1d<LatticeColorMatrix> u;
  // These functions will allow QDP to look into the Chroma gauge field and set
  // the QDP gauge field at each site equal to the one in Chroma. There doesn't
  // seem to be a good way to treat the extra std::vector index of the gauge field,
  
  MdagMSysSolverQOPMG::MdagMSysSolverQOPMG(Handle< LinearOperator<LatticeFermion> > A_,
			Handle< FermState<LatticeFermion,multi1d<LatticeColorMatrix>,multi1d<LatticeColorMatrix > > > state_, 
					   const SysSolverQOPMGParams& invParam_) : 
    A(A_), 
    state(state_), 
    invParam(invParam_),
    Dinv(new LinOpSysSolverQOPMG(A_,state_,MGMdagMInternal::remapParams(invParam_)))
  {            
    QDPIO::cout<<"MdagM multigrid initialized"<<std::endl;
  }
  
  MdagMSysSolverQOPMG::~MdagMSysSolverQOPMG()
  {
    //if (invParam.Levels<0) MGP(finalize)();
  }
  
  
  //! Solve the linear system
  /*!
   * \param psi      solution ( Modify )
   * \param chi      source ( Read )
   * \return syssolver results
   */
  SystemSolverResults_t
  MdagMSysSolverQOPMG::operator() (LatticeFermion& psi, const LatticeFermion& chi,  AbsChronologicalPredictor4D<LatticeFermion>& predictor) const
  {
    START_CODE();
    typedef LatticeFermion T;
    typedef multi1d<LatticeColorMatrix> P;
    typedef multi1d<LatticeColorMatrix> Q;
    SystemSolverResults_t res;

    try {
      AbsTwoStepChronologicalPredictor4D<T>& two_step_predictor = 
	dynamic_cast<AbsTwoStepChronologicalPredictor4D<T>& >(predictor);

      // This is the unpreconditioned operator
      // Whenever we call its op() it will underneath call unprecLinOp()
      Lunprec<T,P,Q> A_unprec(A);
  
      SystemSolverResults_t individual_res;
      StopWatch swatch;
      swatch.reset();
      swatch.start();
      
      // we will solve for g5 D g5 D psi = chi
      // so psi = D^(-1)g5 D^(-1) g5 chi
      
      // Set up the source: g5 chi
      T g5chi = zero;
      g5chi[A->subset()] = Gamma(Nd*Nd-1)*chi ; 
      Double g5chi_norm = sqrt(norm2(g5chi, A->subset()));
      
      T tmpsol_tmp;  // This will hold our temporary solution
      T tmpferm=zero;

      // OK: predictY will predict the solution of M^\dagger Y = b
      // which is the same as Y = g_5 (M^{-1}) g_5 b
      // but actually we will really only want Y = (M^{-1}) g_5 b
      // so after prediction, we have to hit the prediction result
      // with gamma_5. TO AVOID ALLOCATING AN EXTRA FERMION
      // I am using tmpsol_tmp as the source. I am not using chi itself
      // in case for any reason it is dirty outside the target subset.
      tmpsol_tmp = zero;
      tmpsol_tmp[ A->subset() ] = chi;

      // predict with M^\dagger into tmpferm
      two_step_predictor.predictY(tmpferm, A_unprec, tmpsol_tmp);

      // hit tmpferm with g_5 to get the initial guess
      tmpsol_tmp = Gamma(Nd*Nd -1 )*tmpferm;
      Double init_norm = norm2(tmpsol_tmp, all);
      QDPIO::cout << "Initial guess norm_sq = " << init_norm << " norm = " << sqrt(init_norm) << std::endl;
 
      // Now do the solve
      individual_res = (*Dinv)(tmpsol_tmp,g5chi);


      T tmpsol = zero;
      tmpsol[ A->subset() ] = tmpsol_tmp; // Do this in case tmpsol_tmp is dirtied up outside its target
      res.n_count = individual_res.n_count;
    
      Double rel_resid = individual_res.resid / g5chi_norm;
    
      // If we've reached iteration-threshold, then destroy the subspace. The next solve will
      // refresh it.
      if ( individual_res.n_count >= invParam.RefreshThreshold ) {
	
	QDPIO::cout << "QOPMG_MDAGM_SOLVER: RefreshThreshold Iterations (" << invParam.RefreshThreshold <<")" 
		    << " reached in first LinopSolver call."
		    << " Erasing subspace: " << invParam.SubspaceId << std::endl;
	
	// Erase the subspace from Dinv and if saved from the 
	// NamedObjectMap
	(*Dinv).eraseSubspace();
	
	// Re solve if needed
	if ( toBool( rel_resid > invParam.Residual ) ) {
	  QDPIO::cout << "QOPMG_MDAGM_SOLVER: First solve failed with RelResid = " << rel_resid 
		      << " Re-trying with refreshed subspace " << std::endl;
	  
	  // tmpsol_tmp may already be a good guess... 
	  individual_res = (*Dinv)(tmpsol_tmp,g5chi);
	  tmpsol=zero;
	  tmpsol[ A->subset() ] = tmpsol_tmp; // Do this in case tmpsol_tmp is dirtied up outside of its target subset
	  res.n_count += individual_res.n_count;
	  
	  rel_resid = individual_res.resid / g5chi_norm;
	  if ( toBool( rel_resid > invParam.Residual) ) {
	    QDPIO::cout << "QOP_MG_MDAGM_SOLVER: Re-solving of first solve with refreshed space failed with RelResid = " << rel_resid 
			<<" Giving Up! " << std::endl;
	    QDP_abort(1);
	  }
	}
      }
      
    
      // At this point either 
      // -- the first solve worked.
      // -- OR the first solve didn't converge, got refreshed, resolved and worked
      // -- In either case, tmpsol, should be good and tmpsol_tmp holds the 
      // full solution on both subsets. Add it to the predictor.
      two_step_predictor.newYVector(tmpsol_tmp);

      
      T tmpsol2 = zero ;
      tmpsol2[A->subset()] =  Gamma(Nd*Nd-1)*tmpsol ;
      Double tmpsol2_norm = sqrt(norm2(tmpsol2,A->subset()));

      // This is less elaborate. predict X will call with A_unprec
      // no gamma_5 tricks needed
      tmpsol_tmp = zero;
      two_step_predictor.predictX(tmpsol_tmp, A_unprec, tmpsol2);

      init_norm = norm2(tmpsol_tmp, all);
      QDPIO::cout << "Initial guess norm_sq = " << init_norm << " norm = " << sqrt(init_norm) << std::endl;

      individual_res = (*Dinv)(tmpsol_tmp,tmpsol2);

      psi = zero;
      psi[ A->subset() ] = tmpsol_tmp; // Do this in case tmpsol tmp is dirty outside subset

      res.n_count += individual_res.n_count;
      rel_resid = individual_res.resid / tmpsol2_norm;
      
    // If we've reached max iter, then destroy the subspace. The next solve will
      if ( individual_res.n_count >= invParam.RefreshThreshold ) { 
	
	QDPIO::cout << "QOPMG_MDAGM_SOLVER: RefreshThreshold Iterations (" << invParam.RefreshThreshold <<")" 
		    << " reached in second LinopSolver call."
		    << " Erasing subspace: " << invParam.SubspaceId << std::endl;
	
	// Erase the subspace from Dinv and if saved from the 
	// NamedObjectMap
	(*Dinv).eraseSubspace();
	
	
	// Re solve if needed
	if ( toBool( rel_resid > invParam.Residual ) ) {
	  
	  QDPIO::cout << "QOPMG_MDAGM_SOLVER: Second solve failed with RelResid = " << rel_resid 
		      << " Re-trying with refreshed subspace " << std::endl;
	  
	  // tmpsol_tmp may already be a good initial guess
	  individual_res = (*Dinv)(tmpsol_tmp,tmpsol2); // This will internally refresh the subspace
	  psi=zero;
	  psi[ A->subset() ] = tmpsol_tmp; // Do this in case tmpsol_tmp is dirty outside subset
	  res.n_count += individual_res.n_count;
	  
	  rel_resid = individual_res.resid / tmpsol2_norm;
	  if ( toBool( rel_resid > invParam.Residual ) ) {
	    QDPIO::cout << "QOP_MG_MDAGM_SOLVER: Re-solving of second solve with refreshed space failed with RelResid = " << rel_resid 
			<<" Giving Up! " << std::endl;
	    
	    QDP_abort(1);
	  }
	}
      }


      
      // At this point either the second solve will have succeeded and was good
      // OR the second solve will have failed, got subspace refreshed, and resolved if needed.

      // So add tmpsol_tmp to the chrono
      two_step_predictor.newXVector(tmpsol_tmp);
      swatch.stop();
      
      double time = swatch.getTimeInSeconds();
      { 
	T r;
	r[A->subset()] = chi;
	T tmp;
	T tmp2 ;
	(*A)(tmp2, psi , PLUS );
	(*A)(tmp , tmp2, MINUS);
	
	r[A->subset()] -= tmp;
	res.resid = sqrt(norm2(r, A->subset()));
	rel_resid = res.resid / sqrt(norm2(chi, A->subset()));
      }
      
      QDPIO::cout << "QOPMG_MDAGM_SOLVER: "  
		  << " Mass: " << invParam.Mass 
		  << " iterations: " <<  res.n_count
		  << " time: " << time << " sec."
		  << " Rsd: " << res.resid
		  << " Relative Rsd: " << rel_resid
		  << std::endl;
      
      if ( toBool( rel_resid > invParam.RsdToleranceFactor*invParam.Residual ) ) { 
	if( invParam.TerminateOnFail ) { 
	  QDPIO::cout << "***** !!! ERROR !!! NONCONVERGENCE !!! Mass="<<invParam.Mass <<" Aborting !!! *******" << std::endl;
	  QDP_abort(1);
	}
	else { 
	  QDPIO::cout << "***** !!! WARNING!!! NONCONVERGENCE !!! Mass="<<invParam.Mass <<" !!! *******" << std::endl;
	}
      }
    }
    catch(std::bad_cast& bc) { 
      QDPIO::cout << "Couldnt cast predictor to two step predictor" << std::endl;
      QDP_abort(1);
    }
    
    END_CODE();
    return res;
  }
  
  //! Solve the linear system
  /*!
   * \param psi      solution ( Modify )
   * \param chi      source ( Read )
   * \param isign    solve with dagger or not
   * \return syssolver results
   */
  // T is the Lattice Fermion type
  SystemSolverResults_t
  MdagMSysSolverQOPMG::operator() (LatticeFermion& psi, const LatticeFermion& chi) const
  {
    /* Ignoring the predictor for now, since I am solving the UNPREC system anyway as a gateway to solve the PREC system. The prec
       solutions provided by the predictory are good for the PREC system, but may not be good enoug (missing half of) the solution
       of the the UNPREC system */
    Null4DChronoPredictor not_predicting;
    SystemSolverResults_t res  = (*this)(psi,chi, not_predicting);
    return res;
  }//! Solve the linear system
  
  
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
      return new MdagMSysSolverQOPMG(A, state, SysSolverQOPMGParams(xml_in, path));
    }
    
    /*MdagMSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path,
						  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state, 

						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new MdagMSysSolverQOPMG<LatticeFermion>(A, state, SysSolverQOPMGParams(xml_in, path));
    }*/
    
    /*//! Callback function for single precision
    MdagMSystemSolver<LatticeFermionF>*
      createFermF( XMLReader&                                          xml_in,
                   const std::string&                                  path,
                   Handle< FermState< LatticeFermionF, 
                                      multi1d<LatticeColorMatrixF>,
                                      multi1d<LatticeColorMatrixF> > > state,
                   Handle< LinearOperator<LatticeFermionF> >           A)
    {
      return new MdagMSysSolverQOPMG<LatticeFermionF>
                   (A, SysSolverQOPMGParams(xml_in, path));
    }*/

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true;
      if (! registered)
      {  
        success &= Chroma::TheMdagMFermSystemSolverFactory::Instance().registerObject(name, createFerm);
        //success &= Chroma::TheMdagMFFermSystemSolverFactory::Instance().registerObject(name, createFermF);
        registered = true;
      }
      return success;
    }
  }
}

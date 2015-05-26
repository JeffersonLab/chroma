// $Id: syssolver_linop_qdp_mg.cc, v1.0 2013-06-20 22:12 sdcohen $
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
#include "update/molecdyn/predictor/MG_predictor.h"



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
      ret_val.TerminateOnFail = false;
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
  MdagMSysSolverQOPMG::operator() (LatticeFermion& psi, const LatticeFermion& chi) const
  {
    START_CODE();
    typedef LatticeFermion T;
    SystemSolverResults_t res;
    SystemSolverResults_t individual_res;
    StopWatch swatch;
    swatch.reset();
    swatch.start();
    
    // we will solve for g5 D g5 D psi = chi
    // so psi = D^(-1)g5 D^(-1) g5 chi
    
    T g5chi = zero;
    g5chi[A->subset()] = Gamma(Nd*Nd-1)*chi ; 
    Double g5chi_norm = sqrt(norm2(g5chi, A->subset())); //
    
    T tmpsol_tmp = zero ;
    T tmpsol = zero;
    individual_res = (*Dinv)(tmpsol_tmp,g5chi);
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
    // -- In either case, tmpsol, should be good.
    
    T tmpsol2 = zero ;
    tmpsol_tmp = zero;

    tmpsol2[A->subset()] =  Gamma(Nd*Nd-1)*tmpsol ;
    
    Double tmpsol2_norm = sqrt(norm2(tmpsol2,A->subset()));
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
  MdagMSysSolverQOPMG::operator() (LatticeFermion& psi, const LatticeFermion& chi, AbsChronologicalPredictor4D<LatticeFermion>& predictor) const
  {
    /* Ignoring the predictor for now, since I am solving the UNPREC system anyway as a gateway to solve the PREC system. The prec
       solutions provided by the predictory are good for the PREC system, but may not be good enoug (missing half of) the solution
       of the the UNPREC system */
    SystemSolverResults_t res  = (*this)(psi,chi);
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

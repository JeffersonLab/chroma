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

namespace Chroma
{


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
      Dinv(new LinOpSysSolverQOPMG(A_,state_,invParam_))
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
    
    StopWatch swatch;
    swatch.reset();
    swatch.start();

    // we will solve for g5 D g5 D psi = chi
    // so psi = D^(-1)g5 D^(-1) g5 chi

    T g5chi = Gamma(Nd*Nd-1)*chi ; //
    T tmpsol  ;
    T tmpsol2 = zero ;
    (*Dinv)(tmpsol,g5chi);
    tmpsol2[A->subset()] =  Gamma(Nd*Nd-1)*tmpsol ;
    (*Dinv)(psi,tmpsol2);

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
    }

    QDPIO::cout << "MdagM_QOPMG_SOLVER: " //<< res.n_count << " iterations."
                << " Rsd = " << res.resid
                << " Relative Rsd = " << res.resid/sqrt(norm2(chi,A->subset()))
                << std::endl;
    QDPIO::cout << "MdagM_QOPMG_SOLVER_TIME: " << time << " secs" << std::endl;

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
    MG4DChronoPredictor<LatticeFermion>* MG_predictor = dynamic_cast<MG4DChronoPredictor<LatticeFermion>*>(&predictor);
    if( MG_predictor == 0x0 ) {
      QDPIO::cerr << "Unable to downcast AbsChronologicalPredictor4D to MG4DChronoPredictor." << std::endl;
      QDP_abort(1);
    }
    MG_predictor->getSubspace();
    (*this)(psi,chi);
    MG_predictor->resetSubspace(5);
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

/*
 * syssolver_mdagm_mrhs_twisted_proxy.cc
 *
 *  Created on: Mar 15, 2019
 *      Author: bjoo
 */
#include "chromabase.h"
#include "state.h"
#include "handle.h"
#include "actions/ferm/invert/syssolver_mdagm_mrhs_factory.h"
#include "actions/ferm/invert/syssolver_mrhs_twisted_params.h"
#include "actions/ferm/invert/syssolver_mrhs_twisted_proxy.h"

using namespace QDP;
using namespace Chroma;
namespace Chroma
{

  //! CG1 system solver namespace
  namespace MdagMSysSolverMRHSTwistedProxyEnv
  {
    //! Anonymous namespace
    namespace
    {
      //! Name to be used
      const std::string eoprec_name("MULTI_RHS_EOPREC_TWISTED_PROXY_INVERTER");
      const std::string seoprec_name("MULTI_RHS_SEOPREC_TWISTED_PROXY_INVERTER");
      //! Local registration flag
      bool registered = false;
    }


    //! Callback function
    MdagMMRHSSystemSolver<LatticeFermion>* createFermEoprec(XMLReader& xml_in,
						  const std::string& path,
						  Handle< FermAct4D< LatticeFermion,
						  	  	  	  	  	  	 multi1d<LatticeColorMatrix>,
												 multi1d<LatticeColorMatrix> > > S_ferm,
						  Handle< FermState<
						                     LatticeFermion,
						                     multi1d<LatticeColorMatrix>,
						                     multi1d<LatticeColorMatrix>
					 	  	  	  	  	   > > state)
    {
      using T = LatticeFermion;
      using P = multi1d<LatticeColorMatrix>;
      using Q = multi1d<LatticeColorMatrix>;

      return new EvenOddPrecMdagMMRHSSysSolverTwistedProxy<T,P,Q>(SysSolverMRHSTwistedParams(xml_in,path),S_ferm,state);

    }

    MdagMMRHSSystemSolver<LatticeFermion>* createFermSeoprec(XMLReader& xml_in,
						  const std::string& path,
						  Handle< FermAct4D< LatticeFermion,
						  	  	  	  	  	  	 multi1d<LatticeColorMatrix>,
												 multi1d<LatticeColorMatrix> > > S_ferm,
						  Handle< FermState<
						                     LatticeFermion,
						                     multi1d<LatticeColorMatrix>,
						                     multi1d<LatticeColorMatrix>
					 	  	  	  	  	   > > state)
    {
    	using T = LatticeFermion;
    	using P = multi1d<LatticeColorMatrix>;
    	using Q = multi1d<LatticeColorMatrix>;
    	return new SymEvenOddPrecLogDetMdagMMRHSSysSolverTwistedProxy<T,P,Q>(SysSolverMRHSTwistedParams(xml_in,path),S_ferm,state);

    }


    //! Register all the factories
    bool registerAll()
    {
    	bool success = true;
    	if (! registered)
    	{
    		success &= Chroma::TheMdagMFermMRHSSystemSolverFactory::Instance().registerObject(eoprec_name, createFermEoprec);
    		success &= Chroma::TheMdagMFermMRHSSystemSolverFactory::Instance().registerObject(seoprec_name, createFermSeoprec);
    		registered = true;
    	}
    	return success;
    }
  }
}





#include "chromabase.h"
#include "update/molecdyn/integrator/md_integrator_factory.h"
#include "update/molecdyn/integrator/lcm_exp_tdt.h"
#include "update/molecdyn/integrator/lcm_integrator_leaps.h"
#include "io/xmllog_io.h"
#include <string>


namespace Chroma 
{ 
  
  namespace LatColMatExpTdtIntegratorEnv 
  {
    namespace
    {
      AbsComponentIntegrator<multi1d<LatticeColorMatrix>, 
			     multi1d<LatticeColorMatrix> >* 
      createMDIntegrator(
			 XMLReader& xml, 
			 const std::string& path)
      {
	// Read the integrator params
	LatColMatExpTdtIntegratorParams p(xml, path);
    
	return new LatColMatExpTdtIntegrator(p);
      }
      
      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "LCM_EXP_T";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheMDComponentIntegratorFactory::Instance().registerObject(name, createMDIntegrator); 
	registered = true;
      }
      return success;
    }
  }
  
  
  LatColMatExpTdtIntegratorParams::LatColMatExpTdtIntegratorParams(XMLReader& xml_in, const std::string& path) 
  {
    XMLReader paramtop(xml_in, path);
    try {
      read(paramtop, "./n_steps", n_steps);
    }
    catch ( const std::string& e ) { 
      QDPIO::cout << "Error reading XML in LatColMatExpTdtIntegratorParams " << e << endl;
      QDP_abort(1);
    }
  }

  void read(XMLReader& xml, 
	    const std::string& path, 
	    LatColMatExpTdtIntegratorParams& p) {
    LatColMatExpTdtIntegratorParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, 
	     const std::string& path, 
	     const LatColMatExpTdtIntegratorParams& p) {
    push(xml, path);
    write(xml, "n_steps", p.n_steps);
    pop(xml);
  }

  void LatColMatExpTdtIntegrator::operator()( 
					     AbsFieldState<multi1d<LatticeColorMatrix>,
					     multi1d<LatticeColorMatrix> >& s, 
					     const Real& traj_length) const
  {
   
    START_CODE();
    Real dtau = traj_length / n_steps;
    for(int i=0; i < n_steps; i++) { 
      LCMMDIntegratorSteps::leapQ(dtau,s);
    }

    END_CODE();
    

  }


};

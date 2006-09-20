#include "chromabase.h"
#include "update/molecdyn/integrator/md_integrator_factory.h"
#include "update/molecdyn/integrator/lcm_minimum_norm2_integrator.h"

namespace Chroma 
{ 
  
  namespace LatColMatMinimumNorm2IntegratorEnv 
  {
    namespace
    {
      AbsMDIntegrator<multi1d<LatticeColorMatrix>, 
		      multi1d<LatticeColorMatrix> >* createMDIntegrator(
			XMLReader& xml, const std::string& path,  Handle< AbsHamiltonian<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >& H) 
      {
	// Read the integrator params
	LatColMatMinimumNorm2IntegratorParams p(xml, path);
    
	return new LatColMatMinimumNorm2Integrator(p, H);
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "LCM_MINIMUM_NORM_2ND_ORDER_INTEGRATOR";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheMDIntegratorFactory::Instance().registerObject(name, createMDIntegrator); 
	registered = true;
      }
      return success;
    }
  }


  LatColMatMinimumNorm2IntegratorParams::LatColMatMinimumNorm2IntegratorParams(XMLReader& xml_in, const std::string& path) 
  {
    XMLReader paramtop(xml_in, path);
    try {
      read(paramtop, "./n_steps", n_steps);
      if( paramtop.count("./lambda") == 1 ) { 
	read(paramtop, "./lambda", lambda);
      }
      else { 
	// Default lambda from de-forcrand paper
	lambda=Real(0.1931833275037836);
      }
      read(paramtop, "./tau0", tau0);
    }
    catch ( const std::string& e ) { 
      QDPIO::cout << "Error reading XML in LatColMatMinimumNorm2IntegratorParams " << e << endl;
      QDP_abort(1);
    }
  }

  void read(XMLReader& xml, 
	    const std::string& path, 
	    LatColMatMinimumNorm2IntegratorParams& p) {
    LatColMatMinimumNorm2IntegratorParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, 
	     const std::string& path, 
	     const LatColMatMinimumNorm2IntegratorParams& p) {
    push(xml, path);
    write(xml, "n_steps", p.n_steps);
    write(xml, "lambda", p.lambda);
    write(xml, "tau0", p.tau0);
  }


};

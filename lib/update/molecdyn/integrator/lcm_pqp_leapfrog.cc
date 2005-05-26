#include "chromabase.h"
#include "update/molecdyn/integrator/md_integrator_factory.h"
#include "update/molecdyn/integrator/lcm_pqp_leapfrog.h"
#include "io/xmllog_io.h"

#include <string>
#include "util/gauge/taproj.h"
#include "util/gauge/reunit.h"
#include "util/gauge/expmat.h"


namespace Chroma { 
  
  namespace LatColMatPQPLeapfrogIntegratorEnv {

    AbsMDIntegrator<multi1d<LatticeColorMatrix>, 
		    multi1d<LatticeColorMatrix> >* createMDIntegrator(
								       XMLReader& xml, const std::string& path,  Handle< AbsHamiltonian<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >& H) {
      // Read the integrator params
      LatColMatPQPLeapfrogIntegratorParams p(xml, path);
    
      return new LatColMatPQPLeapfrogIntegrator(p, H);
    }

    const std::string name = "LCM_PQP_LEAPFROG_INTEGRATOR";
    const bool registered = TheMDIntegratorFactory::Instance().registerObject(name, createMDIntegrator); 
  };

  LatColMatPQPLeapfrogIntegratorParams::LatColMatPQPLeapfrogIntegratorParams(XMLReader& xml_in, const std::string& path) 
  {
    XMLReader paramtop(xml_in, path);
    try {
      read(paramtop, "./n_steps", n_steps);
      read(paramtop, "./tau0", tau0);
    }
    catch ( const std::string& e ) { 
      QDPIO::cout << "Error reading XML in LatColMatPQPLeapfrogIntegratorParams " << e << endl;
      QDP_abort(1);
    }
  }

  void read(XMLReader& xml, 
	    const std::string& path, 
	    LatColMatPQPLeapfrogIntegratorParams& p) {
    LatColMatPQPLeapfrogIntegratorParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, 
	     const std::string& path, 
	     const LatColMatPQPLeapfrogIntegratorParams& p) {
    push(xml, path);
    write(xml, "n_steps", p.n_steps);
    write(xml, "tau0", p.tau0);
    pop(xml);
  }


};

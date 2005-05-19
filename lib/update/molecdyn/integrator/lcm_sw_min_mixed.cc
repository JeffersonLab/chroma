#include "chromabase.h"
#include "update/molecdyn/integrator/md_integrator_factory.h"
#include "update/molecdyn/integrator/lcm_sw_min_mixed.h"
#include "io/xmllog_io.h"

#include <string>
#include "util/gauge/taproj.h"
#include "util/gauge/reunit.h"
#include "util/gauge/expmat.h"


namespace Chroma { 
  
  namespace LatColMatSexWeinMixedIntegratorEnv {

    AbsMDIntegrator<multi1d<LatticeColorMatrix>, 
		    multi1d<LatticeColorMatrix> >* createMDIntegrator(
								       XMLReader& xml, const std::string& path,  Handle< AbsHamiltonian<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >& H) {
      // Read the integrator params
      LatColMatSexWeinMixedIntegratorParams p(xml, path);
    
      return new LatColMatSexWeinMixedIntegrator(p, H);
    }

    const std::string name = "LCM_TWO_SCALE_MINIMUM_NORM_INTEGRATOR";
    const bool registered = TheMDIntegratorFactory::Instance().registerObject(name, createMDIntegrator); 
  };

  LatColMatSexWeinMixedIntegratorParams::LatColMatSexWeinMixedIntegratorParams(XMLReader& xml_in, const std::string& path) 
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
      read(paramtop, "./n_short_steps", n_short_steps);
      read(paramtop, "./S_short_monomials", S_short_monomials);
      read(paramtop, "./S_long_monomials", S_long_monomials);
    
    }
    catch ( const std::string& e ) { 
      QDPIO::cout << "Error reading XML in LatColMatSexWeinMixedIntegratorParams " << e << endl;
      QDP_abort(1);
    }
  }

  void read(XMLReader& xml, 
	    const std::string& path, 
	    LatColMatSexWeinMixedIntegratorParams& p) {
    LatColMatSexWeinMixedIntegratorParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, 
	     const std::string& path, 
	     const LatColMatSexWeinMixedIntegratorParams& p) {
    push(xml, path);
    write(xml, "n_steps", p.n_steps);
    write(xml, "lambda", p.lambda);
    write(xml, "tau0", p.tau0);
    write(xml, "n_short_steps", p.n_short_steps);
    write(xml, "S_short_monomials", p.S_short_monomials);
    write(xml, "S_long_monomials", p.S_long_monomials);
    pop(xml);
  }


};

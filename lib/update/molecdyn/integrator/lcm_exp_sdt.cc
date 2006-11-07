#include "chromabase.h"
#include "update/molecdyn/integrator/md_integrator_factory.h"
#include "update/molecdyn/integrator/lcm_exp_sdt.h"
#include "update/molecdyn/integrator/lcm_integrator_leaps.h"
#include "io/xmllog_io.h"
#include <string>


namespace Chroma 
{ 
  
  namespace LatColMatExpSdtIntegratorEnv 
  {
    namespace
    {
      AbsComponentIntegrator<multi1d<LatticeColorMatrix>, 
			     multi1d<LatticeColorMatrix> >* 
      createMDIntegrator(
			 XMLReader& xml, 
			 const std::string& path,  
			 Handle< AbsHamiltonian<multi1d<LatticeColorMatrix>, 
			 multi1d<LatticeColorMatrix> > >& H) 
      {
	// Read the integrator params
	LatColMatExpSdtIntegratorParams p(xml, path);
    
	return new LatColMatExpSdtIntegrator(p, H);
      }
      
      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "LCM_EXP_S";

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
  
  
  LatColMatExpSdtIntegratorParams::LatColMatExpSdtIntegratorParams(XMLReader& xml_in, const std::string& path) 
  {
    XMLReader paramtop(xml_in, path);
    try {
      read(paramtop, "./n_steps", n_steps);
      read(paramtop, "./monomial_list", monomial_list);
    }
    catch ( const std::string& e ) { 
      QDPIO::cout << "Error reading XML in LatColMatExpSdtIntegratorParams " << e << endl;
      QDP_abort(1);
    }
  }

  void read(XMLReader& xml, 
	    const std::string& path, 
	    LatColMatExpSdtIntegratorParams& p) {
    LatColMatExpSdtIntegratorParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, 
	     const std::string& path, 
	     const LatColMatExpSdtIntegratorParams& p) {
    push(xml, path);
    write(xml, "n_steps", p.n_steps);
    write(xml, "monomial_list", p.monomial_list);
    pop(xml);
  }

  void LatColMatExpSdtIntegrator::operator()( 
					     AbsFieldState<multi1d<LatticeColorMatrix>,
					     multi1d<LatticeColorMatrix> >& s, 
					     const Real& traj_length) const
  {
   
    START_CODE();
    Real dtau = traj_length / n_steps;
    for(int i=0; i < n_steps; i++) { 
      LCMMDIntegratorSteps::leapP( monomial_list,
				   dtau,
				   (*H_MD), 
				   s  );
    }

    END_CODE();
    

  }


};

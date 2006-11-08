#include "chromabase.h"
#include "update/molecdyn/integrator/md_integrator_factory.h"
#include "update/molecdyn/integrator/lcm_sts_leapfrog_component.h"
#include "update/molecdyn/integrator/lcm_exp_sdt.h"
#include "update/molecdyn/integrator/lcm_exp_tdt.h"
#include "io/xmllog_io.h"
#include <string>


namespace Chroma 
{ 
  
  namespace LatColMatSTSLeapfrogComponentIntegratorEnv 
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
	LatColMatSTSLeapfrogComponentIntegratorParams p(xml, path);
    
	return new LatColMatSTSLeapfrogComponentIntegrator(p, H);
      }
      
      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "LCM_STS_LEAPFROG_COMPONENT";

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
  
  
  LatColMatSTSLeapfrogComponentIntegratorParams::LatColMatSTSLeapfrogComponentIntegratorParams(XMLReader& xml_in, const std::string& path) 
  {
    XMLReader paramtop(xml_in, path);
    try {
      read(paramtop, "./n_steps", n_steps);
      read(paramtop, "./monomial_list", monomial_list);
    }
    catch ( const std::string& e ) { 
      QDPIO::cout << "Error reading XML in LatColMatSTSLeapfrogComponentIntegratorParams " << e << endl;
      QDP_abort(1);
    }
  }

  void read(XMLReader& xml, 
	    const std::string& path, 
	    LatColMatSTSLeapfrogComponentIntegratorParams& p) {
    LatColMatSTSLeapfrogComponentIntegratorParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, 
	     const std::string& path, 
	     const LatColMatSTSLeapfrogComponentIntegratorParams& p) {
    push(xml, path);
    write(xml, "n_steps", p.n_steps);
    write(xml, "monomial_list", p.monomial_list);
    pop(xml);
  }

  void LatColMatSTSLeapfrogComponentIntegrator::operator()( 
					     AbsFieldState<multi1d<LatticeColorMatrix>,
					     multi1d<LatticeColorMatrix> >& s, 
					     const Real& traj_length) const
  {
   
    START_CODE();
    LatColMatExpSdtIntegrator expSdt(1,
				     monomial_list,
				     H_MD);    // Single Step

    LatColMatExpTdtIntegrator expTdt(1);       // Single Step
				     


    Real dtau = traj_length / n_steps;
    Real dtauby2 = dtau/2;

    // Its sts so:
    expSdt(s, dtauby2);  // First half step
    for(int i=0; i < n_steps-1; i++) {  // N-1 full steps
      expTdt(s, dtau);
      expSdt(s, dtau);
    }
    expTdt(s, dtau);     // Last Full Step
    expSdt(s, dtauby2);  // Last Half Step

    END_CODE();
    

  }


};

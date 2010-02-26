#include "chromabase.h"
#include "update/molecdyn/integrator/md_integrator_factory.h"
#include "update/molecdyn/integrator/lcm_exp_sdt_sstdt3.h"
#include "update/molecdyn/integrator/lcm_integrator_leaps.h"
#include "io/xmllog_io.h"
#include <string>


namespace Chroma 
{ 
  
  namespace LatColMatExpSdtMinusSSTdt3IntegratorEnv 
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
	LatColMatExpSdtMinusSSTdt3IntegratorParams p(xml, path);
    
	return new LatColMatExpSdtMinusSSTdt3Integrator(p);
      }
      
      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "LCM_EXP_S_SST";

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
  
  
  LatColMatExpSdtMinusSSTdt3IntegratorParams::LatColMatExpSdtMinusSSTdt3IntegratorParams(XMLReader& xml_in, const std::string& path) 
  {
    XMLReader paramtop(xml_in, path);
    try {
      read(paramtop, "./n_steps", n_steps);
      read(paramtop, "./monomial_list", monomial_list);
    }
    catch ( const std::string& e ) { 
      QDPIO::cout << "Error reading XML in LatColMatExpSdtMinusSSTdt3IntegratorParams " << e << endl;
      QDP_abort(1);
    }
  }

  void read(XMLReader& xml, 
	    const std::string& path, 
	    LatColMatExpSdtMinusSSTdt3IntegratorParams& p) {
    LatColMatExpSdtMinusSSTdt3IntegratorParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, 
	     const std::string& path, 
	     const LatColMatExpSdtMinusSSTdt3IntegratorParams& p) {
    push(xml, path);
    write(xml, "n_steps", p.n_steps);
    write(xml, "monomial_list", p.monomial_list);
    pop(xml);
  }

  void LatColMatExpSdtMinusSSTdt3Integrator::operator()( 
					     AbsFieldState<multi1d<LatticeColorMatrix>,
					     multi1d<LatticeColorMatrix> >& s, 
					     const Real& t1) const
  {
   
    START_CODE();
    
      LCMMDIntegratorSteps::leapFG( monomials,
				    t1,
				    s  );


    END_CODE();
    

  }


};

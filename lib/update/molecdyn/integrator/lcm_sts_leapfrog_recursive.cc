#include "chromabase.h"
#include "update/molecdyn/integrator/md_integrator_factory.h"
#include "update/molecdyn/integrator/lcm_sts_leapfrog_recursive.h"
#include "update/molecdyn/integrator/lcm_exp_sdt.h"
#include "io/xmllog_io.h"

#include <string>
using namespace std;

namespace Chroma 
{ 
  
  namespace LatColMatSTSLeapfrogRecursiveIntegratorEnv 
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
	LatColMatSTSLeapfrogRecursiveIntegratorParams p(xml, path);
    
	return new LatColMatSTSLeapfrogRecursiveIntegrator(p, H);
      }
      
      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "LCM_STS_LEAPFROG_RECURSIVE";

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
  
  
  LatColMatSTSLeapfrogRecursiveIntegratorParams::LatColMatSTSLeapfrogRecursiveIntegratorParams(XMLReader& xml_in, const std::string& path) 
  {
    XMLReader paramtop(xml_in, path);
    try {
      read(paramtop, "./n_steps", n_steps);
      read(paramtop, "./monomial_list", monomial_list);
      XMLReader subint_reader(paramtop, "./SubIntegrator");
      
      std::ostringstream subintegrator_os;
      subint_reader.print(subintegrator_os);
      subintegrator_xml = subintegrator_os.str();
      QDPIO::cout << "Subintegrator XML is: " << endl;
      QDPIO::cout << subintegrator_xml << endl;
    }
    catch ( const std::string& e ) { 
      QDPIO::cout << "Error reading XML in LatColMatSTSLeapfrogRecursiveIntegratorParams " << e << endl;
      QDP_abort(1);
    }
  }
  
  void read(XMLReader& xml, 
	    const std::string& path, 
	    LatColMatSTSLeapfrogRecursiveIntegratorParams& p) {
    LatColMatSTSLeapfrogRecursiveIntegratorParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, 
	     const std::string& path, 
	     const LatColMatSTSLeapfrogRecursiveIntegratorParams& p) {
    push(xml, path);
    write(xml, "n_steps", p.n_steps);
    write(xml, "monomial_list", p.monomial_list);
    xml << p.subintegrator_xml;
    pop(xml);
  }

  
  AbsComponentIntegrator< multi1d<LatticeColorMatrix>,
			  multi1d<LatticeColorMatrix> >* LatColMatSTSLeapfrogRecursiveIntegrator::createSubIntegrator(const std::string& subintegrator_xml, 
														      Handle< AbsHamiltonian< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >& H_) {
    
    std::istringstream is( subintegrator_xml );
    XMLReader top(is);
    
    std::string subint_name;
    try { 
      read(top, "/SubIntegrator/Name", subint_name);
    }
    catch( const std::string ) {
      QDPIO::cerr << "Failed to extract name of subintegrator in LatColMatSTSLeapfrogRecursiveIntegrator" << endl;
      QDP_abort(1);
    }
    std::string root="/SubIntegrator";
    
    AbsComponentIntegrator< multi1d<LatticeColorMatrix>,
      multi1d<LatticeColorMatrix> >* ret_val=
     TheMDComponentIntegratorFactory::Instance().createObject(subint_name, top, root, H_);
    return ret_val;

  }


  void LatColMatSTSLeapfrogRecursiveIntegrator::operator()( 
					     AbsFieldState<multi1d<LatticeColorMatrix>,
					     multi1d<LatticeColorMatrix> >& s, 
					     const Real& traj_length) const
  {
   
    START_CODE();
    LatColMatExpSdtIntegrator expSdt(1,
				     monomial_list,
				     H_MD);    // Single Step

    const AbsComponentIntegrator< multi1d<LatticeColorMatrix>,
      multi1d<LatticeColorMatrix> >& subIntegrator = getSubIntegrator();

				     


    Real dtau = traj_length / n_steps;
    Real dtauby2 = dtau/2;

    // Its sts so:
    expSdt(s, dtauby2);  // First half step
    for(int i=0; i < n_steps-1; i++) {  // N-1 full steps
      subIntegrator(s, dtau); 
      expSdt(s, dtau);
    }
    subIntegrator(s, dtau);     // Last Full Step
    expSdt(s, dtauby2);  // Last Half Step

    END_CODE();
    

  }


};

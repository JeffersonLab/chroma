#include "chromabase.h"
#include "update/molecdyn/integrator/md_integrator_factory.h"
#include "update/molecdyn/integrator/lcm_creutz_gocksch_4_recursive.h"
#include "update/molecdyn/integrator/lcm_exp_sdt.h"
#include "io/xmllog_io.h"

#include <string>
using namespace std;

namespace Chroma 
{ 
  
  namespace LatColMatCreutzGocksch4RecursiveIntegratorEnv 
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
	LatColMatCreutzGocksch4RecursiveIntegratorParams p(xml, path);
    
	return new LatColMatCreutzGocksch4RecursiveIntegrator(p);
      }
      
      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "LCM_CREUTZ_GOCKSCH_4";

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
  
  
  LatColMatCreutzGocksch4RecursiveIntegratorParams::LatColMatCreutzGocksch4RecursiveIntegratorParams(XMLReader& xml_in, const std::string& path) 
  {
    XMLReader paramtop(xml_in, path);
    // Default values for the tuning constants

    try {
      read(paramtop, "./n_steps", n_steps);
      read(paramtop, "./monomial_ids", monomial_ids);

      if( paramtop.count("./SubIntegrator") == 0 ) {
	// BASE CASE: User does not supply sub-integrator 
	//
	// Sneaky way - create an XML document for EXP_T
	XMLBufferWriter subintegrator_writer;
	int one_sub_step=1;

	push(subintegrator_writer, "SubIntegrator");
	write(subintegrator_writer, "Name", "LCM_EXP_T");
	write(subintegrator_writer, "n_steps", one_sub_step);

	pop(subintegrator_writer);

	subintegrator_xml = subintegrator_writer.str();

      }
      else {
	// RECURSIVE CASE: User Does Supply Sub Integrator
	//
	// Read it
	XMLReader subint_reader(paramtop, "./SubIntegrator");
      
	std::ostringstream subintegrator_os;
      
	subint_reader.print(subintegrator_os);
	subintegrator_xml = subintegrator_os.str();
	QDPIO::cout << "Subintegrator XML is: " << endl;
	QDPIO::cout << subintegrator_xml << endl;
      }
    }
    catch ( const std::string& e ) { 
      QDPIO::cout << "Error reading XML in LatColMatCreutzGocksch4RecursiveIntegratorParams " << e << endl;
      QDP_abort(1);
    }
  }
  
  void read(XMLReader& xml, 
	    const std::string& path, 
	    LatColMatCreutzGocksch4RecursiveIntegratorParams& p) {
    LatColMatCreutzGocksch4RecursiveIntegratorParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, 
	     const std::string& path, 
	     const LatColMatCreutzGocksch4RecursiveIntegratorParams& p) {
    push(xml, path);
    write(xml, "n_steps", p.n_steps);
    write(xml, "monomial_ids", p.monomial_ids);

    xml << p.subintegrator_xml;

    pop(xml);
  }

  
  void LatColMatCreutzGocksch4RecursiveIntegrator::operator()( 
					     AbsFieldState<multi1d<LatticeColorMatrix>,
					     multi1d<LatticeColorMatrix> >& state, 
					     const Real& traj_length) const
  {
   
    START_CODE();
    LatColMatExpSdtIntegrator expSdt(1,
				     monomials); // Single Step

    const AbsComponentIntegrator< multi1d<LatticeColorMatrix>,
      multi1d<LatticeColorMatrix> >& subIntegrator = getSubIntegrator();

				     
    Real one_by_three = Real(1)/Real(3);
    Real cube_root_two = pow(Real(2), one_by_three);

    Real dtau = traj_length / params.n_steps;


    Real s= dtau/(Real(2) - cube_root_two);
    Real t = -cube_root_two*dtau/ (Real(2) - cube_root_two);

    Real half_s_plus_t = (s + t)/Real(2);
    Real half_s = s/Real(2);

    expSdt(state, half_s);
    
    subIntegrator(state, s);

    expSdt(state, half_s_plus_t);

    subIntegrator(state, t);
      
    expSdt(state, half_s_plus_t);

    subIntegrator(state, s);

    for(int step=0; step < params.n_steps-1; ++step) {
      // One for the end of the previous and one for the start of current
      expSdt(state, s);

      subIntegrator(state, s);

      expSdt(state, half_s_plus_t);

      subIntegrator(state, t);

      expSdt(state, half_s_plus_t);

      subIntegrator(state, s);
     

    }
    
    expSdt(state, half_s);

    END_CODE();

  }


};

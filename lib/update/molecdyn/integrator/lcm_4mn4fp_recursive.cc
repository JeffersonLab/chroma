#include "chromabase.h"
#include "update/molecdyn/integrator/md_integrator_factory.h"
#include "update/molecdyn/integrator/lcm_4mn4fp_recursive.h"
#include "update/molecdyn/integrator/lcm_exp_sdt.h"
#include "io/xmllog_io.h"

#include <string>
using namespace std;

namespace Chroma 
{ 
  
  namespace LatColMat4MN4FPRecursiveIntegratorEnv 
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
	LatColMat4MN4FPRecursiveIntegratorParams p(xml, path);
    
	return new LatColMat4MN4FPRecursiveIntegrator(p);
      }
      
      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "LCM_4MN4FP";

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
  
  
  LatColMat4MN4FPRecursiveIntegratorParams::LatColMat4MN4FPRecursiveIntegratorParams(XMLReader& xml_in, const std::string& path) 
  {
    XMLReader paramtop(xml_in, path);
    // Default values for the tuning constants
    rho =     Real(0.1786178958448091);
    theta =   Real(-0.06626458266981843);
    lambda =  Real(0.7123418310626056);

    try {
      read(paramtop, "./n_steps", n_steps);
      read(paramtop, "./monomial_ids", monomial_ids);

      // Override default tuning values if desired
      if ( paramtop.count("./theta") > 0 ) { 
	read(paramtop, "./theta", theta);
      }
      if ( paramtop.count("./rho") > 0 ) { 
	read(paramtop, "./rho", rho);
      }
      if ( paramtop.count("./lambda") > 0 ) {
	read(paramtop, "./lambda", lambda);
      }

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
      QDPIO::cout << "Error reading XML in LatColMat4MN4FPRecursiveIntegratorParams " << e << endl;
      QDP_abort(1);
    }
  }
  
  void read(XMLReader& xml, 
	    const std::string& path, 
	    LatColMat4MN4FPRecursiveIntegratorParams& p) {
    LatColMat4MN4FPRecursiveIntegratorParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, 
	     const std::string& path, 
	     const LatColMat4MN4FPRecursiveIntegratorParams& p) {
    push(xml, path);
    write(xml, "n_steps", p.n_steps);
    write(xml, "monomial_ids", p.monomial_ids);
    write(xml, "theta", p.theta);
    write(xml, "rho", p.rho);
    write(xml, "lambda", p.lambda);

    xml << p.subintegrator_xml;

    pop(xml);
  }

  
  void LatColMat4MN4FPRecursiveIntegrator::operator()( 
					     AbsFieldState<multi1d<LatticeColorMatrix>,
					     multi1d<LatticeColorMatrix> >& s, 
					     const Real& traj_length) const
  {
   
    START_CODE();
    LatColMatExpSdtIntegrator expSdt(1,
				     monomials); // Single Step

    const AbsComponentIntegrator< multi1d<LatticeColorMatrix>,
      multi1d<LatticeColorMatrix> >& subIntegrator = getSubIntegrator();

				     

    
    Real dtau = traj_length / params.n_steps;
    
    // Various timestep combinations
    Real theta_dtau = params.theta * dtau;
    Real rho_dtau = params.rho * dtau;
    Real lambda_dtau = params.lambda * dtau;


    Real one_minus_two_lambda_dtau_by2 = (Real(1)-Real(2)*params.lambda)*dtau/Real(2);

    Real one_minus_two_theta_plus_rho_dtau = (Real(1)-Real(2)*(params.theta+params.rho))*dtau;

    Real two_rho_dtau= Real(2)*rho_dtau;

    // From deForcrand-Takaishi paper (eq 24)
    // But unrolled
    // Map exp( V ) -> expS
    // Map exp( T ) -> subIntegrator() (eg at lowest level of recursion)

    subIntegrator(s, rho_dtau);
    expSdt(s, lambda_dtau);
    subIntegrator(s, theta_dtau);
    expSdt(s, one_minus_two_lambda_dtau_by2 );
    subIntegrator(s, one_minus_two_theta_plus_rho_dtau);
    expSdt(s, one_minus_two_lambda_dtau_by2 );
    subIntegrator(s, theta_dtau);
    expSdt(s, lambda_dtau);

    for(int steps = 0; steps < params.n_steps-1; ++steps) {
      subIntegrator(s, two_rho_dtau);
      expSdt(s, lambda_dtau);
      subIntegrator(s, theta_dtau);
      expSdt(s, one_minus_two_lambda_dtau_by2 );
      subIntegrator(s, one_minus_two_theta_plus_rho_dtau);
      expSdt(s, one_minus_two_lambda_dtau_by2 );
      subIntegrator(s, theta_dtau);
      expSdt(s, lambda_dtau);

    }

    subIntegrator(s, rho_dtau);

    END_CODE();

  }


};

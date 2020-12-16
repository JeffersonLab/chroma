#include "chromabase.h"
#include "update/molecdyn/integrator/md_integrator_factory.h"
#include "update/molecdyn/integrator/lcm_4mn5fv_recursive.h"
#include "update/molecdyn/integrator/lcm_exp_sdt.h"
#include "io/xmllog_io.h"

#include <string>

namespace Chroma 
{ 
  
  namespace LatColMat4MN5FVRecursiveIntegratorEnv 
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
	LatColMat4MN5FVRecursiveIntegratorParams p(xml, path);
    
	return new LatColMat4MN5FVRecursiveIntegrator(p);
      }
      
      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "LCM_4MN5FV";

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
  
  
  LatColMat4MN5FVRecursiveIntegratorParams::LatColMat4MN5FVRecursiveIntegratorParams(XMLReader& xml_in, const std::string& path) 
  {
    XMLReader paramtop(xml_in, path);
    // Default values for the tuning constants
    theta =   Real(0.08398315262876693);
    rho =     Real(0.2539785108410595);
    lambda =  Real(0.6822365335719091);
    mu =     Real(-0.03230286765269967);
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
      if ( paramtop.count("./mu") > 0 ) { 
	read(paramtop, "./mu", mu);
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
	QDPIO::cout << "Subintegrator XML is: " << std::endl;
	QDPIO::cout << subintegrator_xml << std::endl;
      }
    }
    catch ( const std::string& e ) { 
      QDPIO::cout << "Error reading XML in LatColMat4MN5FVRecursiveIntegratorParams " << e << std::endl;
      QDP_abort(1);
    }
  }
  
  void read(XMLReader& xml, 
	    const std::string& path, 
	    LatColMat4MN5FVRecursiveIntegratorParams& p) {
    LatColMat4MN5FVRecursiveIntegratorParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, 
	     const std::string& path, 
	     const LatColMat4MN5FVRecursiveIntegratorParams& p) {
    push(xml, path);
    write(xml, "n_steps", p.n_steps);
    write(xml, "monomial_ids", p.monomial_ids);
    write(xml, "theta", p.theta);
    write(xml, "rho", p.rho);
    write(xml, "lambda", p.lambda);
    write(xml, "mu", p.mu);

    xml << p.subintegrator_xml;

    pop(xml);
  }

  
  void LatColMat4MN5FVRecursiveIntegrator::operator()( 
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
    Real two_theta_dtau = Real(2)*theta_dtau;
    Real rho_dtau = params.rho * dtau;
    Real lambda_dtau = params.lambda * dtau;
    Real mu_dtau = params.mu*dtau;

    Real one_minus_two_lambda_plus_theta_dtau_by2 = (Real(1)-Real(2)*(params.lambda+params.theta))*dtau/Real(2);

    Real one_minus_two_mu_plus_rho_dtau = (Real(1)-Real(2)*(params.mu+params.rho))*dtau;


    // From deForcrand-Takaishi paper (eq 24)
    // But unrolled
    // Map exp( V ) -> expS
    // Map exp( T ) -> subIntegrator() (eg at lowest level of recursion)
    expSdt(s, theta_dtau);
    
    subIntegrator(s, rho_dtau);

    expSdt(s, lambda_dtau);
    
    subIntegrator(s, mu_dtau);
    
    expSdt(s, one_minus_two_lambda_plus_theta_dtau_by2);
    
    subIntegrator(s, one_minus_two_mu_plus_rho_dtau);
    
    expSdt(s, one_minus_two_lambda_plus_theta_dtau_by2);
    
    subIntegrator(s, mu_dtau);
    
    expSdt(s, lambda_dtau);
    
    subIntegrator(s, rho_dtau);
     
    for(int step=0; step < params.n_steps-1; ++step) {
      // One for the end of the previous and one for the start of current
      expSdt(s, two_theta_dtau);
      
      subIntegrator(s, rho_dtau);

      expSdt(s, lambda_dtau);

      subIntegrator(s, mu_dtau);

      expSdt(s, one_minus_two_lambda_plus_theta_dtau_by2);

      subIntegrator(s, one_minus_two_mu_plus_rho_dtau);

      expSdt(s, one_minus_two_lambda_plus_theta_dtau_by2);

      subIntegrator(s, mu_dtau);

      expSdt(s, lambda_dtau);

      subIntegrator(s, rho_dtau);
    }

    // Finish off the last one
    expSdt(s, theta_dtau);

    END_CODE();

  }


}

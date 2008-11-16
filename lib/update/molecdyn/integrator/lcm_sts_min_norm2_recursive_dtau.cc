#include "chromabase.h"
#include "update/molecdyn/integrator/md_integrator_factory.h"
#include "update/molecdyn/integrator/lcm_sts_min_norm2_recursive_dtau.h"
#include "update/molecdyn/integrator/lcm_exp_sdt.h"
#include "io/xmllog_io.h"

#include <string>
using namespace std;

namespace Chroma 
{ 
  
  namespace LatColMatSTSMinNorm2DTauRecursiveIntegratorEnv 
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
	LatColMatSTSMinNorm2DTauRecursiveIntegratorParams p(xml, path);
    
	return new LatColMatSTSMinNorm2DTauRecursiveIntegrator(p);
      }
      
      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "LCM_STS_MIN_NORM_2_DTAU";

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
  
  
  LatColMatSTSMinNorm2DTauRecursiveIntegratorParams::LatColMatSTSMinNorm2DTauRecursiveIntegratorParams(XMLReader& xml_in, const std::string& path) 
  {
    XMLReader paramtop(xml_in, path);
    try {
      read(paramtop, "./delta_tau_max", delta_tau_max);
      read(paramtop, "./monomial_ids", monomial_ids);
      if( paramtop.count("./lambda") == 1 ) { 
	read(paramtop, "./lambda", lambda );
      }
      else { 
	lambda = 0.1931833275037836;

	QDPIO::cout << "Warning no lambda param found for minimum norm integrator" << endl;
	QDPIO::cout << "Using default value lambda_c="<< lambda << endl;
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
      QDPIO::cout << "Error reading XML in LatColMatSTSMinNorm2DTauRecursiveIntegratorParams " << e << endl;
      QDP_abort(1);
    }
  }
  
  void read(XMLReader& xml, 
	    const std::string& path, 
	    LatColMatSTSMinNorm2DTauRecursiveIntegratorParams& p) {
    LatColMatSTSMinNorm2DTauRecursiveIntegratorParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, 
	     const std::string& path, 
	     const LatColMatSTSMinNorm2DTauRecursiveIntegratorParams& p) {
    push(xml, path);
    write(xml, "delta_tau_max", p.delta_tau_max);
    write(xml, "monomial_ids", p.monomial_ids);
    write(xml, "lambda", p.lambda);
    xml << p.subintegrator_xml;
    pop(xml);
  }

  

  void LatColMatSTSMinNorm2DTauRecursiveIntegrator::operator()( 
					     AbsFieldState<multi1d<LatticeColorMatrix>,
					     multi1d<LatticeColorMatrix> >& s, 
					     const Real& traj_length) const
  {
   
    START_CODE();
    LatColMatExpSdtIntegrator expSdt(1,
				     monomials);


    const AbsComponentIntegrator< multi1d<LatticeColorMatrix>,
      multi1d<LatticeColorMatrix> >& subIntegrator = getSubIntegrator();

    Real dtau;
    int n_steps;

    // Set up dtau and n_steps to do the minimum amount of 
    // work without exceeding a step-size of delta_tau_max

 
    if( toBool(traj_length > delta_tau_max) ) { 
      // delta_tau_max is SHORTER than traj_length => multiple steps.
      // pick  n = ceiling of traj_length/delta_tau_max.
      // Then find dtau = traj_length / n
      n_steps =toInt(ceil( traj_length/delta_tau_max ));
      dtau = traj_length/Real(n_steps);
    }
    else { 
      // delta_tau_max is BIGGER or EQUAL to traj_length.
      // In either case the best we can do is one step of
      // dtau = traj_length
      n_steps = 1;
      dtau = traj_length;
    }

    //    Real dtau = traj_length / Real(n_steps);
    Real lambda_dt = dtau*lambda;
    Real dtauby2 = dtau / Real(2);
    Real one_minus_2lambda_dt = (Real(1)-Real(2)*lambda)*dtau;
    Real two_lambda_dt = lambda_dt*Real(2);

    // Its sts so:
    expSdt(s, lambda_dt); 
    for(int i=0; i < n_steps-1; i++) {  // N-1 full steps
      // Roll the exp(lambda_dt T) here and start
      // Next iter into one
      subIntegrator(s, dtauby2);
      expSdt(s, one_minus_2lambda_dt);
      subIntegrator(s, dtauby2);
      expSdt(s, two_lambda_dt); 
    }
    // Last step, can't roll the first and last exp(lambda_dt T) 
    // together.
    subIntegrator(s, dtauby2);
    expSdt(s, one_minus_2lambda_dt);
    subIntegrator(s, dtauby2);
    expSdt(s, lambda_dt);


    END_CODE();
    

  }


};

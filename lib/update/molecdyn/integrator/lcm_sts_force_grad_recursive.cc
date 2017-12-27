#include "chromabase.h"
#include "update/molecdyn/integrator/md_integrator_factory.h"
#include "update/molecdyn/integrator/lcm_sts_force_grad_recursive.h"
#include "update/molecdyn/integrator/lcm_integrator_leaps.h"
#include "update/molecdyn/integrator/lcm_exp_sdt.h"
#include "io/xmllog_io.h"

#include <string>

namespace Chroma 
{ 

  namespace LatColMatSTSForceGradRecursiveIntegratorEnv 
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
            LatColMatSTSForceGradRecursiveIntegratorParams p(xml, path);

            return new LatColMatSTSForceGradRecursiveIntegrator(p);
          }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "LCM_STS_FORCE_GRAD";

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


  LatColMatSTSForceGradRecursiveIntegratorParams::LatColMatSTSForceGradRecursiveIntegratorParams(XMLReader& xml_in, const std::string& path) 
  {
    XMLReader paramtop(xml_in, path);
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
        QDPIO::cout << "Subintegrator XML is: " << std::endl;
        QDPIO::cout << subintegrator_xml << std::endl;
      }
    }
    catch ( const std::string& e ) { 
      QDPIO::cout << "Error reading XML in LatColMatSTSForceGradRecursiveIntegratorParams " << e << std::endl;
      QDP_abort(1);
    }
  }

  void read(XMLReader& xml, 
      const std::string& path, 
      LatColMatSTSForceGradRecursiveIntegratorParams& p) {
    LatColMatSTSForceGradRecursiveIntegratorParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, 
      const std::string& path, 
      const LatColMatSTSForceGradRecursiveIntegratorParams& p) {
    push(xml, path);
    write(xml, "n_steps", p.n_steps);
    write(xml, "monomial_ids", p.monomial_ids);
    xml << p.subintegrator_xml;
    pop(xml);
  }


  // Force-gradient evolution
  // The code will update momentum P in state s
  void LatColMatSTSForceGradRecursiveIntegrator::fg_update( 
      AbsFieldState<multi1d<LatticeColorMatrix>,
      multi1d<LatticeColorMatrix> >& s, 
      const Real& traj_length1, const Real& traj_length2) const
  {
    LatColMatExpSdtIntegrator expSdt(1, monomials);

    const AbsComponentIntegrator< multi1d<LatticeColorMatrix>,
          multi1d<LatticeColorMatrix> >& subIntegrator = getSubIntegrator();

    // Temporary state information
    multi1d<LatticeColorMatrix> Q_tmp(Nd);
    multi1d<LatticeColorMatrix> P_tmp(Nd);

    // Save state
    P_tmp = s.getP();
    Q_tmp = s.getQ();

    // Now we will use the state s as a temporary state that will be used 
    // for force-gradient update of eqs. (2.4--5) in arXiv:111.5059

    // Calculate force term and initialize momentum P in s with the force
    // This will evaluate p' <- p + dt1 F 
    // with p=0 it becomes p = + dt1 F 
    s.getP() = zero;
    expSdt(s, traj_length1);

    // Update gauge field with the P calculated above;
    // Not doing this with the sub integrator, as that would 
    // likely mess up the momenta.
    // Instead I go straight for the leapQ.
    LCMMDIntegratorSteps::leapQ(1.0,s);

    // restore momentum P
    s.getP() = P_tmp;

    // Update momentum in s using the U'
    expSdt(s, traj_length2);

    // Restore gauge field back to U
    s.getQ() = Q_tmp;
  }



  void LatColMatSTSForceGradRecursiveIntegrator::operator()( 
      AbsFieldState<multi1d<LatticeColorMatrix>,
      multi1d<LatticeColorMatrix> >& s, 
      const Real& traj_length) const
  {

    START_CODE();
    LatColMatExpSdtIntegrator expSdt(1, monomials);


    const AbsComponentIntegrator< multi1d<LatticeColorMatrix>,
          multi1d<LatticeColorMatrix> >& subIntegrator = getSubIntegrator();

    Real lambda = Real(1)/Real(6);
    Real xi = Real(1)/Real(72);

    Real dtau = traj_length / Real(n_steps);
    Real lambda_dt = dtau*lambda;
    Real dtauby2 = dtau / Real(2);
    Real one_minus_2lambda_dt = (Real(1)-Real(2)*lambda)*dtau;
    Real two_lambda_dt = lambda_dt*Real(2);
    Real xi_dtdt = Real(2)*dtau*dtau*dtau*xi / one_minus_2lambda_dt;

    expSdt(s, lambda_dt); 
    for(int i=0; i < n_steps-1; i++) {  // N-1 full steps
      subIntegrator(s, dtauby2);
      fg_update(s, xi_dtdt, one_minus_2lambda_dt);
      subIntegrator(s, dtauby2);
      expSdt(s, two_lambda_dt); 
    }
    subIntegrator(s, dtauby2);
    fg_update(s, xi_dtdt, one_minus_2lambda_dt);
    subIntegrator(s, dtauby2);
    expSdt(s, lambda_dt);


    END_CODE();

  }


};

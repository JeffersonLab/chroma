#include "chromabase.h"
#include "update/molecdyn/integrator/md_integrator_factory.h"
#include "update/molecdyn/integrator/lcm_sts_force_gradient_recursive.h"
#include "update/molecdyn/integrator/lcm_exp_sdt.h"
#include "update/molecdyn/integrator/lcm_exp_sdt_sstdt3.h"
#include "io/xmllog_io.h"

#include <string>
using namespace std;

namespace Chroma 
{ 
  
  namespace LatColMatSTSForceGradientRecursiveIntegratorEnv 
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
	LatColMatSTSForceGradientRecursiveIntegratorParams p(xml, path);
    
	return new LatColMatSTSForceGradientRecursiveIntegrator(p);
      }
      
      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "LCM_STS_FORCE_GRADIENT";

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
  
  
  LatColMatSTSForceGradientRecursiveIntegratorParams::LatColMatSTSForceGradientRecursiveIntegratorParams(XMLReader& xml_in, const std::string& path) 
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
	QDPIO::cout << "Subintegrator XML is: " << endl;
	QDPIO::cout << subintegrator_xml << endl;
      }
    }
    catch ( const std::string& e ) { 
      QDPIO::cout << "Error reading XML in LatColMatSTSForceGradientRecursiveIntegratorParams " << e << endl;
      QDP_abort(1);
    }
  }
  
  void read(XMLReader& xml, 
	    const std::string& path, 
	    LatColMatSTSForceGradientRecursiveIntegratorParams& p) {
    LatColMatSTSForceGradientRecursiveIntegratorParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, 
	     const std::string& path, 
	     const LatColMatSTSForceGradientRecursiveIntegratorParams& p) {
    push(xml, path);
    write(xml, "n_steps", p.n_steps);
    write(xml, "monomial_ids", p.monomial_ids);
    xml << p.subintegrator_xml;
    pop(xml);
  }

    void LatColMatSTSForceGradientRecursiveIntegrator::getShadow( AbsFieldState<multi1d<LatticeColorMatrix>,
								  multi1d<LatticeColorMatrix> >& s,
								  const Real& dt,
								  Double& H, 
								  Double& Hs) const
  {
    multi1d<Poisson> pb(monomials.size());
    multi1d<Double> shadow(monomials.size());
    
    for(int i =0; i < monomials.size(); i++) {
      ExactMonomial< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& the_monomial
	= dynamic_cast<ExactMonomial< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >&>(*monomials[i]);

      TowerArray<LatticeColorMatrix> F(4);
      TowerArray<LatticeColorMatrix> G(2);
      
      the_monomial.computePBVectorField(F,G,s);
      Poisson p(F,G,s.getP());

      pb[i] = p;
      shadow[i] = Double(0);
      shadow -= Real(36)/Real(155520)*real(pb[i].stsst);
      shadow -= Real(72)/Real(155520)*real(pb[i].sttst);
      shadow -= Real(126.0)/Real(155520)*real(pb[i].ttsst);
      shadow -= Real(54.0)/Real(155520)*real(pb[i].tttst);



      //      shadow += Real(251.0/737280.0)*real(pb[i].sssst);
      // shadow += Real(13.0/15360.0)*real(pb[i].tssst);
    }

    H=Double(0);

    // Kinetic Piece
    for(int mu=0; mu < Nd; mu++) { 
      H += norm2(s.getP()[mu]);
    }
    
    // Action Piece
    for(int i=0; i < monomials.size(); i++) { 
      ExactMonomial< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& the_monomial
	= dynamic_cast<ExactMonomial< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >&>(*monomials[i]);

      H += the_monomial.S(s);
    }

    // Shadow piece
    Hs = H;
    for(int i=0; i < monomials.size(); i++) { 
      Hs += shadow[i]*dt*dt*dt*dt;
    }


  }

  void LatColMatSTSForceGradientRecursiveIntegrator::operator()( 
					     AbsFieldState<multi1d<LatticeColorMatrix>,
					     multi1d<LatticeColorMatrix> >& s, 
					     const Real& traj_length) const
  {
   
    START_CODE();
    LatColMatExpSdtIntegrator expSdt(1,
				     monomials);

    LatColMatExpSdtMinusSSTdt3Integrator expSdtSSTdt3(1,
						      monomials);


    const AbsComponentIntegrator< multi1d<LatticeColorMatrix>,
      multi1d<LatticeColorMatrix> >& subIntegrator = getSubIntegrator();

    				    
    Real dtau = traj_length / Real(n_steps);
    Real dtau_1_6 = dtau / Real(6);
    Real dtau_1_2 = dtau / Real(2);
    Real dtau_48_72 = Real(48)*dtau/Real(72);
    Real m_dtau_cubed_72 = -(dtau*dtau*dtau)/Real(72);

    // It's sts so:
    Double H_0, Hs_0;
    Double H, Hs;
    Double rms_dH=0;
    Double rms_dHs=0;
    getShadow(s,dtau,H_0,Hs_0) ;

    for(int i=0; i < n_steps; i++) { 
      expSdt(s, dtau_1_6);
      subIntegrator(s, dtau_1_2);
      expSdt(s, dtau_48_72);
      expSdtSSTdt3(s, m_dtau_cubed_72);
      subIntegrator(s, dtau_1_2);
      expSdt(s, dtau_1_6);

      getShadow(s,dtau,H,Hs);
      QDPIO::cout << "FORCE_GRADIENT: " << H-H_0 << " " << Hs - Hs_0 << endl;
      rms_dH += (H-H_0)*(H-H_0);
      rms_dHs += (Hs-Hs_0)*(Hs -Hs_0);
    }

    rms_dH /= Double(n_steps);
    rms_dHs /= Double(n_steps);

    QDPIO::cout << "%FORCE_GRADIENT: log(dt) = " << log10(dtau) << " log dH = " << log10(sqrt(rms_dH)) << " log dHs = " << log10(sqrt(rms_dHs)) << endl;
    END_CODE();
    

  }


};

#include "chromabase.h"
#include "update/molecdyn/integrator/md_integrator_factory.h"
#include "update/molecdyn/integrator/lcm_pqp_leapfrog.h"
#include "io/xmllog_io.h"

#include <string>
#include "util/gauge/taproj.h"
#include "util/gauge/reunit.h"
#include "util/gauge/expmat.h"


namespace Chroma { 
  
  namespace LatColMatPQPLeapfrogIntegratorEnv {

    AbsMDIntegrator<multi1d<LatticeColorMatrix>, 
		    multi1d<LatticeColorMatrix> >* createMDIntegrator(
								       XMLReader& xml, const std::string& path,  Handle< AbsHamiltonian<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >& H) {
      // Read the integrator params
      LatColMatPQPLeapfrogIntegratorParams p(xml, path);
    
      return new LatColMatPQPLeapfrogIntegrator(p, H);
    }

    const std::string name = "LCM_PQP_LEAPFROG_INTEGRATOR";
    const bool registered = TheMDIntegratorFactory::Instance().registerObject(name, createMDIntegrator); 
  };

  LatColMatPQPLeapfrogIntegratorParams::LatColMatPQPLeapfrogIntegratorParams(XMLReader& xml_in, const std::string& path) 
  {
    XMLReader paramtop(xml_in, path);
    try {
      read(paramtop, "./dt", dt);
      read(paramtop, "./tau0", tau0);
    }
    catch ( const std::string& e ) { 
      QDPIO::cout << "Error reading XML in LatColMatPQPLeapfrogIntegratorParams " << e << endl;
      QDP_abort(1);
    }
  }

  void read(XMLReader& xml, 
	    const std::string& path, 
	    LatColMatPQPLeapfrogIntegratorParams& p) {
    LatColMatPQPLeapfrogIntegratorParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, 
	     const std::string& path, 
	     const LatColMatPQPLeapfrogIntegratorParams& p) {
    push(xml, path);
    write(xml, "dt", p.dt);
    write(xml, "tau0", p.tau0);
    pop(xml);
  }

#if 0
  void LatColMatPQPLeapfrogIntegrator::leapP(const Real& dt, 
					     AbsFieldState<multi1d<LatticeColorMatrix>,
					     multi1d<LatticeColorMatrix> >& s) {

    AbsHamiltonian<multi1d<LatticeColorMatrix>,
      multi1d<LatticeColorMatrix> >& H = getHamiltonian();

    XMLWriter& xml_out = TheXMLOutputWriter::Instance();
    // Self Description rule
    push(xml_out, "leapP");
    write(xml_out, "dt",dt);

      // Force Term
    multi1d<LatticeColorMatrix> dsdQ(Nd);
    
    // Compute the force 
    H.dsdq(dsdQ, s);
    
    // Zero boundaries ? -- where would this be done then?
    // There is a zero boundary in GaugeBC?
    // H.zero(s.getP());
    
    // This should be done in a one liner..
    // a la s.getP() -= eps*dsdQ;
    // doing it in loop for now
    
    for(int mu =0; mu < Nd; mu++) { 
      (s.getP())[mu] += dt * dsdQ[mu];
      
      // taproj it...
      taproj( (s.getP())[mu] );
    }
    
    pop(xml_out); // pop("leapP");
  }

  //! Leap with Q
  void LatColMatPQPLeapfrogIntegrator::leapQ(const Real& dt, 
					     AbsFieldState<multi1d<LatticeColorMatrix>,
					     multi1d<LatticeColorMatrix> >& s) {
    LatticeColorMatrix tmp_1;
    LatticeColorMatrix tmp_2;

    XMLWriter& xml_out= TheXMLOutputWriter::Instance();
    // Self description rule
    push(xml_out, "leapQ");
    write(xml_out, "dt", dt);

    // Constant
    const multi1d<LatticeColorMatrix>& p_mom = s.getP();
    
    // Mutable
    multi1d<LatticeColorMatrix>& u = s.getQ();
    
    for(int mu = 0; mu < Nd; mu++) { 
      
      //  dt*p[mu]
      tmp_1 = dt*(s.getP())[mu];
      
      // tmp_1 = exp(dt*p[mu])  
      expmat(tmp_1, EXP_TWELFTH_ORDER);
      
      // tmp_2 = exp(dt*p[mu]) u[mu] = tmp_1 * u[mu]
      tmp_2 = tmp_1*(s.getQ())[mu];
      
      // u[mu] =  tmp_1 * u[mu] =  tmp_2 
      (s.getQ())[mu] = tmp_2;
      
      // Reunitarize u[mu]
      int numbad;
      reunit((s.getQ())[mu], numbad, REUNITARIZE_ERROR);
    }

    pop(xml_out);
  }
#endif

};

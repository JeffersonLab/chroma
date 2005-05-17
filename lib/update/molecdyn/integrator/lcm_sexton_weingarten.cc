#include "chromabase.h"
#include "update/molecdyn/integrator/md_integrator_factory.h"
#include "update/molecdyn/integrator/lcm_sexton_weingarten.h"
#include "io/xmllog_io.h"

#include <string>
#include "util/gauge/taproj.h"
#include "util/gauge/reunit.h"
#include "util/gauge/expmat.h"


namespace Chroma { 
  
  namespace LatColMatSextonWeingartenIntegratorEnv {

    AbsMDIntegrator<multi1d<LatticeColorMatrix>, 
		    multi1d<LatticeColorMatrix> >* createMDIntegrator(
								       XMLReader& xml, const std::string& path,  Handle< AbsHamiltonian<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >& H) {
      // Read the integrator params
      LatColMatSextonWeingartenIntegratorParams p(xml, path);
    
      return new LatColMatSextonWeingartenIntegrator(p, H);
    }

    const std::string name = "LCM_SEXTON_WEINGARTEN_INTEGRATOR";
    const bool registered = TheMDIntegratorFactory::Instance().registerObject(name, createMDIntegrator); 
  };

  LatColMatSextonWeingartenIntegratorParams::LatColMatSextonWeingartenIntegratorParams(XMLReader& xml_in, const std::string& path) 
  {
    XMLReader paramtop(xml_in, path);
    try {
      read(paramtop, "./n_steps", n_steps);
      read(paramtop, "./tau0", tau0);
      read(paramtop, "./n_short_steps", n_short_steps);
      read(paramtop, "./S_short_monomials", S_short_monomials);
      read(paramtop, "./S_long_monomials", S_long_monomials);
    
    }
    catch ( const std::string& e ) { 
      QDPIO::cout << "Error reading XML in LatColMatSextonWeingartenIntegratorParams " << e << endl;
      QDP_abort(1);
    }
  }

  void read(XMLReader& xml, 
	    const std::string& path, 
	    LatColMatSextonWeingartenIntegratorParams& p) {
    LatColMatSextonWeingartenIntegratorParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, 
	     const std::string& path, 
	     const LatColMatSextonWeingartenIntegratorParams& p) {
    push(xml, path);
    write(xml, "n_steps", p.n_steps);
    write(xml, "tau0", p.tau0);
    write(xml, "n_short_steps", p.n_short_steps);
    write(xml, "S_short_monomials", p.S_short_monomials);
    write(xml, "S_long_monomials", p.S_long_monomials);
    pop(xml);
  }

  void LatColMatSextonWeingartenIntegrator::leapP(const multi1d<int>& monomial_list,
						  const Real& dt, 
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
    H.dsdq(dsdQ, s, monomial_list);
    
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
  void LatColMatSextonWeingartenIntegrator::leapQ(const Real& dt, 
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

};

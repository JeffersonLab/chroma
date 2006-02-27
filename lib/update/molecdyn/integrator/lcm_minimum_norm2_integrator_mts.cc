// $Id: lcm_minimum_norm2_integrator_mts.cc,v 2.1 2006-02-27 15:34:39 bjoo Exp $

#include "chromabase.h"
#include "update/molecdyn/integrator/md_integrator_factory.h"
#include "update/molecdyn/integrator/lcm_minimum_norm2_integrator_mts.h"
#include "io/xmllog_io.h"

#include <string>
#include "util/gauge/taproj.h"
#include "util/gauge/reunit.h"
#include "util/gauge/expmat.h"


namespace Chroma { 
  
  namespace LatColMatMinimumNorm2IntegratorMtsEnv {
    
    AbsMDIntegrator<multi1d<LatticeColorMatrix>, 
		    multi1d<LatticeColorMatrix> >* createMDIntegrator(
								      XMLReader& xml, const std::string& path,  Handle< AbsHamiltonian<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >& H) {
      // Read the integrator params
      LatColMatMinimumNorm2IntegratorMtsParams p(xml, path);
      
      return new LatColMatMinimumNorm2IntegratorMts(p, H);
    }

    const std::string name = "LCM_MINIMUM_NORM_2ND_ORDER_INTEGRATOR_MTS";
    const bool registered = TheMDIntegratorFactory::Instance().registerObject(name, createMDIntegrator); 
  };

  LatColMatMinimumNorm2IntegratorMtsParams::LatColMatMinimumNorm2IntegratorMtsParams(XMLReader& xml_in, 
									   const std::string& path) 
  {
    XMLReader paramtop(xml_in, path);
    try {
      read(paramtop, "./number_of_timescales", number_of_timescales);
      read(paramtop, "./tau", tau);
      read(paramtop, "./lambda_list", lambda_list);
      read(paramtop, "./n_steps_list", n_steps_list);
      read(paramtop, "./monomial_list", monomial_list);
    }
    catch ( const std::string& e ) { 
      QDPIO::cout << "Error reading XML in LatColMatMinimumNorm2IntegratorMtsParams " << e << endl;
      QDP_abort(1);
    }
  }

  void read(XMLReader& xml, 
	    const std::string& path, 
	    LatColMatMinimumNorm2IntegratorMtsParams& p) {
    LatColMatMinimumNorm2IntegratorMtsParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, 
	     const std::string& path, 
	     const LatColMatMinimumNorm2IntegratorMtsParams& p) {
    push(xml, path);
    write(xml, "n_steps_list", p.n_steps_list);
    write(xml, "tau", p.tau);
    write(xml, "lambda_list", p.lambda_list);
    write(xml, "number_of_timescales", p.number_of_timescales);
    write(xml, "monomial_list", p.monomial_list);
    pop(xml);
  }

#if 0
  void LatColMatMinimumNorm2IntegratorMts::leapP(const multi1d<int>& monomial_list,
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
  void LatColMatMinimumNorm2IntegratorMts::leapQ(const Real& dt, 
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

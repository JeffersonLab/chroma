#include "lcm_integrator_leaps.h"
#include "util/gauge/taproj.h"
#include "util/gauge/reunit.h"
#include "util/gauge/expmat.h"

namespace Chroma { 

  namespace LCMMDIntegratorSteps { 

        //! LeapP for just a selected list of monomials
    void leapP(const Real& dt, 

	       AbsHamiltonian<multi1d<LatticeColorMatrix>,
	                      multi1d<LatticeColorMatrix> >& H,

	       AbsFieldState<multi1d<LatticeColorMatrix>,
	       multi1d<LatticeColorMatrix> >& s) 
    {

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
    
    
    //! LeapP for just a selected list of monomials
    void leapP(const multi1d<int>& monomial_list,
	       const Real& dt, 
	       
	       AbsHamiltonian<multi1d<LatticeColorMatrix>,
	       multi1d<LatticeColorMatrix> >& H,
	       
	       AbsFieldState<multi1d<LatticeColorMatrix>,
	       multi1d<LatticeColorMatrix> >& s) {
      
      
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

    void leapQ(const Real& dt, 
	       AbsFieldState<multi1d<LatticeColorMatrix>,
	       multi1d<LatticeColorMatrix> >& s) 
    {

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
	// expmat(tmp_1, EXP_TWELFTH_ORDER);
	expmat(tmp_1, EXP_EXACT);
	
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
};

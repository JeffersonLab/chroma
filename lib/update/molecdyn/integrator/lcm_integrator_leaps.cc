#include "lcm_integrator_leaps.h"
#include "util/gauge/taproj.h"
#include "util/gauge/reunit.h"
#include "util/gauge/expmat.h"
#include "update/molecdyn/monomial/force_monitors.h"

namespace Chroma 
{ 

  namespace LCMMDIntegratorSteps 
  { 

    //! LeapP for just a selected list of monomials
    void leapP(const multi1d< IntegratorShared::MonomialPair >& monomials,
	                                       
	       const Real& dt, 
	       	       
	       AbsFieldState<multi1d<LatticeColorMatrix>,
	       multi1d<LatticeColorMatrix> >& s)
    {
      START_CODE();
      StopWatch swatch;


      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      // Self Description rule
      push(xml_out, "leapP");
      write(xml_out, "dt", dt);
      multi1d<Real> real_step_size(Nd);

      // Work out the array of step sizes (including all scaling factors
      for(int mu =0; mu < Nd; mu++) 
      {
	real_step_size[mu] = dt * theAnisoStepSizeArray::Instance().getStepSizeFactor(mu);
      }
      write(xml_out, "dt_actual_per_dir", real_step_size);
      
      // Force Term
      multi1d<LatticeColorMatrix> dsdQ(Nd);
      push(xml_out, "AbsHamiltonianForce"); // Backward compatibility
      write(xml_out, "num_terms", monomials.size());
      push(xml_out, "ForcesByMonomial");

      if( monomials.size() > 0 ) { 
	push(xml_out, "elem");
	swatch.reset(); swatch.start();
	monomials[0].mon->dsdq(dsdQ,s);
	swatch.stop();
	QDPIO::cout << "FORCE TIME: " << monomials[0].id <<  " : " << swatch.getTimeInSeconds() << endl;
	pop(xml_out); //elem
	for(int i=1; i < monomials.size(); i++) { 
	  push(xml_out, "elem");
	  multi1d<LatticeColorMatrix> cur_F(Nd);
	  swatch.reset(); swatch.start();
	  monomials[i].mon->dsdq(cur_F, s);
	  swatch.stop();
	  dsdQ += cur_F;

	  QDPIO::cout << "FORCE TIME: " << monomials[i].id << " : " << swatch.getTimeInSeconds() << "\n";
 
	  pop(xml_out); // elem
	}
      }
      pop(xml_out); // ForcesByMonomial
      //monitorForces(xml_out, "TotalForcesThisLevel", dsdQ);
      pop(xml_out); // AbsHamiltonianForce 


      for(int mu =0; mu < Nd; mu++) {

	(s.getP())[mu] += real_step_size[mu] * dsdQ[mu];
	
	// taproj it...
	taproj( (s.getP())[mu] );
      }
      
      pop(xml_out); // pop("leapP");
    
      END_CODE();
    }

    void leapQ(const Real& dt, 
	       AbsFieldState<multi1d<LatticeColorMatrix>,
	       multi1d<LatticeColorMatrix> >& s) 
    {
      START_CODE();

      LatticeColorMatrix tmp_1;
      LatticeColorMatrix tmp_2;
      
      XMLWriter& xml_out= TheXMLLogWriter::Instance();
      // Self description rule
      push(xml_out, "leapQ");
      write(xml_out, "dt", dt);
      multi1d<Real> real_step_size(Nd);

      // Work out the array of step sizes (including all scaling factors
      for(int mu =0; mu < Nd; mu++) 
      {
	real_step_size[mu] = dt * theAnisoStepSizeArray::Instance().getStepSizeFactor(mu);
      }
      write(xml_out, "dt_actual_per_dir", real_step_size);
      
      // Constant
      const multi1d<LatticeColorMatrix>& p_mom = s.getP();
      
      // Mutable
      multi1d<LatticeColorMatrix>& u = s.getQ();
      
      for(int mu = 0; mu < Nd; mu++) 
      {
	//  dt*p[mu]
	tmp_1 = real_step_size[mu]*(s.getP())[mu];
	
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
    
      END_CODE();
    }

  }
}

#include "lcm_integrator_leaps.h"
#include "util/gauge/taproj.h"
#include "util/gauge/reunit.h"
#include "util/gauge/expmat.h"
#include "update/molecdyn/monomial/force_monitors.h"
#include "tower.h"
namespace Chroma 
{ 

  namespace LCMMDIntegratorSteps 
  { 

    //! LeapP for just a selected list of monomials
    void leapP(const multi1d< 
	       Handle<   Monomial< multi1d<LatticeColorMatrix>, 
	                           multi1d<LatticeColorMatrix> > >
	       > monomials,
	                                       
	       const Real& dt, 
	       	       
	       AbsFieldState<multi1d<LatticeColorMatrix>,
	       multi1d<LatticeColorMatrix> >& s) 
    {
      START_CODE();
      
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
      multi1d<Tower<LatticeColorMatrix> > dsdQ(Nd);
      multi1d<Tower<LatticeColorMatrix> > cur_F(Nd);

      // Just need a force term so towers of height 1 are adequate
      for(int mu=0; mu < Nd; mu++) { 
	dsdQ[mu].resize(1);
	cur_F[mu].resize(1);
      }

      push(xml_out, "AbsHamiltonianForce"); // Backward compatibility
      write(xml_out, "num_terms", monomials.size());
      push(xml_out, "ForcesByMonomial");

      if( monomials.size() > 0 ) { 
	push(xml_out, "elem");
	monomials[0]->dsdq(dsdQ,s);
	pop(xml_out); //elem
	for(int i=1; i < monomials.size(); i++) { 
	  push(xml_out, "elem");
	  
	  monomials[i]->dsdq(cur_F, s);
	  for(int mu=0; mu < Nd; mu++) { 
	    dsdQ[mu] += cur_F[mu];
	  }
	  pop(xml_out); // elem
	}
      }
      pop(xml_out); // ForcesByMonomial
      //monitorForces(xml_out, "TotalForcesThisLevel", dsdQ);
      pop(xml_out); // AbsHamiltonianForce 


      for(int mu =0; mu < Nd; mu++) {

	(s.getP())[mu] += real_step_size[mu] * dsdQ[mu][0];
	
	// taproj it...
	taproj( (s.getP())[mu] );
      }
      
      pop(xml_out); // pop("leapP");
    
      END_CODE();
    }

    //! LeapP for just a selected list of monomials
    void leapFG(const multi1d< 
		Handle<   Monomial< multi1d<LatticeColorMatrix>, 
		multi1d<LatticeColorMatrix> > >
		> monomials,
		
		const Real& dt, 
		AbsFieldState<multi1d<LatticeColorMatrix>,
		multi1d<LatticeColorMatrix> >& s) 
    {
      START_CODE();
      
      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      // Self Description rule
      push(xml_out, "leapFG");
      write(xml_out, "dt", dt);

      multi1d<Real> real_step_size(Nd);

      // Work out the array of step sizes (including all scaling factors
      for(int mu =0; mu < Nd; mu++) 
      {
	Real factor = theAnisoStepSizeArray::Instance().getStepSizeFactor(mu);
	real_step_size[mu] = dt * factor*factor*factor;
      }
      write(xml_out, "dt_actual_per_dir", real_step_size);
      
      // Force Gradient Term
      multi1d<Tower<LatticeColorMatrix> > G(Nd);

      {
	multi1d<Tower<LatticeColorMatrix> > cur_F(Nd);
	multi1d<Tower<LatticeColorMatrix> > cur_F2(Nd);

	// Just need a force term so towers of height 1 are adequate
	for(int mu=0; mu < Nd; mu++) { 
	  G[mu].resize(2);
	  
	  cur_F[mu].resize(1);
	  cur_F2[mu].resize(2);

	}
	
	push(xml_out, "AbsHamiltonianForce"); // Backward compatibility
	write(xml_out, "num_terms", monomials.size());
	push(xml_out, "ForcesByMonomial");
	
	if( monomials.size() > 0 ) { 
	  push(xml_out, "elem");

	  monomials[0]->dsdq(cur_F,s);  // Height 1

	  // Pass this in, so it can be lifted from
	  for(int mu=0; mu < Nd; mu++) { 
	    G[mu][0]=cur_F[mu][0];
	  }
	  monomials[0]->dsdq(G, s);  // Height 2
	  pop(xml_out); //elem



	  for(int i=1; i < monomials.size(); i++) { 
	    push(xml_out, "elem");
	    
	    monomials[i]->dsdq(cur_F, s);  // Height 1

	    // Pass this in so it can be lifted from 
	    for(int mu=0; mu < Nd; mu++) { 
	      cur_F2[mu][0] = cur_F[mu][0];
	    }

	    monomials[i]->dsdq(cur_F2,s); // Height 2

	    // Accumulate
	    for(int mu=0; mu < Nd; mu++) { 
	      G[mu] += cur_F2[mu];
	    }
	    pop(xml_out); // elem
	  }
	}
      }

      


      pop(xml_out); // ForcesByMonomial
      //monitorForces(xml_out, "TotalForcesThisLevel", dsdQ);
      pop(xml_out); // AbsHamiltonianForce 


      for(int mu =0; mu < Nd; mu++) {

	(s.getP())[mu] += Real(2)*real_step_size[mu]* G[mu][1];

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

#include "lcm_integrator_leaps.h"
#include "util/gauge/taproj.h"
#include "util/gauge/reunit.h"
#include "util/gauge/expmat.h"
#include "update/molecdyn/monomial/force_monitors.h"

#include "actions/boson/operator/adjoint_derivative.h"
#include "actions/boson/invert/invcg_adj.h"

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
      
      int dir =  LCMMDIntegratorSteps::theHyperPlane::Instance().getDir();
      bool active = LCMMDIntegratorSteps::theHyperPlane::Instance().active();

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
	QDPIO::cout << "FORCE TIME: " << monomials[0].id <<  " : " << swatch.getTimeInSeconds() << std::endl;
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


      for(int mu =0; mu < Nd; mu++)
	if((active && (mu==dir))||(!active)){

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
      int dir =  LCMMDIntegratorSteps::theHyperPlane::Instance().getDir();
      bool active = LCMMDIntegratorSteps::theHyperPlane::Instance().active();
      Real rho(LCMMDIntegratorSteps::theHyperPlane::Instance().getRho());

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

      if(!active){
	for(int mu = 0; mu < Nd; mu++){
	  
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
      }
      else{ // Now do the special thing for the MG
	
	QDPIO::cout<<" leapQ: leaping Q for dir= "<<dir
		   <<" and rho= "<<rho<<std::endl ;
	
	AdjointDerivative D(dir,rho,s.getQ());
	Handle< squaredAdjointDerivative> D2 = new squaredAdjointDerivative(D) ;
	
	
	Real RsdCG = 1.0e-7 ; // hard coded for now
	int MaxCG = 5000 ; // hard coded for now
	
	LatticeColorMatrix P = s.getP()[dir] ;
	LatticeColorMatrix X = zero ;
	SystemSolverResults_t res = InvCG_adj(*D2,P,X,RsdCG,MaxCG) ;

	//  dt*p[mu]
	tmp_1 = real_step_size[dir]*X;
	
	// tmp_1 = exp(dt*p[mu])  
	// expmat(tmp_1, EXP_TWELFTH_ORDER);
	expmat(tmp_1, EXP_EXACT);
	
	// tmp_2 = exp(dt*p[mu]) u[mu] = tmp_1 * u[mu]
	tmp_2 = tmp_1*(s.getQ())[dir];
	
	// u[mu] =  tmp_1 * u[mu] =  tmp_2 
	(s.getQ())[dir] = tmp_2;
	
	// Reunitarize u[mu]
	int numbad;
	reunit((s.getQ())[dir], numbad, REUNITARIZE_ERROR);
	QDPIO::cout<<" reunit: NUMBAD="<<numbad<<std::endl ;
      }
	
      pop(xml_out);
      
      END_CODE();
    }

  }
}

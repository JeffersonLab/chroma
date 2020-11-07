// -*- C++ -*-
/*! \file
 * \brief Abstract HMC trajectory Using the new structure
 *
 * HMC trajectories
 */
#ifndef abs_hmc_new_h
#define abs_hmc_new_h

#include "chromabase.h"
#include "io/xmllog_io.h"
#include "update/molecdyn/field_state.h"
#include "update/molecdyn/hamiltonian//abs_hamiltonian.h"
#include "update/molecdyn/integrator/abs_integrator.h"
#include "update/molecdyn/hmc/global_metropolis_accrej.h"
#include "actions/ferm/invert/mg_solver_exception.h"


namespace Chroma 
{ 

  //! Abstract HMC trajectory
  /*! @ingroup hmc */
  template<typename P, typename Q>
  class AbsHMCTrj {
  public: 
    
    // Virtual destructor
    virtual ~AbsHMCTrj() {};
    

    // Do the HMC trajectory
    virtual void operator()(AbsFieldState<P,Q>& s,
			    const bool WarmUpP, 
			    const bool CheckRevP)
    {
      START_CODE();
      StopWatch swatch;


      AbsMDIntegrator<P,Q>& MD = getMDIntegrator();
      AbsHamiltonian<P,Q>& H_MC = getMCHamiltonian();

      XMLWriter& xml_out = TheXMLOutputWriter::Instance();
      XMLWriter& xml_log = TheXMLLogWriter::Instance();

      // Self encapsulation rule 
      push(xml_out, "HMCTrajectory");
      push(xml_log, "HMCTrajectory");


      write(xml_out, "WarmUpP", WarmUpP);
      write(xml_log, "WarmUpP", WarmUpP);
      
      // HMC Algorithm.
      // 1) Refresh momenta
      //
      swatch.reset(); swatch.start();
      refreshP(s);
      swatch.stop();
      QDPIO::cout << "HMC_TIME: Momentum Refresh Time: " << swatch.getTimeInSeconds() << " \n";

      bool acceptTraj = false; // Default value. Acceptance check can change this
      Double KE_old, PE_old;


      // Try catch block in case MG Solver fails in any pseudofermion refresh.
      // That should cause an abort.
      try { 
	// Refresh Pseudofermions
	swatch.reset(); swatch.start();
	H_MC.refreshInternalFields(s);
	swatch.stop();
	QDPIO::cout << "HMC_TIME: Pseudofermion Refres Time: " << swatch.getTimeInSeconds() << " \n";
      }
      catch ( MGSolverException e ) {
	QDPIO::cout << "ERROR: Caught MG Solver exception in pseudofermion refresh" << std::endl;
	QDPIO::cout << "ERROR: Exception Was: " << e.whatStr() << std::endl;
	QDPIO::cout << "Aborting";
	QDP_abort(2);
      }

      // SaveState -- Perhaps this could be done better?
      Handle< AbsFieldState<P,Q> >  s_old(s.clone());

      
      // Try Catch block in case MG Solver fails in Energy Calculation
      // That should also be an abort
      try {
	// Measure energy of the old state

	push(xml_out, "H_old");
	push(xml_log, "H_old");
	swatch.reset(); swatch.start();
	H_MC.mesE(*s_old, KE_old, PE_old);
	swatch.stop();
	QDPIO::cout << "HMC_TIME: Start Energy Time: " << swatch.getTimeInSeconds() << " \n";
	
	write(xml_out, "KE_old", KE_old);
	write(xml_log, "KE_old", KE_old);
	
	write(xml_out, "PE_old", PE_old);
	write(xml_log, "PE_old", PE_old);
	
	pop(xml_log); // pop H_old
	pop(xml_out); // pop H_old
      }
      catch( MGSolverException e ) { 
	QDPIO::cout << "ERROR: Caught MG Solver exception in Start Energy Calculation" << std::endl;
	QDPIO::cout << "ERROR: Exception Was: " << e.whatStr() << std::endl;
	QDPIO::cout << "Aborting";
	QDP_abort(2);
      }
      
      // Try-Catch Case for MD
      // If solver fails in MD, we reject by setting acceptTraj = false
      // If traj is sccessful (no exception is thrown) acceptTraj = true for WarmUp traj and 
      //    is decided by Accept/Reject step otherwise

      try {	
	swatch.start();
      
	// Copy in fields from the Hamiltonian as needed using the
	// CopyList
	MD.copyFields();

	// Integrate MD trajectory
	MD(s, MD.getTrajLength());
	swatch.stop();
	QDPIO::cout << "HMC_TIME: Traj MD Time: " << swatch.getTimeInSeconds() << " \n";
           
	// Measure the energy of the new state - before reverse check
	// This is because if an MG preconditioner is used, the reverse
	// traj may change it to where it is a poor preconditioner for the 
	// Final energy calculation
	swatch.reset(); 
	swatch.start();
	Double KE, PE;
	push(xml_out, "H_new");
	push(xml_log, "H_new");
	H_MC.mesE(s, KE, PE);
	write(xml_out, "KE_new", KE);
	write(xml_log, "KE_new", KE);
	write(xml_out, "PE_new", PE);
	write(xml_log, "PE_new", PE);
	pop(xml_log);
	pop(xml_out);
	swatch.stop();
	QDPIO::cout << "HMC_TIME: Finish Energy Time: " << swatch.getTimeInSeconds() << " \n";
	// Work out energy differences
	Double DeltaKE = KE - KE_old;
	Double DeltaPE = PE - PE_old;
	Double DeltaH  = DeltaKE + DeltaPE;
	Double AccProb = where(DeltaH < 0.0, Double(1), exp(-DeltaH));
	write(xml_out, "deltaKE", DeltaKE);
	write(xml_log, "deltaKE", DeltaKE);
	
	write(xml_out, "deltaPE", DeltaPE);
	write(xml_log, "deltaPE", DeltaPE);
	
	write(xml_out, "deltaH", DeltaH);
	write(xml_log, "deltaH", DeltaH);
	
	write(xml_out, "AccProb", AccProb);
	write(xml_log, "AccProb", AccProb);
	
	QDPIO::cout << "Delta H = " << DeltaH << std::endl;
	QDPIO::cout << "AccProb = " << AccProb << std::endl;

	// If we intend to do an accept reject step
	// (ie we are not warming up)
	if( !WarmUpP ) {
	 
	  // Measure Acceptance
	  acceptTraj = acceptReject(DeltaH);
	}
	else {
	  // Warm Up: acceptTraj = true -- always
	  acceptTraj = true; 
	}

      }       
      catch( MGSolverException e ) {

	// Exception Handling if the solver didnt converge
	// Automatic rejection
	QDPIO::cout << "WARNING: Caught MG Solver Convergence Exception, during MD or Final Energy Calculation!" << std::endl;
	QDPIO::cout << "WARNING: Exception was: " << e.whatStr() << std::endl;
	QDPIO::cout << "WARNING: Aborting" << std::endl;
	QDP_abort(2);
      }

      write(xml_out, "AcceptP", acceptTraj);
      write(xml_log, "AcceptP", acceptTraj);	  
      QDPIO::cout << "AcceptP = " << acceptTraj << std::endl;


      // If reversebility check is due, do it. It doesn't affect the final state 's'
      if( CheckRevP ) {

	// Copy State
	Handle< AbsFieldState<P,Q> >  s_rev(s.clone());

	// Try and do reverse trajectory 
	// If the MG solver fails in this, in some sense I don't care
	// It doesn't affect acceptance of 's'

	try { 
	  swatch.reset(); swatch.start();
	  
	  QDPIO::cout << "Reversing trajectory for reversability test" <<std::endl;
	
	  // Flip Momenta
	  flipMomenta(*s_rev);
	
	  // Go back
	  MD(*s_rev, MD.getTrajLength());

	  // Flip Momenta back (to original)
	  flipMomenta(*s_rev);


	  Double KE_rev;
	  Double PE_rev;
       
	  H_MC.mesE(*s_rev, KE_rev, PE_rev);

	  Double DeltaDeltaKE = KE_rev - KE_old;
	  Double DeltaDeltaPE = PE_rev - PE_old;
	  Double DeltaDeltaH = DeltaDeltaKE + DeltaDeltaPE;

	
	  Double dq;
	  Double dp;
	  reverseCheckMetrics(dq,dp, *s_rev, *s_old);

	  push(xml_log, "ReversibilityMetrics");
	  write(xml_log, "DeltaDeltaH", fabs(DeltaDeltaH));
	  write(xml_log, "DeltaDeltaKE", fabs(DeltaDeltaKE));
	  write(xml_log, "DeltaDeltaPE", fabs(DeltaDeltaPE));
	  write(xml_log, "DeltaQPerSite", dq);
	  write(xml_log, "DeltaPPerSite", dp);
	  pop(xml_log);
	  
	  QDPIO::cout << "Reversibility: DeltaDeltaH = " << fabs(DeltaDeltaH) <<std::endl;
	  QDPIO::cout << "Reversibility: DeltaQ      = " << dq << std::endl;
	  QDPIO::cout << "Reversibility: DeltaP      = " << dp << std::endl;
	  swatch.stop();
	  QDPIO::cout << "HMC_TIME: Reverse Check Time: " << swatch.getTimeInSeconds() << " \n";
	}
	catch( MGSolverException e ) { 

	  QDPIO::cout << "WARNING: Caught MG Solver Exception in Reverse Trajectory" << std::endl;
	  QDPIO::cout << "WARNING: Exception was: " << e.whatStr() << std::endl;
	  QDPIO::cout << "WARNING: Aborting" << std::endl;
	  QDP_abort(2);
	}

	// s_rev goes away here.

      } // if (CheckRevP)

      // Rejection: acceptTraj is false
      // Either because of exception handling, or 
      // because of the MC Test
      if ( ! acceptTraj ) {

	// restore the old state
	s.getQ() = s_old->getQ();
	s.getP() = s_old->getP();
	
      }
            
      pop(xml_log); // HMCTrajectory
      pop(xml_out); // HMCTrajectory
    
      END_CODE();
    }
    
  protected:
    // Get at the Exact Hamiltonian
    virtual AbsHamiltonian<P,Q>& getMCHamiltonian(void) = 0;
    
    // Get at the TopLevelIntegrator
    virtual AbsMDIntegrator<P,Q>& getMDIntegrator(void)  = 0;
    
    virtual void refreshP(AbsFieldState<P,Q>& state) const = 0;
    
    virtual bool acceptReject(const Double& DeltaH) const = 0;
    

    // These things are for reversebility checking...
    virtual void flipMomenta(AbsFieldState<P,Q>& state) const = 0;
    virtual void reverseCheckMetrics(Double& deltaQ, Double& deltaP,
				     const AbsFieldState<P,Q>& s, 
				     const AbsFieldState<P,Q>& s_old) const = 0;
  };

} // end namespace chroma 



#endif

// -*- C++ -*-
// $Id: abs_hmc.h,v 3.7 2009-01-30 20:57:51 bjoo Exp $
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
      refreshP(s);
      
      // Refresh Pseudofermions
      H_MC.refreshInternalFields(s);
      
      // SaveState -- Perhaps this could be done better?
      Handle< AbsFieldState<P,Q> >  s_old(s.clone());
      
      // Measure energy of the old state
      Double KE_old, PE_old;

      push(xml_out, "H_old");
      push(xml_log, "H_old");

      H_MC.mesE(*s_old, KE_old, PE_old);

      write(xml_out, "KE_old", KE_old);
      write(xml_log, "KE_old", KE_old);

      write(xml_out, "PE_old", PE_old);
      write(xml_log, "PE_old", PE_old);

      pop(xml_log); // pop H_old
      pop(xml_out); // pop H_old
      
      
      // Copy in fields from the Hamiltonian as needed using the
      // CopyList
      MD.copyFields();

      // Integrate MD trajectory
      MD(s, MD.getTrajLength());
           

      // If this is a reverse trajectory
      if( CheckRevP ) { 

	QDPIO::cout << "Reversing trajectory for reversability test" <<endl;

	// Copy state
	Handle< AbsFieldState<P,Q> >  s_rev(s.clone());
	
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

	QDPIO::cout << "Reversibility: DeltaDeltaH = " << fabs(DeltaDeltaH) <<endl;
	QDPIO::cout << "Reversibility: DeltaQ      = " << dq << endl;
	QDPIO::cout << "Reversibility: DeltaP      = " << dp << endl;


	// s_rev goes away... We continue as if nothing happened
	
      }

      //  Measure the energy of the new state
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

      QDPIO::cout << "Delta H = " << DeltaH << endl;
      QDPIO::cout << "AccProb = " << AccProb << endl;

      // If we intend to do an accept reject step
      // (ie we are not warming up)
      if( ! WarmUpP ) 
      { 
	// Measure Acceptance
	bool acceptTestResult = acceptReject(DeltaH);
	write(xml_out, "AcceptP", acceptTestResult);
	write(xml_log, "AcceptP", acceptTestResult);

	QDPIO::cout << "AcceptP = " << acceptTestResult << endl;

	// If rejected restore fields
	// If accepted no need to do anything
	if ( ! acceptTestResult ) 
	{ 
	  s.getQ() = s_old->getQ();
	  
	  // I am going to refresh these so maybe these copies 
	  //  for the momenta and the pseudofermions
	  // are not really needed...
	  //
	  // restore momenta
	  s.getP() = s_old->getP();
	}
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

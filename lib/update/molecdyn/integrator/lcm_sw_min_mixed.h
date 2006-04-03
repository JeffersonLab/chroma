// -*- C++ -*-
// $Id: lcm_sw_min_mixed.h,v 3.0 2006-04-03 04:59:08 edwards Exp $

/*! @file
 * @brief Sexton Weingarten integrator
 *
 * Leapfrog intregator for HMC
 */

#ifndef LCM_SW_MIN_MIXED_H
#define LCM_SW_MIN_MIXED_H


#include "chromabase.h"
#include "update/molecdyn/hamiltonian/abs_hamiltonian.h"
#include "update/molecdyn/integrator/abs_integrator.h"
#include "update/molecdyn/integrator/lcm_integrator_leaps.h"
#include <string>


namespace Chroma 
{

  /*! @ingroup integrator */
  namespace LatColMatSexWeinMixedIntegratorEnv {
    extern const std::string name;
    extern const bool registered;
  };

  /*! @ingroup integrator */
  struct  LatColMatSexWeinMixedIntegratorParams {
    LatColMatSexWeinMixedIntegratorParams();
    LatColMatSexWeinMixedIntegratorParams(XMLReader& xml, const std::string& path);

    int n_steps;
    Real lambda;
    Real tau0;

    int n_short_steps;

    multi1d<int> S_short_monomials; // Indices of short timescale monomials
                                    // which we do more often

    multi1d<int> S_long_monomials; // Indices of long timescale monomials
                                   // which we do less often
  };

  /*! @ingroup integrator */
  void read(XMLReader& xml_in, 
	    const std::string& path,
	    LatColMatSexWeinMixedIntegratorParams& p);

  /*! @ingroup integrator */
  void write(XMLWriter& xml_out,
	     const std::string& path, 
	     const LatColMatSexWeinMixedIntegratorParams& p);

  //! MD integrator interface for PQP leapfrog
  /*! @ingroup integrator
   *  Specialised to multi1d<LatticeColorMatrix>
   */
  class LatColMatSexWeinMixedIntegrator 
    : public AbsMDIntegrator<multi1d<LatticeColorMatrix>,
                             multi1d<LatticeColorMatrix> > {

  public:

    // Construct from params struct and Hamiltonian
    LatColMatSexWeinMixedIntegrator(const LatColMatSexWeinMixedIntegratorParams& p,
				   Handle< AbsHamiltonian< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >& H_) : n_steps(p.n_steps), lambda(p.lambda), tau0(p.tau0), H_MD(H_), n_short_steps(p.n_short_steps), S_short_monomials(p.S_short_monomials), S_long_monomials(p.S_long_monomials) {

      // Check for sanity
      int num_ham_monomials = getHamiltonian().numMonomials();

      int num_sw_monomials = p.S_short_monomials.size() + S_long_monomials.size();

      // Check that the number of short timescale monomials + 
      //   the number of long_timescale monomials equals the total number of 
      //   monomials in the hamiltonian
      if( num_ham_monomials != num_sw_monomials ) {
	QDPIO::cerr << "Error: Wrong number of monomials specified for SW integrator." << endl;
	QDPIO::cerr << "Hamiltonian has " << num_ham_monomials << endl;
	QDPIO::cerr << "SW has " << S_short_monomials.size() << " short timescale monomials and " << endl;
	QDPIO::cerr << S_long_monomials.size() << " long timescale monomials" << endl;
	QDP_abort(1);
      }

      // Check that a monomial has not been used twice and that 
      // no bogus indices are specified
      multi1d<bool> exclusivity_check(num_ham_monomials);
     
      // Initialise filter -- set all elements to false
      for(int i=0; i < num_ham_monomials; i++) { 
	exclusivity_check[i] = false;
      }

      // Go through list of short timescale monomials
      for(int i=0; i < S_short_monomials.size(); i++) {

	// Check for bogus index
	if( S_short_monomials[i] < 0 || 
	    S_short_monomials[i] >= num_ham_monomials ) { 

	  QDPIO::cerr << "Index out of range. S_short_monomials["<<i<<"] has index " << S_short_monomials[i] << " but max index is " << num_ham_monomials-1 << endl;
	  QDP_abort(1);
	}

	// If index is not bogus, check it has not been used already
	// ie only insert it if exclusivity_check[i] is false
	if( exclusivity_check[S_short_monomials[i]] == true ) { 
	  QDPIO::cerr << "Exclusivity check failed. Monomial " << i << " is listed more than once in list of short monomials" << endl;
	  QDP_abort(1);
	}
	else { 
	  // Its not there yet, so insert it
	  exclusivity_check[S_short_monomials[i]] = true;
	}
      }

      // Now go through list of long timescale monomials
      for(int i=0; i < S_long_monomials.size(); i++) {
	// Check for bogus indices
	if( S_long_monomials[i] < 0 || S_long_monomials[i] >= num_ham_monomials ) { 
	  QDPIO::cerr << "Index out of range. S_long_monomials["<<i<<"] has index " << S_long_monomials[i] << " but max index is " << num_ham_monomials-1 << endl;
	  QDP_abort(1);
	}
	
	// Check it hasn't been used already -- if it has it is either
	// a duplicate in long monomials, or has been used in the short
	// monomials 
	if( exclusivity_check[S_long_monomials[i]] == true ) { 
	  QDPIO::cerr << "Exclusivity check failed. Monomial " << i << " is listed more than once in list of long monomials or is duplicated between short and long monomials" << endl;
	  QDP_abort(1);
	}
	else { 
	  exclusivity_check[S_long_monomials[i]] = true;
	}
      }

      // Now go through exclusivity check array. If a member is false,
      // it means that that monomial term has been missed in both
      // the short and long monomial lists
      for(int i=0; i < num_ham_monomials; i++) { 
	if( exclusivity_check[i] == false ) { 
	  QDPIO::cerr << "Monomial " << i << " not used in either long or short monomial list " << endl;
	  QDP_abort(1);
	}
      }

      // At this point we have checked that
      //  -- the number of long and short monomials sums to the correct
      //     number of terms
      //  -- that all  the monomials have been marked one way or the other
      //  -- that there are no duplicates between the lists

      QDPIO::cout << "Long Monomials: " ;
      for(int i=0; i < S_long_monomials.size(); i++) { 
	QDPIO::cout << " " << S_long_monomials[i];
      }
      QDPIO::cout << endl;

      QDPIO::cout << "Short Monomials: " ;
      for(int i=0; i < S_short_monomials.size(); i++) { 
	QDPIO::cout << " " << S_short_monomials[i];
      }
      QDPIO::cout << endl;

      QDPIO::cout << " tau0=" << tau0;
      QDPIO::cout << " lambda = " << lambda;
      QDPIO::cout << " dt = " << getStepSize();
      QDPIO::cout << " n_steps= " << n_steps;
      QDPIO::cout << " n_short_steps= " << n_short_steps;


    }

    // Copy constructor
    LatColMatSexWeinMixedIntegrator(const LatColMatSexWeinMixedIntegrator& l) :
      n_steps(l.n_steps), lambda(l.lambda), tau0(l.tau0), H_MD(l.H_MD), n_short_steps(l.n_short_steps), S_short_monomials(l.S_short_monomials), S_long_monomials(l.S_long_monomials) {}

    // ! Destruction is automagic
    ~LatColMatSexWeinMixedIntegrator(void) {};

    //! Get at the MD Hamiltonian
    AbsHamiltonian<multi1d<LatticeColorMatrix>,
		   multi1d<LatticeColorMatrix> >& getHamiltonian(void) {
      return *H_MD;
    }

    //! Do a trajectory
    void operator()(AbsFieldState<multi1d<LatticeColorMatrix>,
		    multi1d<LatticeColorMatrix> >& s) {

      Real dt = getStepSize();      // Overall step size
      Real dtbyN = dt/n_short_steps;  // length of a short step
      Real lambda_dt = dt*lambda;
      Real middle_dt = (Real(1)-Real(2)*lambda)*dt;
      Real two_lambda_dt = Real(2)*lambda_dt;

      // First half step with long timescale monomials
	U_long(lambda_dt, s);

	// First full step with short timescale monomials
	//  n_short_steps of length dtbyN (reversibility 
	//  and further step subdivision is taken care of in 
	//  U_short
	U_short_1(dtbyN, n_short_steps,  s);


	U_long(middle_dt, s);

	// Do new short timescale step (n_short_step steps of dtbyN
	// (further subdivision is taken care of in U_short)
	U_short_2(dtbyN, n_short_steps, s);
      
      // Now the main loop 
      for(int i=0; i < n_steps-1; i++) { 
	U_long(two_lambda_dt, s);


	// First full step with short timescale monomials
	//  n_short_steps of length dtbyN (reversibility 
	//  and further step subdivision is taken care of in 
	//  U_short
	U_short_1(dtbyN, n_short_steps,  s);


	U_long(middle_dt, s);

	// Do new short timescale step (n_short_step steps of dtbyN
	// (further subdivision is taken care of in U_short)
	U_short_2(dtbyN, n_short_steps, s);
      }

      U_long(lambda_dt, s);      
    }

  protected:

    //! LeapP for just a selected list of monomials
    void leapP(const multi1d<int>& monomial_list,
	       const Real& dt, 
	       AbsFieldState<multi1d<LatticeColorMatrix>,
	       multi1d<LatticeColorMatrix> >& s) {

      LCMMDIntegratorSteps::leapP(monomial_list, 
				  dt, 
				  getHamiltonian(),
				  s);
    }

    //! Leap with Q (with all monomials)
    void leapQ(const Real& dt, 
	       AbsFieldState<multi1d<LatticeColorMatrix>,
	       multi1d<LatticeColorMatrix> >& s) {
      LCMMDIntegratorSteps::leapQ(dt,s);
    }



    //! Does [ U_short(dt/n) ]^n
    inline void U_short_1(const Real& dtbyN, 
			  int N,
			  AbsFieldState<multi1d<LatticeColorMatrix>,
			     multi1d<LatticeColorMatrix> >&  s) {

      Real lambda_dtbyN = dtbyN*lambda;
      Real half_dtbyN = Real(0.5)*dtbyN;
      Real one_minus_two_lambda_dtbyn=Real(0.5)*(Real(1)-Real(2)*lambda)*dtbyN;

      
      for(int i=0; i < N; i++) {
	leapP(S_short_monomials, lambda_dtbyN, s);
	leapQ(half_dtbyN, s);
  	leapP(S_short_monomials, one_minus_two_lambda_dtbyn, s);
      }

    }

    //! Does [ U_short(dt/n) ]^n
    inline void U_short_2(const Real& dtbyN, 
			  int N,
			  AbsFieldState<multi1d<LatticeColorMatrix>,
			     multi1d<LatticeColorMatrix> >&  s) {

      Real lambda_dtbyN = dtbyN*lambda;
      Real half_dtbyN = Real(0.5)*dtbyN;
      Real one_minus_two_lambda_dtbyn=Real(0.5)*(Real(1)-Real(2)*lambda)*dtbyN;

      

      for(int i=0; i < N; i++) {
	leapP(S_short_monomials, one_minus_two_lambda_dtbyn, s);
	leapQ(half_dtbyN, s);
	leapP(S_short_monomials, lambda_dtbyN, s);
      }
    }

    inline void U_long(const Real& dt, 
		AbsFieldState<multi1d<LatticeColorMatrix>,
		multi1d<LatticeColorMatrix> >& s) {
  
      
      leapP(S_long_monomials, dt,s);


    }


    //! Get the trajectory length
    const Real getTrajLength(void) const {
      return tau0;
    }

    //! Get the step size 
    const Real getStepSize(void) const {
      return tau0/Real(n_steps);
    }

  private:
    int n_steps;
    Real lambda;
    Real tau0;
    Handle< AbsHamiltonian<multi1d<LatticeColorMatrix>, 
			    multi1d<LatticeColorMatrix> > > H_MD;
    int n_short_steps;
    multi1d<int> S_short_monomials;
    multi1d<int> S_long_monomials;
    
    
  };

};


#endif

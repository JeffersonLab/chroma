// $Id: lcm_minimum_norm2_integrator_mts.h,v 2.1 2006-02-27 15:34:39 bjoo Exp $
/*! @file
 * @brief Second order minimal norm (2MN) integrator with multiple time scales
 *
 * 2MN intregator with multiple time scales for HMC
 */

#ifndef LCM_MINIMUM_NORM_MTS_H
#define LCM_MINIMUM_NORM_MTS_H


#include "chromabase.h"
#include "update/molecdyn/hamiltonian/abs_hamiltonian.h"
#include "update/molecdyn/integrator/abs_integrator.h"
#include "update/molecdyn/integrator/lcm_integrator_leaps.h"
#include <string>
#include <iostream>

namespace Chroma 
{

  /*! @ingroup integrator */
  namespace LatColMatMinimumNorm2IntegratorMtsEnv {
    extern const std::string name;
    extern const bool registered;
  };

  /*! @ingroup integrator */
  struct  LatColMatMinimumNorm2IntegratorMtsParams {
    LatColMatMinimumNorm2IntegratorMtsParams();
    LatColMatMinimumNorm2IntegratorMtsParams(XMLReader& xml, const std::string& path);

    int number_of_timescales; // Number of timescales
    Real tau; // trajectory length

    multi1d<int> n_steps_list;  // number of integration steps sorted by timescales
    multi1d<Real> lambda_list; // list of lambda parameters
    multi1d< multi1d<int> > monomial_list; // Indices of monomials sorted by timescales

  };

  /*! @ingroup integrator */
  void read(XMLReader& xml_in, 
	    const std::string& path,
	    LatColMatMinimumNorm2IntegratorMtsParams& p);

  /*! @ingroup integrator */
  void write(XMLWriter& xml_out,
	     const std::string& path, 
	     const LatColMatMinimumNorm2IntegratorMtsParams& p);

  //! MD integrator interface for PQP leapfrog
  /*! @ingroup integrator
   *  Specialised to multi1d<LatticeColorMatrix>
   */
  class LatColMatMinimumNorm2IntegratorMts 
    : public AbsMDIntegrator<multi1d<LatticeColorMatrix>,
    multi1d<LatticeColorMatrix> > {
    
    public:
    
    // Construct from params struct and Hamiltonian
    LatColMatMinimumNorm2IntegratorMts(const LatColMatMinimumNorm2IntegratorMtsParams& p,
				  Handle< AbsHamiltonian< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >& H_) : n_steps_list(p.n_steps_list), tau(p.tau), H_MD(H_), monomial_list(p.monomial_list), lambda_list(p.lambda_list), number_of_timescales(p.number_of_timescales) {
      
      // Check the number of timescales
      if(number_of_timescales != monomial_list.size()) {
	QDPIO::cerr << "Error: Wrong number of time scales specified for 2MN MTS integrator." << endl;
	QDPIO::cerr << "Integrator has    " << number_of_timescales << endl;
	QDPIO::cerr << "Monomial list has " << monomial_list.size() << endl;
	QDP_abort(1);
      }
      
      // Check for sanity
      int num_ham_monomials = getHamiltonian().numMonomials();

      // Sum up the number of monomials available
      int num_monomials = 0;
      for(int i = 0; i < p.monomial_list.size1(); i++) {
	num_monomials += p.monomial_list[i].size1();
      }

      // Check that the number of short timescale monomials + 
      //   the number of long_timescale monomials equals the total number of 
      //   monomials in the hamiltonian
      if( num_ham_monomials != num_monomials ) {
	QDPIO::cerr << "Error: Wrong number of monomials specified for 2MN MTS integrator." << endl;
	QDPIO::cerr << "Hamiltonian has " << num_ham_monomials << endl;
	QDPIO::cerr << "Integrator has  " << num_monomials << endl;
	QDP_abort(1);
      }

      // Check that a monomial has not been used twice and that 
      // no bogus indices are specified
      multi1d<bool> exclusivity_check(num_ham_monomials);
     
      // Initialise filter -- set all elements to false
      for(int i=0; i < num_ham_monomials; i++) { 
	exclusivity_check[i] = false;
      }

      // Go through list of monomials
      for(int i=0; i < monomial_list.size1(); i++) {
	for(int j = 0; j < monomial_list[i].size1(); j++) {
	  // Check for bogus index
	  if( monomial_list[i][j] < 0 || 
	      monomial_list[i][j] >= num_ham_monomials ) { 
	    
	    QDPIO::cerr << "Index out of range. monomial_list["<<i<<"]["<<j<<" has index " 
			<< monomial_list[i][j] << " but max index is " << num_ham_monomials-1 << endl;
	    QDP_abort(1);
	  }
	  
	  // If index is not bogus, check it has not been used already
	  // ie only insert it if exclusivity_check[i] is false
	  if( exclusivity_check[monomial_list[i][j]] == true ) { 
	    QDPIO::cerr << "Exclusivity check failed. Monomial " << i << " is listed more than once in list of monomials" << endl;
	    QDP_abort(1);
	  }
	  else { 
	    // Its not there yet, so insert it
	    exclusivity_check[monomial_list[i][j]] = true;
	  }
        }
      }

      // Now go through exclusivity check array. If a member is false,
      // it means that that monomial term has been missed in both
      // the short and long monomial lists
      for(int i=0; i < num_ham_monomials; i++) { 
	if( exclusivity_check[i] == false ) { 
	  QDPIO::cerr << "Monomial " << i << " not used in monomial list " << endl;
	  QDP_abort(1);
	}
      }

      // Check here also the lambda parameter!
      if(lambda_list.size() != number_of_timescales) {
	QDPIO::cerr << "Error! number of lambda parameter does not match the number of timescales" << endl;
	QDPIO::cerr << "number of timesscales " << number_of_timescales << endl;
	QDPIO::cerr << "number of lambda parameter " << lambda_list.size() << endl;
	  QDP_abort(1);
      }

      // At this point we have checked that
      //  -- the number of monomials sums to the correct
      //     number of terms
      //  -- that all the monomials have been marked one way or the other
      //  -- that there are no duplicates between the lists

    }

    // Copy constructor
    LatColMatMinimumNorm2IntegratorMts(const LatColMatMinimumNorm2IntegratorMts& l) :
      n_steps_list(l.n_steps_list), tau(l.tau), H_MD(l.H_MD), monomial_list(l.monomial_list), lambda_list(l.lambda_list), number_of_timescales(l.number_of_timescales)  {}

    // ! Destruction is automagic
    ~LatColMatMinimumNorm2IntegratorMts(void) {};

    //! Get at the MD Hamiltonian
    AbsHamiltonian<multi1d<LatticeColorMatrix>,
      multi1d<LatticeColorMatrix> >& getHamiltonian(void) {
      return *H_MD;
    }

    //! Do a trajectory
    void operator()(AbsFieldState<multi1d<LatticeColorMatrix>,
		    multi1d<LatticeColorMatrix> >& s) {
      int recursive_seed = number_of_timescales-1;

      // The half step at the beginning for all the momenta
      Real dtau  = tau/Real(n_steps_list[recursive_seed]);
      for(int i = number_of_timescales-1; i > 0; i--) {
	leapP(monomial_list[i], lambda_list[i]*dtau, s);
	dtau /= Real(n_steps_list[i-1]*2);
      }
      leapP(monomial_list[0], lambda_list[0]*dtau, s);

      // The rest we do recursively
      recursive_integrator(recursive_seed, tau, true, s);

      // The half step at the end for all the momenta
      dtau = tau/Real(n_steps_list[recursive_seed]);
      for(int i = number_of_timescales-1; i > 0; i--) {
	leapP(monomial_list[i], lambda_list[i]*dtau, s);
	dtau /= Real(n_steps_list[i-1]*2);
      }
      leapP(monomial_list[0], lambda_list[0]*dtau, s);
    }

    protected:

    void recursive_integrator(const int recursion_index, const Real tau0, const bool halfstep,
			      AbsFieldState<multi1d<LatticeColorMatrix>,
			      multi1d<LatticeColorMatrix> >& s) {
      Real dtau = tau0/Real(n_steps_list[recursion_index]);
      Real dtauby2 = dtau/Real(2);
      Real lambda = lambda_list[recursion_index];

      if(recursion_index == 0) {
	for(int i = 1; i < n_steps_list[recursion_index]; i++) {
	  leapQ(dtauby2, s);
	  leapP(monomial_list[recursion_index], 
		(1.-2.*lambda)*dtau, s);
	  leapQ(dtauby2, s);
	  leapP(monomial_list[recursion_index], 
		2*lambda*dtau, s);
	}
	leapQ(dtauby2, s);
	leapP(monomial_list[recursion_index], 
	      (1.-2.*lambda)*dtau, s);
	leapQ(dtauby2, s);
	if(!halfstep) {
	  leapP(monomial_list[recursion_index], 
		2*lambda*dtau, s);
	}
      }
      else {
	for(int i = 1; i < n_steps_list[recursion_index]; i++){
	  recursive_integrator(recursion_index-1, dtauby2, false, s);
	  leapP(monomial_list[recursion_index], 
		(1.-2.*lambda)*dtau, s);
	  recursive_integrator(recursion_index-1, dtauby2, false, s);
	  leapP(monomial_list[recursion_index], 
		2*lambda*dtau, s);
	}
	recursive_integrator(recursion_index-1, dtauby2, false, s);
	leapP(monomial_list[recursion_index], 
	      (1.-2.*lambda)*dtau, s);
	if(recursion_index == (number_of_timescales-1)) {
	  recursive_integrator(recursion_index-1, dtauby2, true, s);
	}
	else {
	  recursive_integrator(recursion_index-1, dtauby2, halfstep, s);
	}
	if(!halfstep) {
	  leapP(monomial_list[recursion_index],
		2*lambda*dtau, s);
	}
      }
    }

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




    //! Get the trajectory length
    const Real getTrajLength(void) const {
      return tau;
    }


    private:
    multi1d<int> n_steps_list;
    multi1d<Real> lambda_list;
    multi1d< multi1d<int> > monomial_list;
    Real tau;
    int number_of_timescales;
    Handle< AbsHamiltonian<multi1d<LatticeColorMatrix>, 
      multi1d<LatticeColorMatrix> > > H_MD;
    
    
  };

};


#endif
